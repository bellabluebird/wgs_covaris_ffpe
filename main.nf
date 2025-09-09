#!/usr/bin/env nextflow

// WGS Quality Control & Analysis Pipeline with BQSR
// Requires: --input_dir, --reference, --known_sites

// enable modular functions
nextflow.enable.dsl = 2

// process modules
include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'
include { FASTP } from './modules/fastp.nf'
//include { BWA_MEM2_INDEX } from './modules/bwa_mem2_index.nf'
include { BWA } from './modules/bwa.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKED } from './modules/samtools_index.nf'
include { QUALIMAP } from './modules/qualimap.nf'
include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates.nf'
include { PICARD_COLLECTINSERTSIZEMETRICS } from './modules/picard_insert_size.nf'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { GATK_BASERECALIBRATOR } from './modules/gatk_baserecalibrator.nf'
include { GATK_APPLYBQSR } from './modules/gatk_applybqsr.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
include { BCFTOOLS_STATS } from './modules/bcftools_stats.nf'
include { MULTIQC } from './modules/multiqc.nf'

// parameters - all data-specific parameters must be provided explicitly
params.outdir = "./results"
params.publish_mode = 'copy'

// required parameters (no defaults - must be specified via CLI or profile)
params.input_dir = null
params.reference = null
params.known_sites = null

// input validation
log.info "Input directory: ${params.input_dir}"
log.info "Reference genome: ${params.reference}"
log.info "Known sites: ${params.known_sites}"

// create input channel from FASTQ files in input directory
ch_input = Channel.fromFilePairs("${params.input_dir}/*_{R1,R2,1,2}.{fastq,fq}{,.gz}", checkIfExists: false)
    .ifEmpty { error "No paired FASTQ files found in ${params.input_dir}. Expected pattern: *_{R1,R2,1,2}.{fastq,fq}{,.gz}" }

// create reference channel from explicit parameter
ch_reference_fasta = Channel.fromPath(params.reference)
    .ifEmpty { error "Reference genome not found at ${params.reference}" }

// create known sites channel for BQSR (use parameter or default compatible sites)
ch_known_sites = Channel.fromPath((params.known_sites ?: params.default_known_sites).split(',').collect { it.trim() })
    .ifEmpty { error "No known sites files found at specified paths" }
    .collect() 

// main workflow
workflow {
    // raw fastqc
    fastqc_raw = FASTQC(ch_input)
    
    // preprocessing with fastp
    fastp_results = FASTP(ch_input)
    
    // post-trim fastqc on cleaned reads
    fastqc_trimmed = FASTQC_TRIMMED(fastp_results.reads)
    
    // use compatible references from config parameters (chr1, chr2 format)
    ch_reference_indexed = Channel.fromPath(params.reference_indexed_path ?: "s3://bp-wgs-covaris-input-data/reference/bwa/*")
        .collect()
    
    // create separate channel with GATK-required reference files (FASTA, .fai, .dict)
    ch_reference_fasta_gatk = Channel.fromPath(params.reference_gatk_path ?: "s3://bp-wgs-covaris-input-data/reference/bwa/GCA_000001405.15_GRCh38_genomic.{fasta,fasta.fai,dict}")
        .collect()
    
    // alignment to reference genome
    bwa_results = BWA(fastp_results.reads, ch_reference_indexed)
    
    // alignment statistics
    samtools_stats = SAMTOOLS_STATS(bwa_results.bam)
    
    // index aligned BAM files
    samtools_index = SAMTOOLS_INDEX(bwa_results.bam)
    
    // mark duplicates (using indexed BAM)
    markduplicates_results = PICARD_MARKDUPLICATES(samtools_index.bam_bai)
    
    // index the marked BAM files
    samtools_index_marked = SAMTOOLS_INDEX_MARKED(markduplicates_results.bam)
    
    // insert size metrics (on duplicate-marked BAM)
    insert_size_metrics = PICARD_COLLECTINSERTSIZEMETRICS(markduplicates_results.bam)
    
    // coverage analysis with mosdepth (on duplicate-marked indexed BAM)
    coverage_results = MOSDEPTH(samtools_index_marked.bam_bai)
    
    // base quality score recalibration with GATK
    bqsr_table = GATK_BASERECALIBRATOR(samtools_index_marked.bam_bai, ch_reference_fasta_gatk, ch_known_sites)
    
    // apply BQSR to get recalibrated BAM
    bqsr_results = GATK_APPLYBQSR(samtools_index_marked.bam_bai, ch_reference_fasta_gatk, bqsr_table.recal_table)
    
    // quality metrics with Qualimap (on recalibrated BAM)
    qualimap_results = QUALIMAP(bqsr_results.bam.join(bqsr_results.bai))

    // collect all reports for multiqc
    multiqc_input = fastqc_raw.zip.map { id, file -> file }
    .mix(fastqc_trimmed.zip.map { id, file -> file })
    .mix(fastp_results.json.map { id, file -> file })
    .mix(samtools_stats.stats.map { id, file -> file })
    .mix(markduplicates_results.metrics.map { id, file -> file })
    .mix(insert_size_metrics.metrics.map { id, file -> file })
    .mix(coverage_results.summary.map { id, file -> file })
    .mix(qualimap_results.genome_results.map { id, file -> file })
    .collect()

    // multiqc report
    MULTIQC(multiqc_input)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Results directory: ${params.outdir}"
}