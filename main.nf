#!/usr/bin/env nextflow

// WGS Quality Control & Analysis Pipeline with BQSR
// Requires: --input_dir, --known_sites
// Example: nextflow run main.nf --input_dir "s3://bucket/samples" --known_sites "s3://bucket/dbsnp.vcf.gz,s3://bucket/mills.vcf.gz"

// enable modular functions
nextflow.enable.dsl = 2

// process modules
include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'
include { FASTP } from './modules/fastp.nf'
include { BWA_MEM2_INDEX } from './modules/bwa_mem2_index.nf'
include { BWA_MEM2 } from './modules/bwa_mem2.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKED } from './modules/samtools_index.nf'
// new modules i'm working on implementing - uncomment when ready
include { QUALIMAP } from './modules/qualimap.nf'
include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates.nf'
include { PICARD_COLLECTINSERTSIZEMETRICS } from './modules/picard_insert_size.nf'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { GATK_BASERECALIBRATOR } from './modules/gatk_baserecalibrator.nf'
include { GATK_APPLYBQSR } from './modules/gatk_applybqsr.nf'
// include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
// include { BCFTOOLS_STATS } from './modules/bcftools_stats.nf'
include { MULTIQC } from './modules/multiqc.nf'

// parameters
// input_dir comes from command line --input_dir parameter
params.outdir = "./results"
params.publish_mode = 'copy'
// BQSR known sites parameter (required - comma-separated list of VCF files)
// Example: --known_sites "s3://bucket/dbsnp_138.hg38.vcf.gz,s3://bucket/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.known_sites = null

// input validation with debugging
log.info "DEBUG: Received input_dir parameter: '${params.input_dir}'"
log.info "DEBUG: Parameter type: ${params.input_dir?.getClass()}"
log.info "DEBUG: Parameter length: ${params.input_dir?.length()}"

if (!params.input_dir) {
    error "Please provide input directory with --input_dir (received: '${params.input_dir}')"
}

if (!params.known_sites) {
    error "Please provide known sites VCF files with --known_sites (required for BQSR)"
}

// creating input channel from paired fastq files with extensive debugging
log.info "=== FILE DISCOVERY DEBUGGING ==="
log.info "Input directory: ${params.input_dir}"
log.info "AWS region from config: ${params.aws_region ?: 'not set'}"
log.info "Current working directory: ${workflow.launchDir}"
log.info "Work directory: ${workflow.workDir}"

// using the current pattern
log.info "Using fromFilePairs with pattern *_R{1,2}.fastq.gz"
def pattern1 = "${params.input_dir}/*_R{1,2}.fastq.gz"
log.info "Full pattern: ${pattern1}"

// create S3 file pairs manually
log.info "Creating S3 file channels manually for ERR008539"
ch_input = Channel.of([
    'ERR008539', 
    ["${params.input_dir}/ERR008539_R1.fastq.gz", "${params.input_dir}/ERR008539_R2.fastq.gz"]
])

// create reference channel from input directory (just FASTA for indexing)
ch_reference_fasta = Channel.fromPath("${params.input_dir}/*.fasta")
    .ifEmpty { error "No reference genome (.fasta) found in ${params.input_dir}" }
    .first()

// create known sites channel for BQSR
ch_known_sites = Channel.fromPath(params.known_sites.split(',').collect { it.trim() })
    .ifEmpty { error "No known sites VCF files found at specified paths" }
    .collect() 

// main workflow
workflow {
    // raw fastqc
    fastqc_raw = FASTQC(ch_input)
    
    // preprocessing with fastp
    fastp_results = FASTP(ch_input)
    
    // post-trim fastqc on cleaned reads
    fastqc_trimmed = FASTQC_TRIMMED(fastp_results.reads)
    
    // create BWA-MEM2 index from reference
    bwa_index = BWA_MEM2_INDEX(ch_reference_fasta)
    
    // combine indexed reference files for alignment
    ch_reference_indexed = bwa_index.fasta.mix(bwa_index.index).collect()
    
    // alignment to reference genome
    bwa_results = BWA_MEM2(fastp_results.reads, ch_reference_indexed)
    
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
    bqsr_table = GATK_BASERECALIBRATOR(samtools_index_marked.bam_bai, ch_reference_indexed, ch_known_sites)
    
    // apply BQSR to get recalibrated BAM
    bqsr_results = GATK_APPLYBQSR(samtools_index_marked.bam_bai, ch_reference_indexed, bqsr_table.recal_table)
    
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