#!/usr/bin/env nextflow

// basic qc pipeline
// BP MAKE alterations to handle 30x sequencing, this will only work on paired-end reads

// enable modular functions
nextflow.enable.dsl = 2

// process modules
include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'
include { FASTP } from './modules/fastp.nf'
include { BWA_MEM2_INDEX } from './modules/bwa_mem2_index.nf'
include { BWA_MEM2 } from './modules/bwa_mem2.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
// new modules i'm working on implementing - uncomment when ready
// include { QUALIMAP } from './modules/qualimap.nf'
// include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates.nf'
// include { PICARD_COLLECTINSERTSIZEMETRICS } from './modules/picard_insert_size.nf'
// include { MOSDEPTH } from './modules/mosdepth.nf'
// include { GATK_BASERECALIBRATOR } from './modules/gatk_baserecalibrator.nf'
// include { GATK_APPLYBQSR } from './modules/gatk_applybqsr.nf'
// include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
// include { BCFTOOLS_STATS } from './modules/bcftools_stats.nf'
include { MULTIQC } from './modules/multiqc.nf'

// parameters
// input_dir comes from command line --input_dir parameter
params.outdir = "./results"
params.publish_mode = 'copy'
// PHASE 5: BQSR parameters (uncomment when implementing BQSR)
// params.known_sites = null  // optional known sites for BQSR (e.g., dbSNP, 1000G, Mills)

// input validation with debugging
log.info "DEBUG: Received input_dir parameter: '${params.input_dir}'"
log.info "DEBUG: Parameter type: ${params.input_dir?.getClass()}"
log.info "DEBUG: Parameter length: ${params.input_dir?.length()}"

if (!params.input_dir) {
    error "Please provide input directory with --input_dir (received: '${params.input_dir}')"
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

// PHASE 5: Known sites channel for BQSR (uncomment when implementing BQSR)
// ch_known_sites = params.known_sites ? 
//     Channel.fromPath(params.known_sites.split(',').collect { it.trim() })
//         .collect() : 
//     Channel.empty()

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
    
    // PHASE 2: Mark duplicates and get metrics (uncomment when implementing)
    // markduplicates_results = PICARD_MARKDUPLICATES(bwa_results.bam)
    // ch_marked_bam = markduplicates_results.bam.join(markduplicates_results.bai)
    
    // PHASE 3: Insert size metrics (uncomment when implementing)  
    // insert_size_metrics = PICARD_COLLECTINSERTSIZEMETRICS(markduplicates_results.bam)
    
    // PHASE 1: Quality metrics with Qualimap (uncomment when implementing)
    // qualimap_results = QUALIMAP(ch_marked_bam)
    
    // PHASE 4: Coverage analysis (uncomment when implementing)
    // coverage_results = MOSDEPTH(ch_marked_bam)
    
    // PHASE 5: Base quality score recalibration (uncomment when implementing)
    // if (params.known_sites) {
    //     bqsr_table = GATK_BASERECALIBRATOR(ch_marked_bam, ch_reference_indexed, ch_known_sites)
    //     bqsr_results = GATK_APPLYBQSR(ch_marked_bam, ch_reference_indexed, bqsr_table.recal_table)
    //     ch_final_bam = bqsr_results.bam.join(bqsr_results.bai)
    // } else {
    //     ch_final_bam = ch_marked_bam
    // }
    
    // PHASE 6: Variant calling (uncomment when implementing)
    // variant_results = GATK_HAPLOTYPECALLER(ch_final_bam, ch_reference_indexed)
    // variant_stats = BCFTOOLS_STATS(variant_results.vcf)
    
    // collect all reports for multiqc
    multiqc_input = fastqc_raw.zip.map { id, file -> file }
    .mix(fastqc_trimmed.zip.map { id, file -> file })
    .mix(fastp_results.json.map { id, file -> file })
    .mix(samtools_stats.stats.map { id, file -> file })
    // PHASE 2: .mix(markduplicates_results.metrics.map { id, file -> file })
    // PHASE 3: .mix(insert_size_metrics.metrics.map { id, file -> file })
    // PHASE 1: .mix(qualimap_results.genome_results.map { id, file -> file })
    // PHASE 4: .mix(coverage_results.summary.map { id, file -> file })
    // PHASE 6: .mix(variant_stats.stats.map { id, file -> file })
    .collect()

    // multiqc report
    MULTIQC(multiqc_input)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Results directory: ${params.outdir}"
}