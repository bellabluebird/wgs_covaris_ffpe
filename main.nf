#!/usr/bin/env nextflow

// basic qc pipeline
// BP MAKE alterations to handle 30x sequencing, this will only work on paired-end reads

// enable modular functions
nextflow.enable.dsl = 2

// process modules
include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'
include { FASTP } from './modules/fastp.nf'
include { MULTIQC } from './modules/multiqc.nf'

// parameters
// input_dir comes from command line --input_dir parameter
params.outdir = "./results"
params.publish_mode = 'copy'

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

// Using the current pattern
log.info "Using fromFilePairs with pattern *_R{1,2}.fastq.gz"
def pattern1 = "${params.input_dir}/*_R{1,2}.fastq.gz"
log.info "Full pattern: ${pattern1}"

// Create S3 file pairs manually
log.info "Creating S3 file channels manually for ERR008539"
ch_input = Channel.of([
    'ERR008539', 
    ["${params.input_dir}/ERR008539_R1.fastq.gz", "${params.input_dir}/ERR008539_R2.fastq.gz"]
])

// main workflow
workflow {
    // raw fastqc
    fastqc_raw = FASTQC(ch_input)
    
    // preprocessing with fastp
    fastp_results = FASTP(ch_input)
    
    // post-trim fastqc on cleaned reads
    fastqc_trimmed = FASTQC_TRIMMED(fastp_results.reads)
    
    // collect all reports for multiqc
    multiqc_input = fastqc_raw.zip
        .mix(fastqc_trimmed.zip)
        .mix(fastp_results.json)
        .collect()
    
    // multiqc report
    MULTIQC(multiqc_input)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Results directory: ${params.outdir}"
}