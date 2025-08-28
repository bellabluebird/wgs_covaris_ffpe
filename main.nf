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

// test different approaches to find files
log.info "Testing different file discovery methods:"

// Method 1: Try the current pattern
log.info "Method 1: Using fromFilePairs with pattern *_R{1,2}.fastq.gz"
def pattern1 = "${params.input_dir}/*_R{1,2}.fastq.gz"
log.info "Full pattern: ${pattern1}"

try {
    def test_ch1 = Channel.fromPath(pattern1)
    test_ch1.view { "Method 1 found file: ${it}" }
    log.info "Method 1: Pattern matching attempted"
} catch (Exception e) {
    log.error "Method 1 failed with error: ${e.getMessage()}"
}

// Method 2: Try alternative pattern
log.info "Method 2: Using fromPath with pattern *_R*.fastq.gz"
def pattern2 = "${params.input_dir}/*_R*.fastq.gz"
log.info "Full pattern: ${pattern2}"

try {
    def test_ch2 = Channel.fromPath(pattern2)
    test_ch2.view { "Method 2 found file: ${it}" }
} catch (Exception e) {
    log.error "Method 2 failed with error: ${e.getMessage()}"
}

// Method 3: Try listing everything in the directory
log.info "Method 3: Listing all files in directory"
def pattern3 = "${params.input_dir}/*"
try {
    def test_ch3 = Channel.fromPath(pattern3)
    test_ch3.view { "Method 3 found file: ${it}" }
} catch (Exception e) {
    log.error "Method 3 failed with error: ${e.getMessage()}"
}

log.info "=== END FILE DISCOVERY DEBUGGING ==="

// Create S3 file pairs manually (more reliable for S3)
log.info "Creating S3 file channels manually for ERR008539"
ch_input = Channel.of([
    'ERR008539', 
    ["${params.input_dir}/ERR008539_R1.fastq.gz", "${params.input_dir}/ERR008539_R2.fastq.gz"]
])

// main workflow
workflow {
    // First, run S3 debugging process
    DEBUG_S3_ACCESS() | view { "S3 DEBUG OUTPUT: $it" }
    
    // debug: show what files were found
    ch_input.view { sample_id, files -> "Found sample: ${sample_id} with files: ${files}" }
    
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