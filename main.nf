#!/usr/bin/env nextflow

// basic qc pipeline
// BP MAKE alterations to handle 30x sequencing, this will only work on paired-end reads

// enable modular functions
nextflow.enable.dsl = 2

// process modules
include { FASTQC } from './modules/fastqc.nf'
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

// creating input channel from paired fastq files
// try multiple patterns to be more flexible with file naming
ch_input = Channel
    .fromFilePairs("${params.input_dir}/*_{R1,R2}.fastq.gz", size: 2)
    .ifEmpty { 
        log.error "Cannot find any paired FASTQ files in: ${params.input_dir}"
        log.error "Expected pattern: *_{R1,R2}.fastq.gz"
        log.error "Your files should be named like: ERR008539_R1.fastq.gz, ERR008539_R2.fastq.gz"
        log.error "Full pattern being searched: ${params.input_dir}/*_{R1,R2}.fastq.gz"
        error "No paired FASTQ files found"
    }

// main workflow
workflow {
    // debug: show what files were found
    ch_input.view { sample_id, files -> "Found sample: ${sample_id} with files: ${files}" }
    
    // raw fastqc
    FASTQC(ch_input)
    
    // preprocessing with fastp
    FASTP(ch_input)
    
    // post-trim fastqc - create alias to make FASTQC reusable
    FASTQC_TRIMMED = FASTQC
    FASTQC_TRIMMED(FASTP.out.reads)
    
    // collect all reports for multiqc
    multiqc_input = FASTQC.out.zip
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(FASTP.out.json)
        .collect()
    
    // multiqc report
    MULTIQC(multiqc_input)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Results directory: ${params.outdir}"
}