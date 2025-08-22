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
params.input_dir = ""
params.outdir = "./results"
params.publish_mode = 'copy'

// input validation
if (!params.input_dir) {
    error "Please provide input directory with --input_dir"
}

// creating input channel from paired fastq files
ch_input = Channel
    .fromFilePairs("${params.input_dir}/*_{R1,R2,1,2}*.{fastq,fq}.gz", size: 2)
    .ifEmpty { error "Cannot find any paired FASTQ files in: ${params.input_dir}" }

// main workflow
workflow {
    // raw fastqc
    FASTQC(ch_input)
    
    // preprocessing with fastp
    FASTP(ch_input)
    
    // post-trim fastqc
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