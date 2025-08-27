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
// try multiple patterns to be more flexible with file naming
ch_input = Channel
    .fromFilePairs([
        "${params.input_dir}/*_{R1,R2}*.{fastq,fq}.gz",
        "${params.input_dir}/*_{1,2}*.{fastq,fq}.gz",
        "${params.input_dir}/*{R1,R2}*.{fastq,fq}.gz",
        "${params.input_dir}/*{1,2}.{fastq,fq}.gz"
    ], size: 2)
    .ifEmpty { 
        log.error "Cannot find any paired FASTQ files in: ${params.input_dir}"
        log.error "Tried patterns:"
        log.error "  - *_{R1,R2}*.{fastq,fq}.gz"
        log.error "  - *_{1,2}*.{fastq,fq}.gz" 
        log.error "  - *{R1,R2}*.{fastq,fq}.gz"
        log.error "  - *{1,2}.{fastq,fq}.gz"
        log.error "Make sure your files follow one of these naming conventions"
        error "No paired FASTQ files found - see patterns above"
    }

// main workflow
workflow {
    // debug: show what files were found
    ch_input.view { sample_id, files -> "Found sample: ${sample_id} with files: ${files}" }
    
    // raw fastqc
    FASTQC(ch_input)
    
    // preprocessing with fastp
    FASTP(ch_input)
    
    // post-trim fastqc
    FASTQC(FASTP.out.reads)
    
    // collect all reports for multiqc
    multiqc_input = FASTQC.out.zip
        .mix(FASTQC.out.zip)
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