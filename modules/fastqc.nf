#!/usr/bin/env nextflow

process FASTQC {
    // quality control analysis on individual fastq files
    // label 'process_low' = low resource requirements from config
    // container has fastqc tool pre-installed
    label 'process_low'
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    // from transposed channel: ↓ ["sample1", sample1_R1.fastq]

    output:
    tuple val(sample_id), path("*.html"), emit: html    // visual reports
    tuple val(sample_id), path("*.zip"), emit: zip      // data for multiqc
    // maintains sample_id for tracking through pipeline

    shell:
    """
    fastqc -t $task.cpus $reads
    """
    // runs fastqc with allocated cpu threads
    // $task.cpus gets cpu allocation from config
    // generates html report + zip file with metrics
}