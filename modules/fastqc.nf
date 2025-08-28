// modules/fastqc.nf

process FASTQC {
    // label for job while running
    tag "$sample_id"
    // 
    publishDir "${params.outdir}/fastqc", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::fastqc=0.12.1'
    // docker profile option; these are published biocontainers
    container 'biocontainers/fastqc:v0.11.9_cv8'
    
    // input: paired reads from fastp
    input:
    tuple val(sample_id), path(reads)
    
    // output: fastqc html 
    // output: zip files available to multiqc
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip
    
    script:
    """
    # run fastqc on the reads
    fastqc --threads ${task.cpus} ${reads}
    """
}