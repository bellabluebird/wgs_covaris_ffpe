// modules/fastqc.nf

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::fastqc=0.12.1'
    
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