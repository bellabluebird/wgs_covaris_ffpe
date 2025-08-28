// modules/multiqc.nf

process MULTIQC {
    // where to publish results
    // mode: set as copy
    publishDir "${params.outdir}/multiqc", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::multiqc=1.19'
    // docker profile option; these are published biocontainers
    container 'staphb/multiqc:1.19'
    
    // path to file inputs
    input:
    path(reports)
    
    // output: multiqc report html
    // output: raw data directory w csv, json etc
    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
    
    // run multiqc in the current working directory
    script:
    """
    multiqc .
    """
}