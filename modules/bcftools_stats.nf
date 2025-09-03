// modules/bcftools_stats.nf

process BCFTOOLS_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_stats", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::bcftools=1.17'
    
    // input: VCF file
    input:
    tuple val(sample_id), path(vcf)
    
    // output: variant statistics for multiqc
    output:
    tuple val(sample_id), path("${sample_id}.bcftools_stats.txt"), emit: stats
    
    script:
    """
    # generate variant statistics with bcftools
    bcftools stats ${vcf} > ${sample_id}.bcftools_stats.txt
    """
}