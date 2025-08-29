// modules/samtools_stats.nf

process SAMTOOLS_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment_stats", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::samtools=1.17'
    
    // input: sorted BAM file
    input:
    tuple val(sample_id), path(bam)
    
    // output: alignment statistics for multiqc
    output:
    tuple val(sample_id), path("${sample_id}.samtools_stats.txt"), emit: stats
    
    script:
    """
    # generate alignment statistics
    samtools stats ${bam} > ${sample_id}.samtools_stats.txt
    """
}