// modules/samtools_index.nf
process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/markduplicates", mode: params.publish_mode
    
    conda 'bioconda::samtools=1.17'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${bam}.bai"), emit: bai
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam_bai

    script:
    """
    samtools index ${bam}
    """
}