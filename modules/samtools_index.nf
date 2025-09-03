// modules/samtools_index.nf
process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/markduplicates", mode: params.publish_mode
    
    container = 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.marked.bam.bai"), emit: bai

    script:
    """
    samtools index ${bam}
    """
}