// modules/picard_markduplicates.nf
process PICARD_MARKDUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/markduplicates", mode: params.publish_mode
    
    container 'quay.io/biocontainers/picard:2.27.4--hdfd78af_0'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.marked.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.marked_dup_metrics.txt"), emit: metrics

    script:
    """
    mkdir -p tmp
    picard MarkDuplicates \
        INPUT=${bam} \
        OUTPUT=${sample_id}.marked.bam \
        METRICS_FILE=${sample_id}.marked_dup_metrics.txt \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$PWD/tmp
    """
}
