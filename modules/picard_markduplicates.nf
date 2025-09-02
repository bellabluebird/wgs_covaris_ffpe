// modules/picard_markduplicates.nf

process PICARD_MARKDUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/markduplicates", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::picard=2.27.4'
    
    // input: sorted BAM file
    input:
    tuple val(sample_id), path(bam)
    
    // output: marked BAM file, index, and metrics
    output:
    tuple val(sample_id), path("${sample_id}.marked.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.marked.bam.bai"), emit: bai
    tuple val(sample_id), path("${sample_id}.marked_dup_metrics.txt"), emit: metrics
    
    script:
    """
    # mark duplicates with Picard
    picard MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}.marked.bam \\
        METRICS_FILE=${sample_id}.marked_dup_metrics.txt \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=SILENT \\
        TMP_DIR=\$PWD/tmp
    
    # create tmp directory if it doesn't exist
    mkdir -p tmp
    """
}