// modules/picard_insert_size.nf

process PICARD_COLLECTINSERTSIZEMETRICS {
    tag "$sample_id"
    publishDir "${params.outdir}/insert_size", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::picard=2.27.4'
    
    // input: marked BAM file
    input:
    tuple val(sample_id), path(bam)
    
    // output: insert size metrics and histogram
    output:
    tuple val(sample_id), path("${sample_id}.insert_size_metrics.txt"), emit: metrics
    tuple val(sample_id), path("${sample_id}.insert_size_histogram.pdf"), emit: histogram
    
    script:
    """
    # collect insert size metrics with Picard
    picard CollectInsertSizeMetrics \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}.insert_size_metrics.txt \\
        HISTOGRAM_FILE=${sample_id}.insert_size_histogram.pdf \\
        VALIDATION_STRINGENCY=SILENT \\
        TMP_DIR=\$PWD/tmp
    
    # create tmp directory if it doesn't exist
    mkdir -p tmp
    """
}