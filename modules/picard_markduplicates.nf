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
    # create tmp directory if it doesn't exist
    mkdir -p tmp
    
    # mark duplicates with Picard (without indexing)
    java -jar /usr/local/share/picard-*/picard.jar MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}.marked.bam \\
        METRICS_FILE=${sample_id}.marked_dup_metrics.txt \\
        VALIDATION_STRINGENCY=SILENT \\
        TMP_DIR=\$PWD/tmp
    
    # create index with samtools (more reliable in mulled container)
    samtools index ${sample_id}.marked.bam
    
    # verify that both output files were created
    if [[ ! -f "${sample_id}.marked.bam" ]]; then
        echo "ERROR: BAM file ${sample_id}.marked.bam was not created"
        exit 1
    fi
    
    if [[ ! -f "${sample_id}.marked.bam.bai" ]]; then
        echo "ERROR: Index file ${sample_id}.marked.bam.bai was not created"
        exit 1
    fi
    
    echo "SUCCESS: Both BAM and index files created successfully"
    ls -la ${sample_id}.marked.bam*
    """
}