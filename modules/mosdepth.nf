// modules/mosdepth.nf

process MOSDEPTH {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::mosdepth=0.3.8'
    
    // input: marked BAM file and index
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    // output: coverage files for analysis and multiqc
    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}.per-base.bed.gz"), emit: per_base, optional: true
    
    script:
    """
    # verify BAM file is indexed
    if [[ ! -f "${bai}" ]]; then
        echo "ERROR: BAM index file not found: ${bai}"
        exit 1
    fi
    
    # run mosdepth for coverage analysis with larger windows for efficiency
    mosdepth \\
        --threads ${task.cpus} \\
        --by 5000 \\
        ${sample_id} \\
        ${bam}
    """
}