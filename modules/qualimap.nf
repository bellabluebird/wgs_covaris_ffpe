// modules/qualimap.nf

process QUALIMAP {
    tag "$sample_id"
    publishDir "${params.outdir}/qualimap", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::qualimap=2.2.2'
    
    // input: sorted BAM file and index
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    // output: qualimap results for multiqc
    output:
    tuple val(sample_id), path("${sample_id}_qualimap"), emit: results
    tuple val(sample_id), path("${sample_id}_qualimap/genome_results.txt"), emit: genome_results
    
    script:
    """
    # run qualimap bamqc for genome analysis
    qualimap bamqc \\
        -bam ${bam} \\
        -outdir ${sample_id}_qualimap \\
        -outformat HTML \\
        -nt ${task.cpus} \\
        --java-mem-size=${task.memory.toGiga()}G
    """
}