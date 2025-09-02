// modules/gatk_baserecalibrator.nf

process GATK_BASERECALIBRATOR {
    tag "$sample_id"
    publishDir "${params.outdir}/bqsr", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: marked BAM file, index, reference files, and known sites
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_files
    path known_sites
    
    // output: recalibration table
    output:
    tuple val(sample_id), path("${sample_id}.recal_data.table"), emit: recal_table
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    def known_sites_args = known_sites ? known_sites.collect { "--known-sites ${it}" }.join(' ') : ""
    """
    # generate base recalibration table
    gatk BaseRecalibrator \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        ${known_sites_args} \\
        -O ${sample_id}.recal_data.table
    """
}