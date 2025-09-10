// modules/gatk_haplotypecaller.nf

process GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "${params.outdir}/gvcf", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: recalibrated BAM file, index, and reference files
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_files
    
    // output: GVCF files only (for joint genotyping)
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcf
    tuple val(sample_id), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf_index
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    """
    # call variants with GATK HaplotypeCaller in GVCF mode
    gatk HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -O ${sample_id}.g.vcf.gz \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus}
    """
}