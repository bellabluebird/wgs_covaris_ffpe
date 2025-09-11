// modules/gatk_genotypegvcfs.nf

process GATK_GENOTYPEGVCFS {
    tag "joint_genotyping"
    publishDir "${params.outdir}/variants", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: combined GVCF file and reference files
    input:
    path combined_gvcf
    path gvcf_index
    path reference_files
    
    // output: joint-called VCF file
    output:
    path "joint_genotyped.vcf.gz", emit: vcf
    path "joint_genotyped.vcf.gz.tbi", emit: vcf_index
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    """
    # joint genotyping with GATK GenotypeGVCFs using combined GVCF
    gatk GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V ${combined_gvcf} \\
        -O joint_genotyped.vcf.gz
    """
}