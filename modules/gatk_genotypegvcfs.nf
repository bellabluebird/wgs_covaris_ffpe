// modules/gatk_genotypegvcfs.nf

process GATK_GENOTYPEGVCFS {
    tag "joint_genotyping"
    publishDir "${params.outdir}/variants", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: collection of GVCF files and reference files
    input:
    path gvcfs
    path reference_files
    
    // output: joint-called VCF file
    output:
    path "joint_genotyped.vcf.gz", emit: vcf
    path "joint_genotyped.vcf.gz.tbi", emit: vcf_index
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    def gvcf_inputs = gvcfs.collect { "--variant ${it}" }.join(' ')
    """
    # joint genotyping with GATK GenotypeGVCFs
    gatk GenotypeGVCFs \\
        -R ${ref_fasta} \\
        ${gvcf_inputs} \\
        -O joint_genotyped.vcf.gz
    """
}