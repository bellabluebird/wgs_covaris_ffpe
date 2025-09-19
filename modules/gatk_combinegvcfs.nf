// modules/gatk_combinegvcfs.nf

process GATK_COMBINEGVCFS {
    tag "combine_gvcfs"
    publishDir "${params.outdir}/gvcf", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: collection of GVCF files, their indices, and reference files
    input:
    path gvcfs
    path gvcf_indices
    path reference_files
    
    // output: combined GVCF file for joint genotyping
    output:
    path "cohort_combined.g.vcf.gz", emit: gvcf
    path "cohort_combined.g.vcf.gz.tbi", emit: gvcf_index
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    def gvcf_inputs = gvcfs.collect { "--variant ${it}" }.join(' ')
    """
    # combine individual GVCFs into single cohort GVCF
    gatk CombineGVCFs \\
        -R ${ref_fasta} \\
        ${gvcf_inputs} \\
        -O cohort_combined.g.vcf.gz
    """
}