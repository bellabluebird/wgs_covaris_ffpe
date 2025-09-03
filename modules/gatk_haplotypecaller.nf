// modules/gatk_haplotypecaller.nf

process GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: marked BAM file, index, and reference files
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_files
    
    // output: GVCF and VCF files
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcf
    tuple val(sample_id), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf_index
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}.vcf.gz.tbi"), emit: vcf_index
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    """
    # call variants with GATK HaplotypeCaller
    gatk HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -O ${sample_id}.g.vcf.gz \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus}
    
    # genotype GVCF to get regular VCF
    gatk GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V ${sample_id}.g.vcf.gz \\
        -O ${sample_id}.vcf.gz
    """
}