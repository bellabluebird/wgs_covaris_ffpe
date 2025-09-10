// modules/gatk_variantfiltration.nf

process GATK_VARIANTFILTRATION {
    tag "variant_filtering"
    publishDir "${params.outdir}/variants", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0'
    
    // input: raw VCF file and reference files
    input:
    path vcf
    path vcf_index
    path reference_files
    
    // output: filtered VCF file
    output:
    path "filtered_variants.vcf.gz", emit: vcf
    path "filtered_variants.vcf.gz.tbi", emit: vcf_index
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    """
    # apply hard filtering to variants using GATK best practices
    gatk VariantFiltration \\
        -R ${ref_fasta} \\
        -V ${vcf} \\
        -O filtered_variants.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
        --filter-expression "FS > 60.0" --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
    """
}