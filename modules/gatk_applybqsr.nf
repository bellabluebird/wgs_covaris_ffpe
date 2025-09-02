// modules/gatk_applybqsr.nf

process GATK_APPLYBQSR {
    tag "$sample_id"
    publishDir "${params.outdir}/bqsr", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::gatk4=4.4.0.0 bioconda::samtools=1.17'
    
    // input: marked BAM file, index, reference files, and recalibration table
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_files
    tuple val(sample_id), path(recal_table)
    
    // output: BQSR-recalibrated BAM file and index
    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.recal.bam.bai"), emit: bai
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    """
    # apply base quality score recalibration
    gatk ApplyBQSR \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${sample_id}.recal.bam
    
    # index the recalibrated BAM file
    samtools index ${sample_id}.recal.bam
    """
}