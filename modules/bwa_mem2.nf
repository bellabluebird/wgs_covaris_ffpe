// modules/bwa_mem2.nf

process BWA_MEM2 {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment", mode: params.publish_mode
    
    // conda option
    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.17'
    
    // input: paired trimmed reads from fastp
    input:
    tuple val(sample_id), path(reads)
    path reference_files
    
    // output: sorted BAM file and index
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai
    
    script:
    def ref_fasta = reference_files.find { it.toString().endsWith('.fasta') }
    """
    # align reads to reference genome and sort
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${ref_fasta} \\
        ${reads[0]} \\
        ${reads[1]} | \\
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${sample_id}.sorted.bam
    
    # index the BAM file
    samtools index ${sample_id}.sorted.bam
    """
}