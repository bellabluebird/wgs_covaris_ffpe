// modules/fastp.nf

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: params.publish_mode
    
    // input: paired fastq files
    input:
    tuple val(sample_id), path(reads)
    
    // output: sample name + trimmed files
    // output: sample name + fastp json report
    // output: sample + fastp html report
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    
    script:
    """
    # run fastp with paired-end reads
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}_R1_trimmed.fastq.gz \\
        --out2 ${sample_id}_R2_trimmed.fastq.gz \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe
    """
}