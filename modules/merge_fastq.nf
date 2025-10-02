// modules/merge_fastq.nf

process MERGE_FASTQ {
    tag "$sample_id"
    publishDir "${params.outdir}/merged_fastq", mode: params.publish_mode

    // conda option
    conda 'conda-forge::pigz=2.8'

    // input: grouped FASTQ files by sample ID
    input:
    tuple val(sample_id), path(r1_files), path(r2_files)

    // output: merged FASTQ files (as list to match FASTP input expectations)
    output:
    tuple val(sample_id), path("${sample_id}_merged_{R1,R2}.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_merge_log.txt"), emit: log

    script:
    """
    # verify we have matching numbers of R1 and R2 files
    R1_COUNT=\$(echo "${r1_files}" | wc -w)
    R2_COUNT=\$(echo "${r2_files}" | wc -w)

    if [ "\$R1_COUNT" != "\$R2_COUNT" ]; then
        echo "ERROR: Mismatched R1 (\$R1_COUNT) and R2 (\$R2_COUNT) file counts for sample ${sample_id}" > ${sample_id}_merge_log.txt
        exit 1
    fi

    echo "Merging \$R1_COUNT lane files for sample ${sample_id}" > ${sample_id}_merge_log.txt
    echo "R1 files: ${r1_files}" >> ${sample_id}_merge_log.txt
    echo "R2 files: ${r2_files}" >> ${sample_id}_merge_log.txt

    # sort files by name to ensure consistent lane order
    R1_SORTED=\$(echo "${r1_files}" | tr ' ' '\\n' | sort | tr '\\n' ' ')
    R2_SORTED=\$(echo "${r2_files}" | tr ' ' '\\n' | sort | tr '\\n' ' ')

    echo "Sorted R1 files: \$R1_SORTED" >> ${sample_id}_merge_log.txt
    echo "Sorted R2 files: \$R2_SORTED" >> ${sample_id}_merge_log.txt

    # concatenate R1 files
    echo "Concatenating R1 files..." >> ${sample_id}_merge_log.txt
    cat \$R1_SORTED > ${sample_id}_merged_R1.fastq.gz

    # concatenate R2 files
    echo "Concatenating R2 files..." >> ${sample_id}_merge_log.txt
    cat \$R2_SORTED > ${sample_id}_merged_R2.fastq.gz

    # verify output files exist and are non-empty
    if [ ! -s "${sample_id}_merged_R1.fastq.gz" ] || [ ! -s "${sample_id}_merged_R2.fastq.gz" ]; then
        echo "ERROR: Merged FASTQ files are empty or missing" >> ${sample_id}_merge_log.txt
        exit 1
    fi

    echo "Successfully merged files for sample ${sample_id}" >> ${sample_id}_merge_log.txt
    echo "Final R1 size: \$(du -h ${sample_id}_merged_R1.fastq.gz | cut -f1)" >> ${sample_id}_merge_log.txt
    echo "Final R2 size: \$(du -h ${sample_id}_merged_R2.fastq.gz | cut -f1)" >> ${sample_id}_merge_log.txt
    """
}