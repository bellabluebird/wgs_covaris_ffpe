// modules/bwa_mem2_index.nf

process BWA_MEM2_INDEX {
    tag "$fasta.baseName"

    // publish alongside the fasta file for referencing in future runs
    publishDir { fasta.parent }, mode: 'copy'

    // conda option
    conda 'bioconda::bwa-mem2=2.2.1'

    // input: reference fasta file
    input:
    path fasta

    // output: indexed reference files
    output:
    path "${fasta}", emit: fasta
    path "${fasta}.*", emit: index

    script:
    """
    # create BWA-MEM2 index
    bwa-mem2 index ${fasta}
    """
}
