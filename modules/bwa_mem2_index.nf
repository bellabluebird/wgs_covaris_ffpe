process BWA_MEM2_INDEX {
    tag "$fasta.baseName"

    // publish alongside reference FASTA for future reuse
    publishDir "s3://bp-wgs-covaris-input-data/reference", mode: 'copy'

    conda 'bioconda::bwa-mem2=2.2.1'

    input:
    path fasta

    output:
    path "${fasta}", emit: fasta
    path "${fasta}.*", emit: index

    script:
    def reference_name = fasta.getName()
    
    """
    # Check if all index files exist locally
    if [ -f "${reference_name}.amb" ] && \
       [ -f "${reference_name}.ann" ] && \
       [ -f "${reference_name}.bwt.2bit.64" ] && \
       [ -f "${reference_name}.pac" ] && \
       [ -f "${reference_name}.0123" ]; then
        echo "All BWA-MEM2 index files exist locally - skipping indexing."
    else
        echo "Index files missing - creating BWA-MEM2 index."
        bwa-mem2 index ${fasta}
    fi
    """
}
