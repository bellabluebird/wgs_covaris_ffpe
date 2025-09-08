// modules/bwa_mem2_index.nf

process BWA_MEM2_INDEX {
    tag "$fasta.baseName"

    // publish alongside reference FASTA for future reuse
    publishDir "s3://bp-wgs-covaris-input-data/reference", mode: 'copy'

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
    def reference_name = fasta.getName()
    def index_exists = ['.amb', '.ann', '.bwt.2bit.64', '.pac', '.0123'].every { ext ->
        file("s3://bp-wgs-covaris-input-data/reference/${reference_name}${ext}").exists()
    }
    
    if (index_exists) {
        """
        echo "BWA-MEM2 index files already exist - using existing files"
        ln -s s3://bp-wgs-covaris-input-data/reference/${reference_name}.amb ${reference_name}.amb
        ln -s s3://bp-wgs-covaris-input-data/reference/${reference_name}.ann ${reference_name}.ann  
        ln -s s3://bp-wgs-covaris-input-data/reference/${reference_name}.bwt.2bit.64 ${reference_name}.bwt.2bit.64
        ln -s s3://bp-wgs-covaris-input-data/reference/${reference_name}.pac ${reference_name}.pac
        ln -s s3://bp-wgs-covaris-input-data/reference/${reference_name}.0123 ${reference_name}.0123
        """
    } else {
        """
        echo "BWA-MEM2 index files missing - creating new index"
        bwa-mem2 index ${fasta}
        """
    }
}
