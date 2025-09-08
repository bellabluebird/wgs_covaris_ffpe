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
        echo "BWA-MEM2 index files already exist - copying from S3"
        aws s3 cp s3://bp-wgs-covaris-input-data/reference/${reference_name}.amb .
        aws s3 cp s3://bp-wgs-covaris-input-data/reference/${reference_name}.ann .
        aws s3 cp s3://bp-wgs-covaris-input-data/reference/${reference_name}.bwt.2bit.64 .
        aws s3 cp s3://bp-wgs-covaris-input-data/reference/${reference_name}.pac .
        aws s3 cp s3://bp-wgs-covaris-input-data/reference/${reference_name}.0123 .
        """
    } else {
        """
        echo "BWA-MEM2 index files missing - creating new index"
        bwa-mem2 index ${fasta}
        """
    }
}
