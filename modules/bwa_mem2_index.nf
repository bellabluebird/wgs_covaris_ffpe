// modules/bwa_mem2_index.nf

process BWA_MEM2_INDEX {
    tag "$fasta.baseName"

    // publish alongside reference FASTA for future reuse
    publishDir "s3://bp-wgs-covaris-input-data/reference", mode: 'copy'

    // conda option (make sure AWS CLI is available, see note below)
    conda 'bioconda::bwa-mem2=2.2.1,bioconda::awscli'

    // input: reference fasta file
    input:
    path fasta

    // output: indexed reference files
    output:
    path "${fasta}", emit: fasta
    path "${fasta}.*", emit: index

    script:
    def reference_name = fasta.getName()
    def index_extensions = ['.amb', '.ann', '.bwt.2bit.64', '.pac', '.0123']
    
    """
    # Check if all index files exist locally
    if [ ! -f "${reference_name}.amb" ] || [ ! -f "${reference_name}.ann" ] || [ ! -f "${reference_name}.bwt.2bit.64" ] || [ ! -f "${reference_name}.pac" ] || [ ! -f "${reference_name}.0123" ]; then
        echo "BWA-MEM2 index files not found locally - checking S3..."
        
        # Download missing files from S3
        for ext in ${index_extensions.join(' ')}; do
            if [ ! -f "${reference_name}${ext}" ]; then
                echo "Downloading ${reference_name}${ext} from S3..."
                aws s3 cp s3://bp-wgs-covaris-input-data/reference/${reference_name}${ext} ./
            fi
        done

        # If still missing, create the index
        missing_files=false
        for ext in ${index_extensions.join(' ')}; do
            if [ ! -f "${reference_name}${ext}" ]; then
                missing_files=true
            fi
        done

        if [ "$missing_files" = true ]; then
            echo "Some index files missing after S3 copy - creating new index..."
            bwa-mem2 index ${fasta}
        else
            echo "All index files successfully copied from S3."
        fi
    else
        echo "All BWA-MEM2 index files exist locally - skipping download."
    fi
    """
}
