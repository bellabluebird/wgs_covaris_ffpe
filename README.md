# nextflow wgs pipeline

a nextflow pipeline for whole genome sequencing analysis on aws batch

## what it does

this pipeline performs comprehensive wgs analysis through these modules:

### quality control & preprocessing
1. **fastqc** - analyzes raw fastq files → generates quality reports (html/zip)
2. **fastp** - trims adapters and low-quality bases from fastq files → produces cleaned fastq files + json reports
3. **fastqc** - re-analyzes cleaned fastq files → generates post-trim quality reports

### alignment & processing  
4. **bwa_mem2_index** - indexes reference fasta file → creates bwa-mem2 index files
5. **bwa_mem2** - aligns cleaned fastq pairs to reference → produces sorted bam files
6. **samtools_index** - indexes aligned bam files → creates bai index files
7. **samtools_stats** - analyzes bam alignment statistics → generates stats files

### duplicate marking & metrics
8. **picard_markduplicates** - marks pcr duplicates in bam files → produces duplicate-marked bam + metrics
9. **samtools_index** - indexes duplicate-marked bam → creates new bai files
10. **picard_insert_size** - calculates insert size metrics from bam → generates insert size reports

### base quality recalibration
11. **gatk_baserecalibrator** - analyzes base quality patterns using known sites vcf → creates recalibration table
12. **gatk_applybqsr** - applies recalibration to bam file → produces bqsr-corrected bam + index

### final quality assessment
13. **mosdepth** - calculates coverage depth from final bam → generates coverage summaries
14. **qualimap** - performs comprehensive bam quality analysis → creates alignment quality reports
15. **multiqc** - aggregates all reports from steps 1-14 → produces final combined html report

## what you need

- paired-end fastq files (like `sample_R1.fastq.gz` and `sample_R2.fastq.gz`)
- reference genome fasta file
- known sites vcf file for bqsr (like dbsnp)
- aws account with batch setup (see `AWS_SETUP.md`)

## how to run

1. follow the setup guide in `AWS_SETUP.md`
2. upload your files to s3:
   - fastq files → `s3://your-bucket/samples/`
   - reference genome → `s3://your-bucket/reference/`
   - known sites vcf → `s3://your-bucket/known_sites/`
3. go to github actions and click "run workflow"
4. fill in your s3 paths:
   - input s3 path: `s3://your-bucket/samples`
   - reference s3 path: `s3://your-bucket/reference/genome.fasta`
   - known sites: `s3://your-bucket/known_sites/dbsnp.vcf.gz`
   - output s3 path: `s3://your-bucket/results`
5. wait for completion and check results in s3

## files in this folder

- `main.nf` - main pipeline script
- `nextflow.config` - aws batch configuration
- `modules/` - individual process definitions
- `.github/workflows/aws-nextflow.yml` - github actions workflow
- `aws/aws-infrastructure.yaml` - cloudformation template
- `aws/AWS_SETUP.md` - detailed setup instructions

## results

your results will be in your s3 output bucket with:
- `fastqc/` - quality reports
- `fastp/` - trimmed files and reports
- `alignment/` - bam files
- `qualimap/` - alignment quality metrics
- `mosdepth/` - coverage analysis
- `multiqc/` - combined summary report (start here!)

## troubleshooting

- make sure fastq files are paired: `*_R1.fastq.gz` and `*_R2.fastq.gz`
- verify all s3 paths are correct in github actions
- check aws batch job status in aws console
- look at github actions logs for errors
- see `AWS_SETUP.md` for detailed troubleshooting

## need help?

- aws setup issues: check `aws/AWS_SETUP.md`
- pipeline issues: check github actions logs
- nextflow docs: https://nextflow.io/docs/latest/


9/8/25: this version has successful indexing and saving to references! from here, we will switch to using the index that i've stored for grch38