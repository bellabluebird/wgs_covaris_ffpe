# nextflow wgs pipeline

a simple nextflow pipeline for quality control of whole genome sequencing data

## what it does

this pipeline takes your raw dna sequencing files (fastq) and:
1. checks quality with fastqc
2. trims low quality bases with fastp  
3. checks quality again after trimming
4. creates a summary report with multiqc

## what you need

- paired-end fastq files (like `sample_R1.fastq.gz` and `sample_R2.fastq.gz`)
- nextflow installed
- either conda, docker, or aws batch set up

## how to run locally

```bash
# with conda
nextflow run main.nf -profile conda --input_dir /path/to/fastq/files

# with docker  
nextflow run main.nf -profile docker --input_dir /path/to/fastq/files
```

## how to run on aws

1. follow the setup guide in `AWS_SETUP.md`
2. push your code to github
3. go to github actions and click "run workflow"
4. fill in your s3 bucket paths
5. wait for it to finish and download the reports

## files in this folder

- `main.nf` - the main pipeline script
- `nextflow.config` - configuration for different environments
- `modules/` - individual process definitions
  - `fastqc.nf` - quality control checks
  - `fastp.nf` - adapter trimming and filtering
  - `multiqc.nf` - summary report generation
- `.github/workflows/aws-nextflow.yml` - github actions workflow for aws
- `aws-infrastructure.yaml` - cloudformation template for aws setup
- `AWS_SETUP.md` - detailed aws setup instructions

## results

your results will be in the `results/` folder (or s3 bucket) with:
- `fastqc/` - quality reports before trimming
- `fastp/` - trimmed files and trimming reports
- `multiqc/` - combined summary report (start here!)

## troubleshooting

- make sure your fastq files are paired and follow naming conventions
- check that you have enough memory/cpu allocated
- look at the nextflow log files if something fails
- for aws issues, check the setup guide and aws console

## need help?

- check the nextflow documentation: https://nextflow.io/docs/latest/
- for aws batch help: see `AWS_SETUP.md`
- for pipeline issues: check the `.nextflow.log` file