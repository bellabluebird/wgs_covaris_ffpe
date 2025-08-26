# AWS Setup Guide for Nextflow WGS Pipeline

This guide will help you set up and run the Nextflow WGS pipeline on AWS using GitHub Actions.

## Prerequisites

1. AWS Account with appropriate permissions
2. GitHub repository with this pipeline code
3. AWS CLI configured locally (for initial setup)

## Step 1: Create S3 Buckets

Create three S3 buckets for your pipeline:

```bash
# Replace 'your-unique-prefix' with your own prefix
aws s3 mb s3://your-unique-prefix-nextflow-workdir
aws s3 mb s3://your-unique-prefix-nextflow-results  
aws s3 mb s3://your-unique-prefix-input-data
```

## Step 2: Deploy AWS Infrastructure

Deploy the CloudFormation stack to set up AWS Batch:

```bash
aws cloudformation create-stack \
  --stack-name nextflow-infrastructure \
  --template-body file://aws-infrastructure.yaml \
  --capabilities CAPABILITY_IAM \
  --parameters \
    ParameterKey=InstanceTypes,ParameterValue="c5.large,c5.xlarge,c5.2xlarge" \
    ParameterKey=MaxvCpus,ParameterValue=1000
```

Monitor the stack creation:
```bash
aws cloudformation wait stack-create-complete --stack-name nextflow-infrastructure
```

## Step 3: Update Configuration

Edit `nextflow.config` and replace the placeholder S3 paths:

```groovy
// Replace these with your actual bucket names
workDir = 's3://your-unique-prefix-nextflow-workdir/work'
params {
    outdir = 's3://your-unique-prefix-nextflow-results/results'
    input_dir = 's3://your-unique-prefix-input-data/samples'
}
```

## Step 4: Set up GitHub Secrets

### Option A: Using IAM User (Simpler but less secure)

1. Create an IAM user with the following policies:
   - `AWSBatchFullAccess`
   - `AmazonS3FullAccess`
   - `CloudWatchLogsFullAccess`

2. Add these secrets to your GitHub repository:
   ```
   AWS_ACCESS_KEY_ID: <your-access-key>
   AWS_SECRET_ACCESS_KEY: <your-secret-key>
   ```

### Option B: Using OIDC (Recommended for production)

1. Create an OIDC provider in AWS IAM
2. Create a role that can be assumed by GitHub Actions
3. Add this secret to your GitHub repository:
   ```
   AWS_ROLE_ARN: arn:aws:iam::ACCOUNT_ID:role/github-actions-role
   ```

## Step 5: Upload Input Data

Upload your FASTQ files to the input S3 bucket:

```bash
# Upload paired-end FASTQ files
aws s3 sync /path/to/your/fastq/files/ s3://your-unique-prefix-input-data/samples/
```

Ensure your FASTQ files follow the naming convention:
- `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`
- or `sample1_1.fastq.gz` and `sample1_2.fastq.gz`

## Step 6: Run the Pipeline

1. Go to your GitHub repository
2. Navigate to Actions → "Run Nextflow Pipeline on AWS"
3. Click "Run workflow"
4. Fill in the parameters:
   - **Input S3 Path**: `s3://your-unique-prefix-input-data/samples`
   - **Output S3 Path**: `s3://your-unique-prefix-nextflow-results/results`
   - **AWS Region**: `us-east-1` (or your preferred region)

## Step 7: Monitor and Results

- Monitor progress in GitHub Actions
- Check AWS Batch console for job status
- Results will be stored in your output S3 bucket
- Download the Nextflow reports from GitHub Actions artifacts

## Cost Optimization Tips

1. Use Spot instances by adding to `aws-infrastructure.yaml`:
   ```yaml
   BidPercentage: 50  # Use spot instances at 50% of on-demand price
   ```

2. Set appropriate instance limits:
   ```yaml
   MaxvCpus: 100  # Adjust based on your needs
   ```

3. Clean up resources when not in use:
   ```bash
   aws cloudformation delete-stack --stack-name nextflow-infrastructure
   ```

## Troubleshooting

### Common Issues:

1. **Permission Errors**: Check IAM roles and S3 bucket policies
2. **Queue Not Available**: Verify AWS Batch compute environment is ENABLED
3. **Container Pull Errors**: Check internet connectivity in your VPC
4. **S3 Access Denied**: Verify bucket names and IAM permissions

### Debug Commands:

```bash
# Check Batch compute environments
aws batch describe-compute-environments

# Check job queues
aws batch describe-job-queues

# Monitor running jobs
aws batch list-jobs --job-queue nextflow-compute-queue --job-status RUNNING
```

## Cleanup

To avoid ongoing costs, delete resources when finished:

```bash
# Delete CloudFormation stack (removes Batch, VPC, etc.)
aws cloudformation delete-stack --stack-name nextflow-infrastructure

# Empty and delete S3 buckets
aws s3 rm s3://your-unique-prefix-nextflow-workdir --recursive
aws s3 rb s3://your-unique-prefix-nextflow-workdir

aws s3 rm s3://your-unique-prefix-nextflow-results --recursive  
aws s3 rb s3://your-unique-prefix-nextflow-results

aws s3 rm s3://your-unique-prefix-input-data --recursive
aws s3 rb s3://your-unique-prefix-input-data
```