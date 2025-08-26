# AWS Setup Guide for Nextflow WGS Pipeline

This guide contains the complete, tested instructions for setting up your Nextflow WGS pipeline to run on AWS using GitHub Actions.

## Overview

This setup creates:
- **AWS Batch compute environment** that auto-scales EC2 instances
- **GitHub Actions workflow** that runs your pipeline from GitHub
- **OIDC authentication** for secure, keyless AWS access
- **S3 buckets** for input data, work files, and results

## Prerequisites

1. **AWS Account** with billing set up
2. **GitHub repository** with this pipeline code
3. **AWS CLI installed** and configured on your computer
4. **Basic familiarity** with AWS Console

## Part 1: Deploy AWS Infrastructure

### Step 1: Choose Your AWS Region

Pick a region close to you. This guide uses `us-east-2` (Ohio) but you can use any region:
- `us-east-1` (Virginia) - usually cheapest
- `us-east-2` (Ohio) - good alternative
- `us-west-2` (Oregon) - west coast
- `eu-west-1` (Ireland) - Europe

**Important**: Use the same region for all steps!

### Step 2: Deploy CloudFormation Stack

**Option A: AWS Console (Recommended for beginners)**

1. **Log into AWS Console**: https://aws.amazon.com/console/
2. **Select your region** in the top-right corner (e.g., us-east-2)
3. **Go to CloudFormation**: Search "CloudFormation" in the services search
4. **Create Stack**:
   - Click **Create stack** → **With new resources (standard)**
   - Choose **Upload a template file**
   - Select your `aws-infrastructure.yaml` file
   - Click **Next**
5. **Stack Details**:
   - **Stack name**: `bp-nextflow-infrastructure`
   - **MaxvCpus**: `8` (or leave default)
   - Click **Next**
6. **Configure Options**: 
   - Leave defaults, click **Next**
7. **Review**:
   - Check **"I acknowledge that AWS CloudFormation might create IAM resources"**
   - Click **Submit**
8. **Wait for completion** (~10-15 minutes)
   - Watch the **Events** tab for progress
   - Final status should be **CREATE_COMPLETE**

**Option B: AWS CLI**
```bash
aws cloudformation create-stack \
  --stack-name bp-nextflow-infrastructure \
  --template-body file://aws-infrastructure.yaml \
  --capabilities CAPABILITY_IAM \
  --region us-east-2
```

### Step 3: Create S3 Buckets

Create three S3 buckets in the **same region** as your CloudFormation stack:

```bash
# Replace us-east-2 with your chosen region
aws s3 mb s3://bp-wgs-covaris-nextflow-workdir --region us-east-2
aws s3 mb s3://bp-wgs-covaris-nextflow-results --region us-east-2  
aws s3 mb s3://bp-wgs-covaris-input-data --region us-east-2
```

**Verify buckets were created:**
```bash
aws s3 ls | grep bp-wgs-covaris
```

### Step 4: Verify Infrastructure

**Check that everything was created successfully:**

```bash
# Check CloudFormation stack
aws cloudformation describe-stacks --stack-name bp-nextflow-infrastructure --region us-east-2

# Check Batch compute environment
aws batch describe-compute-environments --region us-east-2

# Check job queue  
aws batch describe-job-queues --region us-east-2
```

**Expected outputs:**
- Stack status: `CREATE_COMPLETE`
- Compute environment state: `ENABLED` and `VALID`
- Job queue state: `ENABLED` and `VALID`

## Part 2: Set Up GitHub Authentication

### Step 5: Create OIDC Provider

**In AWS Console:**
1. **Go to IAM service**
2. **Click Identity providers** (left sidebar)
3. **Click Add provider**
4. **Configure provider**:
   - **Provider type**: `OpenID Connect`
   - **Provider URL**: `https://token.actions.githubusercontent.com`
   - **Audience**: `sts.amazonaws.com`
5. **Click Get thumbprint** (auto-fills)
6. **Click Add provider**

### Step 6: Create IAM Role for GitHub Actions

**Create the role:**
1. **In IAM**, click **Roles** → **Create role**
2. **Trusted entity**: Select **Web identity**
3. **Identity provider**: Choose `token.actions.githubusercontent.com`
4. **Audience**: Select `sts.amazonaws.com`
5. **GitHub organization**: Enter your GitHub username (e.g., `bpfeiffer`)
6. **GitHub repository**: Enter `wgs_covaris_ashg`
7. **Click Next**

**Skip policy selection:**
8. **Don't attach any policies yet**
9. **Click Next**
10. **Role name**: `GitHubActions-NextflowRole`
11. **Description**: `Role for running Nextflow pipelines from GitHub Actions`
12. **Click Create role**

### Step 7: Add Custom Policy to Role

1. **Find your new role** in IAM → Roles → `GitHubActions-NextflowRole`
2. **Go to Permissions tab** → **Add permissions** → **Create inline policy**
3. **Click JSON tab** and paste this policy:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "batch:DescribeComputeEnvironments",
        "batch:DescribeJobQueues",
        "batch:DescribeJobs",
        "batch:ListJobs",
        "batch:SubmitJob",
        "batch:TerminateJob",
        "batch:CancelJob"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:DeleteObject",
        "s3:ListBucket",
        "s3:GetBucketLocation"
      ],
      "Resource": [
        "arn:aws:s3:::bp-wgs-covaris-*",
        "arn:aws:s3:::bp-wgs-covaris-*/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "logs:CreateLogGroup",
        "logs:CreateLogStream",
        "logs:PutLogEvents",
        "logs:DescribeLogStreams"
      ],
      "Resource": "arn:aws:logs:*:*:log-group:/aws/batch/*"
    }
  ]
}
```

4. **Click Next**
5. **Policy name**: `NextflowPipelinePolicy`
6. **Click Create policy**

### Step 8: Get Role ARN and Add to GitHub

**Get the role ARN:**
1. **In your role summary**, copy the **ARN** (looks like: `arn:aws:iam::123456789012:role/GitHubActions-NextflowRole`)

**Add to GitHub secrets:**
1. **Go to your GitHub repository**
2. **Click Settings tab** → **Secrets and variables** → **Actions**
3. **Click New repository secret**
4. **Name**: `AWS_ROLE_ARN`
5. **Secret**: Paste the ARN you copied
6. **Click Add secret**

## Part 3: Test the Setup

### Step 9: Test GitHub Actions Workflow

**Run a test workflow:**
1. **Go to your GitHub repository**
2. **Click Actions tab**
3. **Click "Run Nextflow Pipeline on AWS"** (left sidebar)
4. **Click "Run workflow"** (green button)
5. **Fill in the form**:
   - **Input S3 path**: `s3://bp-wgs-covaris-input-data/samples`
   - **Output S3 path**: `s3://bp-wgs-covaris-nextflow-results/results`
   - **AWS region**: `us-east-2` (or your chosen region)
   - **Nextflow version**: `23.10.1`
6. **Click "Run workflow"**

**Expected results:**
- ✅ Checkout repository
- ✅ Configure AWS credentials  
- ✅ Setup Java
- ✅ Install Nextflow
- ✅ Validate AWS setup
- ❌ Run Nextflow pipeline (fails because no input data - this is expected!)
- ⚠️ Upload reports (shows "No files found" - this is normal)

**If the "Validate AWS setup" step passes, your entire setup is working!**

## Part 4: Running Your Pipeline

### Step 10: Upload Data and Run Pipeline

**Upload your FASTQ files:**
```bash
# Example: upload paired-end sequencing files
aws s3 cp sample1_R1.fastq.gz s3://bp-wgs-covaris-input-data/samples/
aws s3 cp sample1_R2.fastq.gz s3://bp-wgs-covaris-input-data/samples/
aws s3 cp sample2_R1.fastq.gz s3://bp-wgs-covaris-input-data/samples/
aws s3 cp sample2_R2.fastq.gz s3://bp-wgs-covaris-input-data/samples/
```

**Required file naming:**
- Files must be paired-end and compressed (.gz)
- Supported patterns: `*_R1.fastq.gz` & `*_R2.fastq.gz` OR `*_1.fastq.gz` & `*_2.fastq.gz`

**Run the pipeline:**
1. **Go to GitHub Actions** again
2. **Run the workflow** with your data
3. **Monitor progress** in the Actions tab
4. **Download reports** from the artifacts section when complete
5. **Check results** in your S3 results bucket

## Cost Management

### Expected Costs

**Typical costs for small datasets (2-4 samples):**
- **Compute**: $1-5 per run (automatically scales down when idle)
- **Storage**: $0.10-1 per month (S3 storage for data and results)
- **Infrastructure**: $0 (no charges when not running)

### Cost Controls Built In

1. **MaxvCpus**: Limited to 8 CPU cores maximum
2. **Auto-scaling**: Instances automatically terminate when jobs finish
3. **Optimal instances**: AWS chooses cheapest available instance types
4. **No persistent resources**: Only pay for what you use

### Monitor Costs

1. **AWS Billing Console**: Check costs daily for the first week
2. **Set up billing alerts**: Get notified if costs exceed $10/month
3. **Tag resources**: All resources are tagged for easy cost tracking

## Troubleshooting

### Common Issues

**1. CloudFormation stack fails:**
- Check the Events tab for specific error messages
- Try a different AWS region
- Verify your account has proper permissions

**2. "ComputeEnvironment did not stabilize":**
- Delete the failed stack completely
- Wait 5 minutes, then redeploy
- Try using a different region

**3. GitHub Actions authentication fails:**
- Verify the OIDC provider URL is exactly: `https://token.actions.githubusercontent.com`
- Check that AWS_ROLE_ARN secret matches your role ARN exactly
- Ensure the role trust policy includes your exact GitHub username/repo

**4. "Cannot find any paired FASTQ files":**
- Check file naming conventions (`*_R1.fastq.gz`, `*_R2.fastq.gz`)
- Verify files are uploaded to the correct S3 path
- Ensure files are gzip compressed

### Getting Help

**Check logs in this order:**
1. **GitHub Actions logs** - for workflow and authentication issues
2. **AWS Batch console** - for job execution issues  
3. **CloudWatch logs** - for detailed container logs
4. **S3 bucket contents** - to verify file uploads/downloads

**Useful AWS CLI commands:**
```bash
# Check batch job status
aws batch list-jobs --job-queue nextflow-compute-queue --region us-east-2

# Monitor running jobs
aws batch describe-jobs --jobs JOB-ID --region us-east-2

# Check S3 bucket contents  
aws s3 ls s3://bp-wgs-covaris-input-data/samples/ --recursive
```

## Cleanup

### When You're Done

**To avoid ongoing costs, clean up resources:**

```bash
# Delete S3 bucket contents (can't delete non-empty buckets)
aws s3 rm s3://bp-wgs-covaris-nextflow-workdir --recursive
aws s3 rm s3://bp-wgs-covaris-nextflow-results --recursive  
aws s3 rm s3://bp-wgs-covaris-input-data --recursive

# Delete S3 buckets
aws s3 rb s3://bp-wgs-covaris-nextflow-workdir
aws s3 rb s3://bp-wgs-covaris-nextflow-results
aws s3 rb s3://bp-wgs-covaris-input-data

# Delete CloudFormation stack (removes all AWS Batch, VPC, etc.)
aws cloudformation delete-stack --stack-name bp-nextflow-infrastructure --region us-east-2
```

**Keep these for future use:**
- OIDC provider (reusable for other projects)
- IAM role (reusable, no ongoing costs)
- GitHub repository and workflow

## Success Checklist

By the end of this setup, you should have:

- ✅ CloudFormation stack with status `CREATE_COMPLETE`
- ✅ 3 S3 buckets created
- ✅ OIDC provider configured
- ✅ IAM role with proper permissions
- ✅ GitHub secret `AWS_ROLE_ARN` configured
- ✅ Successful test run of GitHub Actions workflow
- ✅ "Validate AWS setup" step passing

**You now have a production-ready, scalable genomics pipeline running on AWS!** 🎉
