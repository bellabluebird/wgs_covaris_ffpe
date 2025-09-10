#!/usr/bin/env nextflow

// WGS Quality Control & Analysis Pipeline with BQSR
// Requires: --input_dir, --reference, --known_sites

// enable modular functions
nextflow.enable.dsl = 2

// process modules
include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'
include { FASTP } from './modules/fastp.nf'
//include { BWA_MEM2_INDEX } from './modules/bwa_mem2_index.nf'
include { BWA } from './modules/bwa.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKED } from './modules/samtools_index.nf'
include { QUALIMAP } from './modules/qualimap.nf'
include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates.nf'
include { PICARD_COLLECTINSERTSIZEMETRICS } from './modules/picard_insert_size.nf'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { GATK_BASERECALIBRATOR } from './modules/gatk_baserecalibrator.nf'
include { GATK_APPLYBQSR } from './modules/gatk_applybqsr.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BQSR } from './modules/samtools_index.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
include { GATK_GENOTYPEGVCFS } from './modules/gatk_genotypegvcfs.nf'
include { GATK_VARIANTFILTRATION } from './modules/gatk_variantfiltration.nf'
include { BCFTOOLS_STATS } from './modules/bcftools_stats.nf'
include { MULTIQC } from './modules/multiqc.nf'

// parameters - all data-specific parameters must be provided explicitly
params.outdir = "./results"
params.publish_mode = 'copy'

// required parameters (no defaults - must be specified via CLI or profile)
params.input_dir = null
params.reference = null
params.known_sites = null

// input validation
log.info "Input directory: ${params.input_dir}"
log.info "Reference genome: ${params.reference}"
log.info "Known sites: ${params.known_sites}"

// create input channel from FASTQ files in input directory
ch_input = Channel.fromFilePairs("${params.input_dir}/*_{R1,R2,1,2}.{fastq,fq}{,.gz}", checkIfExists: false)
    .ifEmpty { error "No paired FASTQ files found in ${params.input_dir}. Expected pattern: *_{R1,R2,1,2}.{fastq,fq}{,.gz}" }

// create reference channel from explicit parameter
ch_reference_fasta = Channel.fromPath(params.reference)
    .ifEmpty { error "Reference genome not found at ${params.reference}" }

// create known sites channel for BQSR (use compatible known sites from hg38)
ch_known_sites = (params.known_sites ?
    Channel.fromPath(params.known_sites.split(',').collect { it.trim() }) :
    Channel.fromPath("s3://bp-wgs-covaris-input-data/reference/known_sites_hg38/Homo_sapiens_assembly38.dbsnp138.vcf,s3://bp-wgs-covaris-input-data/reference/known_sites_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.idx,s3://bp-wgs-covaris-input-data/reference/known_sites_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,s3://bp-wgs-covaris-input-data/reference/known_sites_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi".split(',').collect { it.trim() }))
ch_known_sites = ch_known_sites
    .ifEmpty { error "No known sites files found at specified paths" }
    .collect() 

// main workflow
workflow {
    // raw fastqc
    fastqc_raw = FASTQC(ch_input)
    
    // preprocessing with fastp
    fastp_results = FASTP(ch_input)
    
    // post-trim fastqc on cleaned reads
    fastqc_trimmed = FASTQC_TRIMMED(fastp_results.reads)
    
    // use compatible hg38_broad references (chr1, chr2 format - matches VCF)
    ch_reference_indexed = Channel.fromPath("s3://bp-wgs-covaris-input-data/reference/hg38_broad/Homo_sapiens_assembly38.{fasta,fasta.64.*}")
        .collect()
    
    // create separate channel with GATK-required reference files (FASTA, .fai, .dict)
    ch_reference_fasta_gatk = Channel.fromPath("s3://bp-wgs-covaris-input-data/reference/hg38_broad/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}")
        .collect()
    
    // alignment to reference genome
    bwa_results = BWA(fastp_results.reads, ch_reference_indexed)
    
    // alignment statistics
    samtools_stats = SAMTOOLS_STATS(bwa_results.bam)
    
    // index aligned BAM files
    samtools_index = SAMTOOLS_INDEX(bwa_results.bam)
    
    // mark duplicates (using indexed BAM)
    markduplicates_results = PICARD_MARKDUPLICATES(samtools_index.bam_bai)
    
    // index the marked BAM files
    samtools_index_marked = SAMTOOLS_INDEX_MARKED(markduplicates_results.bam)
    
    // insert size metrics (on duplicate-marked BAM)
    insert_size_metrics = PICARD_COLLECTINSERTSIZEMETRICS(markduplicates_results.bam)
    
    // coverage analysis with mosdepth (on duplicate-marked indexed BAM)
    coverage_results = MOSDEPTH(samtools_index_marked.bam_bai)
    
    // base quality score recalibration with GATK
    bqsr_table = GATK_BASERECALIBRATOR(samtools_index_marked.bam_bai, ch_reference_fasta_gatk, ch_known_sites)
    
    // apply BQSR to get recalibrated BAM (join BAM and BQSR table by sample_id)
    bam_bqsr_joined = samtools_index_marked.bam_bai.join(bqsr_table.recal_table)
    bqsr_results = GATK_APPLYBQSR(bam_bqsr_joined, ch_reference_fasta_gatk)
    
    // index the recalibrated BAM files
    samtools_index_bqsr = SAMTOOLS_INDEX_BQSR(bqsr_results.bam)
    
    // quality metrics with Qualimap (on recalibrated indexed BAM)
    qualimap_results = QUALIMAP(samtools_index_bqsr.bam_bai)
    
    // variant calling with GATK HaplotypeCaller (generate GVCFs for joint genotyping)
    haplotypecaller_results = GATK_HAPLOTYPECALLER(samtools_index_bqsr.bam_bai, ch_reference_fasta_gatk)
    
    // debug to see what's being passed to the new channel
    haplotypecaller_results.gvcf.view { row -> "HaplotypeCaller Output: $row" }

    // collect all GVCFs for joint genotyping
    all_gvcfs = haplotypecaller_results.gvcf.map { id, gvcf -> gvcf }.unique().collect()
    
    // debug to see the collected GVCFs 
    all_gvcfs.view { list -> "All GVCFs for joint genotyping → ${list}" }

    // joint genotyping across all samples
    joint_vcf = GATK_GENOTYPEGVCFS(all_gvcfs, ch_reference_fasta_gatk)
    
    // apply hard filtering to variants
    filtered_vcf = GATK_VARIANTFILTRATION(joint_vcf.vcf, joint_vcf.vcf_index, ch_reference_fasta_gatk)
    
    // generate variant statistics on filtered VCF
    variant_stats = BCFTOOLS_STATS(filtered_vcf.vcf.map { vcf -> ["joint_variants", vcf] })

    // collect all reports for multiqc
    multiqc_input = fastqc_raw.zip.map { id, file -> file }
    .mix(fastqc_trimmed.zip.map { id, file -> file })
    .mix(fastp_results.json.map { id, file -> file })
    .mix(samtools_stats.stats.map { id, file -> file })
    .mix(markduplicates_results.metrics.map { id, file -> file })
    .mix(insert_size_metrics.metrics.map { id, file -> file })
    .mix(coverage_results.summary.map { id, file -> file })
    .mix(qualimap_results.results.map { id, dir -> dir })
    .mix(variant_stats.stats.map { id, file -> file })
    .collect()

    // multiqc report
    MULTIQC(multiqc_input)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Results directory: ${params.outdir}"
}