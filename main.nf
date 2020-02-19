#!/usr/bin/env nextflow

/*
Developed by the Genome Institute of Singapore for
the National Precision Medicine Programme

Copyright: 2020 Genome Institute of Singapore
License: The MIT License (MIT)

See LICENSE for more copyright information
*/

/*
----------------------------------------------------------------------
FUNCTIONS
----------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    Usage: nextflow run main.nf -config nextflow.config -params-file sample_params.yaml 
                                -profile nscc -work-dir ./ --outdir ./
                                [-resume] [--keep_workdir] [--help]
    
    Options:
    -config           Generic workflow settings
    -params-file      Sample-specific settings
    -profile          Job launch settings
    -resume           Re-use existing results (optional, omit to re-run from scratch)
    --keep_workdir    Keep work directory (optional, omit for auto-deletion)
    --help            Print this help message
    """.stripIndent()
}

def nextflowMessage() {
    log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def minimalInformationMessage() {
    log.info "User name    : " + workflow.userName
    log.info "Command Line : " + workflow.commandLine
    log.info "Project Dir  : " + workflow.projectDir
    log.info "Launch Dir   : " + workflow.launchDir
    log.info "Work Dir     : " + workflow.workDir
    log.info "Results Dir  : " + params.publishdir
    log.info "Info Dir     : " + params.infodir
    log.info "Profile      : " + workflow.profile
}

def startMessage() {
    this.nextflowMessage()
    this.minimalInformationMessage()
}

/*
----------------------------------------------------------------------
USAGE
----------------------------------------------------------------------
*/

if (params.help) exit 0, helpMessage()

/*
----------------------------------------------------------------------
LAUNCH INFO
----------------------------------------------------------------------
*/

startMessage()

/*
----------------------------------------------------------------------
INPUT CHANNELS
----------------------------------------------------------------------
*/

Channel
    .fromPath(params.cram)
    .into { cram_ch_cram_to_bam
          ; cram_ch_samtools_flagstat }

Channel
    .fromPath(params.vcf)
    .into { vcf_ch_bcftools_stats
          ; vcf_ch_bcftools_gtcheck
          ; vcf_ch_picard_collect_variant_calling_metrics_vcf
          ; vcf_ch_picard_collect_variant_calling_metrics_gvcf }

if (params.pst_vcf) {
    Channel
        .fromPath(params.pst_vcf)
        .set { pst_vcf_ch_bcftools_gtcheck }
} else {
    Channel
        .empty()
        .set { pst_vcf_ch_bcftools_gtcheck }
}

Channel
    .fromPath(params.ref_fa)
    .into { ref_fa_ch_cram_to_bam
          ; ref_fa_ch_picard_collect_alignment_summary_metrics
          ; ref_fa_ch_picard_collect_wgs_metrics
          ; ref_fa_ch_picard_collect_gc_bias_metrics
          ; ref_fa_ch_verifybamid2 }

Channel
    .fromPath(params.dbsnp_vcf)
    .into { dbsnp_vcf_ch_picard_collect_variant_calling_metrics_vcf
          ; dbsnp_vcf_ch_picard_collect_variant_calling_metrics_gvcf }
/*
----------------------------------------------------------------------
PROCESSES
---------------------------------------------------------------------
*/

process cram_to_bam {

    input:
    file ref_fa from ref_fa_ch_cram_to_bam
    file cram from cram_ch_cram_to_bam

    output:
    file "*" into cram_to_bam_ch_samtools_stats \
                , cram_to_bam_ch_picard_collect_quality_yield_metrics \
                , cram_to_bam_ch_picard_collect_alignment_summary_metrics \
                , cram_to_bam_ch_picard_collect_wgs_metrics \
                , cram_to_bam_ch_sg10k_cov_062017 \
                , cram_to_bam_ch_picard_collect_insert_size_metrics \
                , cram_to_bam_ch_picard_collect_gc_bias_metrics \
                , cram_to_bam_ch_verifybamid2

    script:
    """
    samtools view -h -T ${ref_fa} -b ${cram} -o ${params.sample_id}.bam
    samtools index ${params.sample_id}.bam ${params.sample_id}.bam.bai
    """

}

process samtools_stats {

    publishDir "${params.publishdir}/samtools", mode: "copy"

    input:
    file "*" from cram_to_bam_ch_samtools_stats

    output:
    file "*" into samtools_stats_ch \
                , samtools_stats_ch_plotbamstats

    script:
    """
    samtools stats ${params.sample_id}.bam > ${params.sample_id}.stats
    """

}

process samtools_flagstat {

    publishDir "${params.publishdir}/samtools", mode: "copy"

    input:
    file cram from cram_ch_samtools_flagstat

    output:
    file "*" into samtools_flagstat_ch

    script:
    """
    samtools flagstat ${cram} > ${params.sample_id}.flagstat
    """

}

process bcftools_stats {

    publishDir "${params.publishdir}/bcftools", mode: "copy"

    input:
    file vcf from vcf_ch_bcftools_stats

    output:
    file "*.bcftools_stats.txt" into bcftools_stats_ch

    script:
    """
    bcftools index --tbi ${vcf}
    bcftools stats ${vcf} > ${params.sample_id}.bcftools_stats.txt
    """

}

process bcftools_gtcheck {

    publishDir "${params.publishdir}/bcftools", mode: "copy"

    input:
    file vcf from vcf_ch_bcftools_gtcheck
    file pst_vcf from pst_vcf_ch_bcftools_gtcheck

    output:
    file "*.bcftools_gtcheck.txt" into bcftools_gtcheck_ch

    when:
    params.pst_vcf != null

    script:
    """
    bcftools index --tbi ${vcf}
    bcftools index --tbi ${pst_vcf}
    bcftools gtcheck --GTs-only 1 --genotypes ${pst_vcf} ${vcf} > ${params.sample_id}.bcftools_gtcheck.txt
    """
}

process picard_collect_quality_yield_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file "*" from cram_to_bam_ch_picard_collect_quality_yield_metrics

    output:
    file "*" into picard_collect_quality_yield_metrics_ch

    script:
    """
    picard CollectQualityYieldMetrics \
        I=${params.sample_id}.bam \
        O=${params.sample_id}.quality_yield_metrics.txt
    """

}

process picard_collect_alignment_summary_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_picard_collect_alignment_summary_metrics
    file "*" from cram_to_bam_ch_picard_collect_alignment_summary_metrics

    output:
    file "*" into picard_collect_alignment_summary_metrics_ch

    script:
    """
    picard CollectAlignmentSummaryMetrics \
        R=${ref_fa} \
        I=${params.sample_id}.bam \
        O=${params.sample_id}.alignment_summary_metrics.txt
    """

}

process picard_collect_wgs_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_picard_collect_wgs_metrics
    file "*" from cram_to_bam_ch_picard_collect_wgs_metrics

    output:
    file "*" into picard_collect_wgs_metrics_ch

    script:
    """
    picard CollectWgsMetrics \
        R=${ref_fa} \
        I=${params.sample_id}.bam \
        O=${params.sample_id}.wgs_metrics.txt
    """

}

process sg10k_cov_062017 {

    publishDir "${params.publishdir}/sg10k_cov_062017", mode: "copy"

    input:
    file "*" from cram_to_bam_ch_sg10k_cov_062017

    output:
    file "*" into sg10k_cov_062017_ch

    script:
    """
    sg10k-cov-062017.sh ${params.sample_id}.bam > ${params.sample_id}.sg10k_cov_062017.txt
    """

}

process picard_collect_insert_size_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file "*" from cram_to_bam_ch_picard_collect_insert_size_metrics

    output:
    file "*" into picard_collect_insert_size_metrics_ch

    script:
    """
    picard CollectInsertSizeMetrics \
        I=${params.sample_id}.bam \
        O=${params.sample_id}.insert_size_metrics.txt \
        H=${params.sample_id}.insert_size_histogram.pdf \
        M=0.5
    """

}

process picard_collect_gc_bias_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_picard_collect_gc_bias_metrics
    file "*" from cram_to_bam_ch_picard_collect_gc_bias_metrics

    output:
    file "*" into picard_collect_gc_bias_metrics_ch

    script:
    """
    picard CollectGcBiasMetrics \
        I=${params.sample_id}.bam \
        O=${params.sample_id}.gc_bias_metrics.txt \
        CHART=${params.sample_id}.gc_bias_metrics.pdf \
        S=${params.sample_id}.gc_bias_summary_metrics.txt \
        R=${ref_fa}
    """

}

process picard_collect_variant_calling_metrics_vcf {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file vcf from vcf_ch_picard_collect_variant_calling_metrics_vcf
    file dbsnp_vcf from dbsnp_vcf_ch_picard_collect_variant_calling_metrics_vcf

    output:
    file "*metrics" into picard_collect_variant_calling_metrics_vcf_ch

    when:
    vcf.name =~ /.*\.vcf.gz$/

    script:
    """
    bcftools index --tbi ${vcf}
    picard CollectVariantCallingMetrics \
        I=${vcf} \
        O=${params.sample_id} \
        DBSNP=${dbsnp_vcf}
    """

}

process picard_collect_variant_calling_metrics_gvcf {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file vcf from vcf_ch_picard_collect_variant_calling_metrics_gvcf
    file dbsnp_vcf from dbsnp_vcf_ch_picard_collect_variant_calling_metrics_gvcf

    output:
    file "*metrics" into picard_collect_variant_calling_metrics_gvcf_ch

    when:
    vcf.name =~ /.*\.gvcf.gz$/

    script:
    """
    bcftools view -v snps,indels -O z -o ${params.sample_id}.vcf.gz ${vcf}
    bcftools index --tbi ${params.sample_id}.vcf.gz
    picard CollectVariantCallingMetrics \
        I=${params.sample_id}.vcf.gz \
        O=${params.sample_id} \
        DBSNP=${dbsnp_vcf}
    """
}

process verifybamid2 {

    publishDir "${params.publishdir}/verifybamid2", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_verifybamid2
    file "*" from cram_to_bam_ch_verifybamid2

    output:
    file "*" into verifybamid2_ch

    script:
    """
    VerifyBamID --SVDPrefix ${params.vbi2_svdprefix} --Reference ${params.ref_fa} --BamFile ${params.sample_id}.bam
    """

}

process plot_bamstats {

    publishDir "${params.publishdir}/plot_bamstats", mode: "copy"

    input:
    file "samtools/*" from samtools_stats_ch_plotbamstats

    output:
    file "${params.sample_id}*" into plot_bamstats_ch

    script:
    """
    plot-bamstats -p ${params.sample_id} samtools/${params.sample_id}.stats
    """

}

process multiqc {

    publishDir "${params.publishdir}/multiqc", mode: "copy"

    input:
    file "samtools/*" from samtools_stats_ch
    file "samtools/*" from samtools_flagstat_ch
    file "bcftools/*" from bcftools_stats_ch
    file "bcftools/*" from bcftools_gtcheck_ch.collect().ifEmpty([])
    file "picard/*" from picard_collect_quality_yield_metrics_ch
    file "picard/*" from picard_collect_alignment_summary_metrics_ch
    file "picard/*" from picard_collect_wgs_metrics_ch
    file "picard/*" from picard_collect_insert_size_metrics_ch
    file "picard/*" from picard_collect_gc_bias_metrics_ch
    file "picard/*" from picard_collect_variant_calling_metrics_vcf_ch.collect().ifEmpty([])
    file "picard/*" from picard_collect_variant_calling_metrics_gvcf_ch.collect().ifEmpty([])
    file "verifybamid2/*" from verifybamid2_ch
    file "sg10k_cov_062017/*" from sg10k_cov_062017_ch

    output:
    file "multiqc_data/*" into multiqc_ch

    script:
    """
    fix_multiqc_picard_wgsmetrics.py --input picard/${params.sample_id}.wgs_metrics.txt
    multiqc . --data-format json --enable-npm-plugin
    """

}

process compile_metrics {

    publishDir "${params.publishdir}", mode: "copy"

    input:
    file "*" from multiqc_ch

    output:
    file "*" into compile_metrics_ch

    script:
    """
    compile_metrics.py --multiqc_json multiqc_data.json --output_json ${params.sample_id}.metrics.json
    """

}


/*
----------------------------------------------------------------------
COMPLETION INFO
----------------------------------------------------------------------
*/

workflow.onComplete {
    log.info "Started     : " + workflow.start
    log.info "Completed   : " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Status      : " + workflow.success
    log.info "Publish dir : " + params.publishdir
}

workflow.onError {
    log.info "Workflow execution stopped with the following message:"
    log.info "Exit status   : " + workflow.exitStatus
    log.info "Error message : " + workflow.errorMessage
    log.info "Error report  : " + (workflow.errorReport ?: '-')
}
