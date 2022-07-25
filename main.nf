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
fcbam = file( params.bam_cram )
fvcf = file( params.vcf )
assert params.bam_cram  != null: "Missing CRAM / BAM input param"
assert params.vcf  != null: "Missing VCF input param"

ftype = fcbam.getExtension()

if ( fcbam.getExtension() == 'cram') {
    println "Input file type and name: CRAM : ${params.cram} \n                               ${params.vcf}"
}
else if ( fcbam.getExtension() == 'bam') {
    println "Input file type and name: BAM : ${params.cram} \n                                ${params.vcf}"
}

Channel
    .fromPath(params.bam_cram)
    .set { input_ch_cram_bam }

Channel
    .fromPath(params.vcf)
    .into { vcf_ch_count_variants
          ; vcf_ch_bcftools_stats
          ; vcf_ch_bcftools_gtcheck
          ; vcf_ch_picard_collect_variant_calling_metrics_vcf
          ; vcf_ch_picard_collect_variant_calling_metrics_gvcf
          }

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
    .into { ref_fa_ch_cram_bam_index
          ; ref_fa_ch_picard_collect_wgs_metrics
          ; ref_fa_ch_picard_collect_multiple_metrics
          ; ref_fa_ch_verifybamid2
          ; ref_fa_ch_mosdepth }

Channel
    .fromPath(params.dbsnp_vcf)
    .into { dbsnp_vcf_ch_picard_collect_variant_calling_metrics_vcf
          ; dbsnp_vcf_ch_picard_collect_variant_calling_metrics_gvcf }

Channel
    .fromPath(params.autosomes_bed)
    .set { autosomes_bed_ch_mosdepth }

Channel
    .fromPath(params.n_regions_bed)
    .set { n_regions_bed_ch_mosdepth }

Channel
    .fromPath(params.vbi2_ud)
    .set { ud_path_ch_verifybamid2 }

Channel
    .fromPath(params.vbi2_bed)
    .set { bed_path_ch_verifybamid2 }

Channel
    .fromPath(params.vbi2_mean)
    .set { mean_path_ch_verifybamid2 }

Channel
    .fromPath(params.version_info)
    .set { version_info_ch }

/*
----------------------------------------------------------------------
PROCESSES
---------------------------------------------------------------------
*/

includeConfig 'conf/base.config'

process cram_bam_index {

    input:
    file ref_fa from ref_fa_ch_cram_bam_index
    file bam_cram from input_ch_cram_bam

    output:
    file "*" into cram_bam_ch_samtools_stats \
                , cram_bam_ch_samtools_flagstat \
                , cram_bam_ch_picard_collect_wgs_metrics \
                , cram_bam_ch_picard_collect_multiple_metrics \
                , cram_bam_ch_verifybamid2 \
                , cram_bam_ch_mosdepth

    script:
        if ( fcbam.getExtension() == 'cram')
            """
            samtools faidx ${ref_fa} -o ${ref_fa}.fai
            mv ${bam_cram} ${params.sample_id}.qc.cram
            samtools index ${params.sample_id}.qc.cram ${params.sample_id}.qc.cram.crai
            """
        else if ( fcbam.getExtension() == 'bam')
            """
            mv ${bam_cram} ${params.sample_id}.qc.bam
            samtools index ${params.sample_id}.qc.bam ${params.sample_id}.qc.bam.bai
            """
}

process samtools_stats {

    publishDir "${params.publishdir}/samtools", mode: "copy"

    input:
    file "*" from cram_bam_ch_samtools_stats

    output:
    file "*" into samtools_stats_ch \
                , samtools_stats_ch_plotbamstats

    script:
    """
      samtools stats ${params.sample_id}.qc.${ftype} > ${params.sample_id}.stats
    """

}

process samtools_flagstat {

    publishDir "${params.publishdir}/samtools", mode: "copy"

    input:
    file bam_cram from cram_bam_ch_samtools_flagstat

    output:
    file "*" into samtools_flagstat_ch

    script:
    """
    samtools flagstat ${params.sample_id}.qc.${ftype} > ${params.sample_id}.flagstat
    """

}

process count_variants {

    publishDir "${params.publishdir}/count_variants", mode: "copy"

    input:
    file vcf from vcf_ch_count_variants

    output:
    file "*" into count_variants_ch

    script:
    """
    count_variants.py \
        --input_vcf ${vcf} \
        --output_json ${params.sample_id}.variant_counts.json \
        --loglevel DEBUG
    """

}

process bcftools_stats {

    publishDir "${params.publishdir}/bcftools", mode: "copy"

    input:
    file vcf from vcf_ch_bcftools_stats

    output:
    file "*" into bcftools_stats_ch

    script:
    """
    bcftools stats -f PASS ${vcf} > ${params.sample_id}.pass.stats
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

process mosdepth {

    publishDir "${params.publishdir}/mosdepth", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_mosdepth
    file "*" from cram_bam_ch_mosdepth
    file autosomes_bed from autosomes_bed_ch_mosdepth
    file n_regions_bed from n_regions_bed_ch_mosdepth

    output:
    file "*" into mosdepth_ch

    script:
    """
    run_mosdepth.sh \
        --input_bam=${params.sample_id}.qc.${ftype} \
        --ref_fasta=${ref_fa} \
        --autosomes_bed=${autosomes_bed} \
        --n_regions_bed=${n_regions_bed} \
        --output_csv=${params.sample_id}.mosdepth.csv \
        --work_dir=.
    """

}

process picard_collect_wgs_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_picard_collect_wgs_metrics
    file "*" from cram_bam_ch_picard_collect_wgs_metrics

    output:
    file "*" into picard_collect_wgs_metrics_ch

    script:
    """
        picard CollectWgsMetrics \
        R=${ref_fa} \
        I=${params.sample_id}.qc.${ftype} \
        O=${params.sample_id}.wgs_metrics.txt
    """

}

process picard_collect_multiple_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file ref_fa from ref_fa_ch_picard_collect_multiple_metrics
    file "*" from cram_bam_ch_picard_collect_multiple_metrics

    output:
    file "*" into picard_collect_multiple_metrics_ch

    script:
    """
    picard CollectMultipleMetrics  \
        I=${params.sample_id}.qc.${ftype} \
        O=${params.sample_id} \
        ASSUME_SORTED=true \
        PROGRAM=null \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectQualityYieldMetrics \
        PROGRAM=CollectGcBiasMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=ALL_READS \
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
    file vbi2_ud from ud_path_ch_verifybamid2
    file vbi2_bed from bed_path_ch_verifybamid2
    file vbi2_mean from mean_path_ch_verifybamid2
    file "*" from cram_bam_ch_verifybamid2

    output:
    file "*" into verifybamid2_ch

    script:
    """
    VerifyBamID --UDPath ${vbi2_ud} --BedPath ${vbi2_bed} --MeanPath ${vbi2_mean} --Reference ${ref_fa} --BamFile ${params.sample_id}.qc.${ftype}
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
    file "count_variants/*" from count_variants_ch
    file "bcftools/*" from bcftools_stats_ch
    file "bcftools/*" from bcftools_gtcheck_ch.collect().ifEmpty([])
    file "picard/*" from picard_collect_wgs_metrics_ch
    file "picard/*" from picard_collect_variant_calling_metrics_vcf_ch.collect().ifEmpty([])
    file "picard/*" from picard_collect_variant_calling_metrics_gvcf_ch.collect().ifEmpty([])
    file "picard/*" from picard_collect_multiple_metrics_ch
    file "verifybamid2/*" from verifybamid2_ch
    file "mosdepth/*" from mosdepth_ch

    output:
    file "multiqc_data/*" into multiqc_ch

    script:
    """
    multiqc . --data-format json --enable-npm-plugin
    """

}

process compile_metrics {

    publishDir "${params.publishdir}", mode: "copy"

    input:
    file "*" from multiqc_ch
    file version_info from version_info_ch

    output:
    file "*" into compile_metrics_ch

    script:
    """
    compile_metrics.py \
        --multiqc_json multiqc_data.json \
        --output_json ${params.sample_id}.metrics.json \
        --version_info ${version_info}
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
