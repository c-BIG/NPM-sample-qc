#!/usr/bin/env nextflow

/*
Developed by the Genome Institute of Singapore for
the National Precision Medicine Program 

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
    Usage: nextflow run main.nf -config nextflow.config -params-file sample_params.yaml -profile nscc [-resume] [--keep_workdir] [--help]
    Options:
    -config           Generic workflow settings
    -params-file      Sample-specific settings
    -profile          Job launch settings
    -resume           Re-use existing results (optional, omit to re-run from scratch)
    --keep_workdir    Keep work directory (optional, omit for auto-deletion)
    --help            Print this help message
    """.stripIndent()
}

def minimalInformationMessage() {
    log.info "User name    : " + workflow.userName
    log.info "Command Line : " + workflow.commandLine
    log.info "Profile      : " + workflow.profile
    log.info "Project Dir  : " + workflow.projectDir
    log.info "Launch Dir   : " + workflow.launchDir
    log.info "Work Dir     : " + workflow.workDir
}

def nextflowMessage() {
    log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
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

params.help = null
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
    .into { cram_ch_cram_to_bam; cram_ch_samtools_stats; cram_ch_samtools_flagstat }

Channel
    .fromPath(params.vcf)
    .set { vcf_ch } 

Channel
    .fromPath(params.ref_fa)
    .set { ref_fa_ch_cram_to_bam }

/*
----------------------------------------------------------------------
PROCESSES
---------------------------------------------------------------------
*/
//process samtools_stats {
//
//    publishDir "${params.publishdir}/samtools-stats"
//
//    input:
//    file cram from cram_ch_samtools_stats
//
//    output:
//    file "*" into samtools_stats_ch
//
//    script:
//    """
//    samtools stats ${cram} > ${cram}.stats
//    """
//
//}

process cram_to_bam {

    input:
    file ref_fa from ref_fa_ch_cram_to_bam
    file cram from cram_ch_cram_to_bam

    output:
    file "*" into cram_to_bam_ch

    script:
    """
    samtools view -h -T ${ref_fa} -b ${cram} -o sample.bam
    samtools index sample.bam sample.bam.bai
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
    samtools flagstat ${cram} > ${cram}.flagstat
    """

}

process bcftools_stats {

    publishDir "${params.publishdir}/bcftools", mode: "copy"

    input:
    file vcf from vcf_ch

    output:
    file "*" into bcftools_stats_ch

    script:
    """
    bcftools stats ${vcf} > ${vcf}.stats
    """

}

process picard_collect_quality_yield_metrics {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file "*" from cram_to_bam_ch

    output:
    file "quality_yield_metrics.txt" into picard_collect_quality_yield_metrics_ch

    script:
    """
    picard CollectQualityYieldMetrics I=sample.bam O=quality_yield_metrics.txt
    """

}


process multiqc {

    publishDir "${params.publishdir}/multiqc", mode: "copy"

    input:
    file "samtools/*flagstat" from samtools_flagstat_ch
    file "bcftools/*stats" from bcftools_stats_ch
    file "picard/*quality_yield_metrics.txt" from picard_collect_quality_yield_metrics_ch

    output:
    file "multiqc_data/*" into multiqc_ch

    script:
    """
    multiqc . --data-format json --module samtools --module bcftools
    """

}

process compile_metrics {

    publishDir "${params.publishdir}", mode: "copy"

    input:
    file "multiqc_data/*" from multiqc_ch
    
    output:
    file "metrics.json" into compile_metrics_ch

    script:
    """
    compile_metrics.py --multiqc_json multiqc_data/multiqc_data.json --output_json metrics.json
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
