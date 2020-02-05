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
    .into { cram_ch_samtools_stats; cram_ch_samtools_flagstat }

Channel
    .fromPath(params.vcf)
    .set { vcf_ch } 


/*
----------------------------------------------------------------------
PROCESSES
---------------------------------------------------------------------
*/
//process samtools_stats {
//
//    tag "samtools stats"
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

process samtools_flagstat {

    tag "samtools flagstat"

    publishDir "${params.publishdir}/samtools-flagstat", mode: "copy"

    input:
    file cram from cram_ch_samtools_flagstat

    output:
    file "*" into samtools_flagstat_ch

    script:
    """
    samtools flagstat ${cram} > ${cram}.stats
    """

}

process bcftools_stats {

    tag "bcftools stats"

    publishDir "${params.publishdir}/bcftools-stats", mode: "copy"

    input:
    file vcf from vcf_ch

    output:
    file "*" into bcftools_stats_ch

    script:
    """
    bcftools stats ${vcf} > ${vcf}.stats
    """

}

process multiqc {

    tag "multiqc"

    publishDir "${params.publishdir}/multiqc", mode: "copy"

    input:
    file "samtools-flagstat/*" from samtools_flagstat_ch
    file "bcftools-stats/*" from bcftools_stats_ch

    output:
    file "multiqc_data/*" into multiqc_ch

    script:
    """
    multiqc . \
        --data-format json \
        --module samtools \
        --module bcftools
    """

}

process calculate_qc_metrics {

    tag "calculate_qc_metrics"

    publishDir "${params.publishdir}", mode: "copy"

    input:
    file "multiqc_data/*" from multiqc_ch
    
//    output:

    script:
    """
    calculate_qc_metrics.py
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
