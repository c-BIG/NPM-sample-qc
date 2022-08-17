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
assert params.bam_cram  != null: "Missing CRAM / BAM input param"

ftype = fcbam.getExtension()

if ( fcbam.getExtension() == 'cram') {
    println "Input file type and name: CRAM : ${params.bam_cram}"
}
else if ( fcbam.getExtension() == 'bam') {
    println "Input file type and name: BAM : ${params.bam_cram}"
}

Channel
    .fromPath(params.bam_cram)
    .set { input_ch_cram_bam }

Channel
    .fromPath(params.ref_fa)
    .into { ref_fa_ch_cram_bam_index
          ; ref_fa_ch_picard_collect_multiple_metrics
          ; ref_fa_ch_mosdepth }

Channel
    .fromPath(params.autosomes_bed)
    .set { autosomes_bed_ch_mosdepth }

Channel
    .fromPath(params.n_regions_bed)
    .set { n_regions_bed_ch_mosdepth }


/*
----------------------------------------------------------------------
PROCESSES
---------------------------------------------------------------------
*/

process cram_bam_index {

    input:
    file ref_fa from ref_fa_ch_cram_bam_index
    file bam_cram from input_ch_cram_bam

    output:
    file "*" into cram_bam_ch_samtools_stats \
                , cram_bam_ch_picard_collect_multiple_metrics \
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
    file "*" into samtools_stats_ch

    script:
    """
      samtools stats ${params.sample_id}.qc.${ftype} > ${params.sample_id}.stats
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
        FILE_EXTENSION=".txt" \
        PROGRAM=null \
        PROGRAM=CollectQualityYieldMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=ALL_READS \
        R=${ref_fa}
    """

}

process multiqc {

    publishDir "${params.publishdir}/multiqc", mode: "copy"

    input:
    file "samtools/*" from samtools_stats_ch
    file "picard/*" from picard_collect_multiple_metrics_ch
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

    output:
    file "*" into compile_metrics_ch

    script:
    """
    compile_metrics.py \
        --multiqc_json multiqc_data.json \
        --output_json ${params.sample_id}.metrics.json \
        --sample_id ${params.sample_id}
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
