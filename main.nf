#!/usr/bin/env nextflow

nextflow.enable.dsl=2
version = "0.5" // nf qc workflow version

/*
Developed by the Genome Institute of Singapore for
the National Precision Medicine Programme

Copyright: 2022 Genome Institute of Singapore
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
                                -profile docker -work-dir ./ --outdir ./
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

def version_message() {
    println("NPM-sample-qc  ~  version ${version}")
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
    this.version_message()
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
PROCESSES
---------------------------------------------------------------------
*/

process samtools_stats {
    tag "${sample_id}"
    publishDir "${params.publishdir}/samtools", mode: "copy"

    input:
    tuple val(sample_id), file(cbam), file(idx), file(fa), file(fai)

    output:
    path "${sample_id}.stats"

    script:
    """
      samtools stats ${cbam} > ${sample_id}.stats
    """

}

process mosdepth {
    tag "${sample_id}"
    publishDir "${params.publishdir}/mosdepth", mode: "copy"

    input:
    tuple val(sample_id), file(cbam), file(idx), file(fa), file(fai)
    path gap_regions

    output:
    path "${sample_id}.*"

    script:
    """
    mosdepth --no-per-base --by 1000 --mapq 20 --threads 4 --fasta ${fa} ${sample_id} ${cbam}
    run_datamash.sh \
        --ref_fasta=${fa} \
        --gap_regions=${gap_regions} \
        --output_csv=${sample_id}.mosdepth.csv \
        --sample_id=${params.sample_id}
    """

}

process picard_collect_multiple_metrics {
    tag "${sample_id}"
    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    tuple val(sample_id), file(cbam), file(idx), file(fa), file(fai)

    output:
    path "${sample_id}.*"

    script:
    """
    picard CollectMultipleMetrics  \
        I=${cbam} \
        O=${sample_id} \
        ASSUME_SORTED=true \
        FILE_EXTENSION=".txt" \
        PROGRAM=null \
        PROGRAM=CollectQualityYieldMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=ALL_READS \
        R=${fa}
    """

}

process multiqc {

    publishDir "${params.publishdir}/multiqc", mode: "copy"
	
    input:
    path "*"

    output:
    path  "multiqc_data/*", emit: multiqc_ch

    script:
    """
    multiqc . --data-format json --enable-npm-plugin
    """

}

process compile_metrics {
    publishDir "${params.publishdir}", mode: "copy"

    input:
    path multiqc

    output:
    path "${params.sample_id}.metrics.json", emit: compile_metrics_out

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
WORKFLOW
---------------------------------------------------------------------
*/
	
// input channels

reference = channel.fromPath(params.reference, checkIfExists: true)
    .map{ fa -> tuple(fa, fa + ".fai") }

input_file = file ( params.bam_cram )
index_type = input_file.getExtension()

if (index_type == "bam")
    cbam = channel.fromPath(params.bam_cram, checkIfExists: true)
        .map{ cbam -> tuple(cbam.simpleName, cbam, cbam + ".bai") }
else if (index_type == "cram")
    cbam = channel.fromPath(params.bam_cram, checkIfExists: true)
        .map{ cbam -> tuple(cbam.simpleName, cbam, cbam + ".crai") }

gap_regions = channel.fromPath(params.gap_regions, checkIfExists: true)

inputs = cbam.combine(reference)

// main
workflow {
    samtools_stats(inputs)
    picard_collect_multiple_metrics(inputs)
    mosdepth( inputs, gap_regions )
    multiqc( samtools_stats.out.mix( picard_collect_multiple_metrics.out, mosdepth.out ).collect() )
    compile_metrics(multiqc.out)
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
