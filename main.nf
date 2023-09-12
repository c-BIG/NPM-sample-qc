#!/usr/bin/env nextflow

import groovy.yaml.YamlSlurper

nextflow.enable.dsl=2

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
    Usage: nextflow run main.nf -config nextflow.config -params-file params.yaml 
                                -work-dir ./ --outdir ./
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
    println("NPM-sample-qc  ~  version ${workflow.manifest.version}")
}

def minimalInformationMessage() {
    log.info "User name    : " + workflow.userName
    log.info "Command Line : " + workflow.commandLine
    log.info "Project Dir  : " + workflow.projectDir
    log.info "Launch Dir   : " + workflow.launchDir
    log.info "Work Dir     : " + workflow.workDir
    log.info "Results Dir  : " + params.publish_dir
    log.info "Info Dir     : " + params.info_dir
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

include { samtools_stats as samtools_stats_bam } from './modules/samtools'
include { samtools_stats as samtools_stats_cram } from './modules/samtools'
include { mosdepth as mosdepth_bam } from './modules/mosdepth'
include { mosdepth as mosdepth_cram } from './modules/mosdepth'
include { mosdepth_datamash } from './modules/mosdepth_datamash'
include { verifybamid2 as verifybamid2_bam } from './modules/verifybamid2'
include { verifybamid2 as verifybamid2_cram } from './modules/verifybamid2'
include { picard_collect_multiple_metrics as picard_collect_multiple_metrics_bam } from './modules/CollectMultipleMetrics'
include { picard_collect_multiple_metrics as picard_collect_multiple_metrics_cram } from './modules/CollectMultipleMetrics'
include { multiqc } from './modules/multiqc'
include { compile_metrics } from './modules/compile_metrics'

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/

// main

// params.inputs_list = "inputs.yaml"

workflow {

    ref_fasta = file( params.reference )
    ref_fasta_idx = file( params.reference + ".fai" )
    autosomes_non_gap_regions = file( params.autosomes_non_gap_regions )
    vbi2_ud = file( params.vbi2_ud )
    vbi2_bed = file( params.vbi2_bed )
    vbi2_mean = file( params.vbi2_mean )

    inputs = new YamlSlurper().parse(params.inputs_list as File)

    Channel
        .fromList(inputs['samples'])
        .ifEmpty { ['biosample_id': params.biosample_id, 'aln': params.aln] }
        .set { samples }

    Channel
        samples.branch { rec ->
            def aln_file = rec.aln ? file( rec.aln ) : null

            bam: rec.biosample_id && aln_file?.extension == 'bam'
                def bam_idx = file( "${rec.aln}.bai" )

                return tuple( rec.biosample_id, aln_file, bam_idx )

            cram: rec.biosample_id && aln_file?.extension == 'cram'
                def cram_idx = file( "${rec.aln}.crai" )

                return tuple( rec.biosample_id, aln_file, cram_idx )
        }
        .set { aln_inputs }


    samtools_stats_bam( aln_inputs.bam, [] )
    samtools_stats_cram( aln_inputs.cram, ref_fasta )

    verifybamid2_bam( aln_inputs.bam, ref_fasta, vbi2_ud, vbi2_bed, vbi2_mean )
    verifybamid2_cram( aln_inputs.cram, ref_fasta, vbi2_ud, vbi2_bed, vbi2_mean )

    picard_collect_multiple_metrics_bam( aln_inputs.bam, [], [] )
    picard_collect_multiple_metrics_cram( aln_inputs.cram, ref_fasta, ref_fasta_idx )

    mosdepth_bam( aln_inputs.bam, [] )
    mosdepth_cram( aln_inputs.cram, ref_fasta )

    Channel
        .empty()
        .mix( mosdepth_bam.out.regions )
        .mix( mosdepth_cram.out.regions )
        .set { mosdepth_regions }

    mosdepth_datamash( mosdepth_regions, autosomes_non_gap_regions )
//    mosdepth_datamash( autosomes_non_gap_regions, mosdepth_bam.out.regions.mix( mosdepth_cram.out.regions ) )


    Channel
        .empty()
        .mix( mosdepth_bam.out.dists )
        .mix( mosdepth_bam.out.summary )
        .mix( mosdepth_cram.out.dists )
        .mix( mosdepth_cram.out.summary )
        .mix( mosdepth_datamash.out.coverage )
        .mix( verifybamid2_bam.out.freemix )
        .mix( verifybamid2_cram.out.freemix )
        .mix( verifybamid2_bam.out.ancestry )
        .mix( verifybamid2_cram.out.ancestry )
        .mix( picard_collect_multiple_metrics_bam.out.insert_size )
        .mix( picard_collect_multiple_metrics_cram.out.insert_size )
        .mix( picard_collect_multiple_metrics_bam.out.quality )
        .mix( picard_collect_multiple_metrics_cram.out.quality )
        .mix( samtools_stats_bam.out )
        .mix( samtools_stats_cram.out )
        .map { sample, files -> files }
        .collect()
        .set { log_files }

    multiqc( log_files )


    Channel
        samples.map { it.biosample_id }
        .set { sample_ids }

    compile_metrics ( sample_ids, multiqc.out.json_data )    
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
    log.info "Publish dir : " + params.publish_dir
}

workflow.onError {
    log.info "Workflow execution stopped with the following message:"
    log.info "Exit status   : " + workflow.exitStatus
    log.info "Error message : " + workflow.errorMessage
    log.info "Error report  : " + (workflow.errorReport ?: '-')
}
