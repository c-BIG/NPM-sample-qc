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
include { verifybamid2 as verifybamid2_bam } from './modules/verifybamid2'
include { verifybamid2 as verifybamid2_cram } from './modules/verifybamid2'
include { picard_collect_multiple_metrics as picard_collect_multiple_metrics_bam } from './modules/CollectMultipleMetrics'
include { picard_collect_multiple_metrics as picard_collect_multiple_metrics_cram } from './modules/CollectMultipleMetrics'
include { picard_collect_wgs_metrics as picard_collect_wgs_metrics_bam } from './modules/CollectWgsMetrics'
include { picard_collect_wgs_metrics as picard_collect_wgs_metrics_cram } from './modules/CollectWgsMetrics'
include { multiqc } from './modules/multiqc'

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/

// main

workflow {

    ref_fasta = file( params.reference )
    ref_fasta_idx = file( params.reference + ".fai" )
    autosomes_non_gap_regions = file( params.autosomes_non_gap_regions )
    vbi2_ud = file( params.vbi2_ud )
    vbi2_bed = file( params.vbi2_bed )
    vbi2_mean = file( params.vbi2_mean )

    inputs = new YamlSlurper().parse(file(params.inputs_list))

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

    picard_collect_wgs_metrics_bam( aln_inputs.bam, autosomes_non_gap_regions, ref_fasta, ref_fasta_idx )
    picard_collect_wgs_metrics_cram( aln_inputs.cram, autosomes_non_gap_regions, ref_fasta, ref_fasta_idx )

// channel for samplelist input file type bam processed outputs
    Channel
        .empty()
        samtools_stats_bam.out.stats
        .join( picard_collect_multiple_metrics_bam.out.insert_size )
        .join( picard_collect_multiple_metrics_bam.out.quality )
        .join( picard_collect_wgs_metrics_bam.out.wgs_coverage )
        .join( verifybamid2_bam.out.freemix, remainder: true )
        .set { ch_bam }

// channel for samplelist input file type cram processed outputs
    Channel
        .empty()
        samtools_stats_cram.out.stats
        .join( picard_collect_multiple_metrics_cram.out.insert_size )
        .join( picard_collect_multiple_metrics_cram.out.quality )
        .join( picard_collect_wgs_metrics_cram.out.wgs_coverage )
        .join( verifybamid2_cram.out.freemix, remainder: true )
        .set { ch_cram }

// channel to mix the bam/cram process outputs and map the verifybamid2 'null' to '[]' if the verifybamid2 process output is empty
    ch_bam.mix(ch_cram)
        .map { sample, stats, insertsize, quality, wgs_coverage, freemix -> [ sample, stats, insertsize, quality, wgs_coverage, freemix ?: [] ] }
        .set { multiqc_in }

    multiqc( multiqc_in )

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
