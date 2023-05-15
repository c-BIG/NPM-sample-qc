#!/usr/bin/env nextflow

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

include { bcftools_stats } from './modules/bcftools'
include { mosdepth as mosdepth_bam } from './modules/mosdepth'
include { mosdepth as mosdepth_cram } from './modules/mosdepth'
include { multiqc } from './modules/multiqc'
include { compile_metrics } from './modules/compile_metrics'
include { samtools_stats as samtools_stats_bam } from './modules/samtools'
include { samtools_stats as samtools_stats_cram } from './modules/samtools'

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/

// main

workflow {

    ref_fasta = file( params.reference )

    Channel
        .fromList( params.samples )
        .branch { rec ->
            def aln_file = file( rec.bam )

            bam: aln_file.extension == 'bam'
                def bam_idx = file( "${rec.bam}.bai" )

                return tuple( rec.id, aln_file, bam_idx )

            cram: aln_file.extension == 'cram'
                def cram_idx = file( "${rec.bam}.crai" )

                return tuple( rec.id, aln_file, cram_idx )
        }
        .set { aln_inputs }

    Channel
        .fromList( params.samples )
        .map { rec ->
            def vcf_file = file( rec.vcf )
            def vcf_idx = file( "${rec.vcf}.tbi" )

            tuple( rec.id, vcf_file, vcf_idx )
        }
        .set { vcf_inputs }


    Channel
        .fromList( params.samples )
        .branch { rec ->
            def aln_file = file( rec.bam )

            bam: aln_file.extension == 'bam'
                def bam_idx = file( "${rec.bam}.bai" )

                return tuple( rec.id, aln_file, bam_idx )

            cram: aln_file.extension == 'cram'
                def cram_idx = file( "${rec.bam}.crai" )

                return tuple( rec.id, aln_file, cram_idx )
        }
        .set { compile_metrics_inputs }


    mosdepth_bam( aln_inputs.bam, [] )
    mosdepth_cram( aln_inputs.cram, ref_fasta )

    samtools_stats_bam( aln_inputs.bam, [] )
    samtools_stats_cram( aln_inputs.cram, ref_fasta )

    bcftools_stats( vcf_inputs )

    Channel
        .empty()
        .mix( mosdepth_bam.out.dists )
        .mix( mosdepth_bam.out.summary )
        .mix( mosdepth_cram.out.dists )
        .mix( mosdepth_cram.out.summary )
        .mix( samtools_stats_bam.out )
        .mix( samtools_stats_cram.out )
        .map { sample, files -> files }
        .collect()
        .set { log_files }

    multiqc( log_files )

   compile_metrics ( compile_metrics_inputs.mix(multiqc.out).collect() )    
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
