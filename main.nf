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
include { picard_collect_variant_calling_metrics_vcf } from './modules/CollectVariantCallingMetrics'
include { bcftools_stats } from './modules/bcftools'
include { count_variants } from './modules/count_variants'
include { count_aln } from './modules/count_aln'
include { count_aln_vcf } from './modules/count_aln_vcf'
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
    autosomes_non_gap_regions_bed = file( params.autosomes_non_gap_regions_bed )
    vbi2_ud = file( params.vbi2_ud )
    vbi2_bed = file( params.vbi2_bed )
    vbi2_mean = file( params.vbi2_mean )
    //ref_dbsnp = file( params.ref_dbsnp )

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


    Channel
        samples.branch { rec ->
            def aln_file = rec.aln ? file( rec.aln ) : null
            def vcf_file = rec.vcf ? file( rec.vcf ) : null

            count: rec.biosample_id && aln_file && vcf_file

                return tuple( rec.biosample_id )
        }
        .view()
        .set { inputs_cond }


    samtools_stats_bam( aln_inputs.bam, [] )
    samtools_stats_cram( aln_inputs.cram, ref_fasta )

    verifybamid2_bam( aln_inputs.bam, ref_fasta, vbi2_ud, vbi2_bed, vbi2_mean )
    verifybamid2_cram( aln_inputs.cram, ref_fasta, vbi2_ud, vbi2_bed, vbi2_mean )

    picard_collect_multiple_metrics_bam( aln_inputs.bam, [], [] )
    picard_collect_multiple_metrics_cram( aln_inputs.cram, ref_fasta, ref_fasta_idx )

    picard_collect_wgs_metrics_bam( aln_inputs.bam, autosomes_non_gap_regions, ref_fasta, ref_fasta_idx )
    picard_collect_wgs_metrics_cram( aln_inputs.cram, autosomes_non_gap_regions, ref_fasta, ref_fasta_idx )


    Channel
        samples.branch { rec ->
            def vcf_file = rec.vcf ? file( rec.vcf ) : null

            output: rec.biosample_id && vcf_file
                def vcf_idx = file( "${rec.vcf}.tbi" )

                return tuple( rec.biosample_id, vcf_file, vcf_idx )
        }
        .set { vcf_inputs }

    bcftools_stats( vcf_inputs )
    //picard_collect_variant_calling_metrics_vcf( vcf_inputs, ref_dbsnp )
    count_variants ( vcf_inputs, autosomes_non_gap_regions_bed )

    Channel
        samples.map { it.biosample_id }
        .set { sample_ids }

// channel for samplelist vcf input file processed outputs
    Channel
        .empty()
        sample_ids
        .join( count_variants.out )
        //count_variants.out
        .join( bcftools_stats.out )
        //.view()
        .set { vcf_qc }

// channel for samplelist input file type bam processed outputs
    Channel
        .empty()
        sample_ids
        .join( samtools_stats_bam.out.metrics )
        .join( picard_collect_multiple_metrics_bam.out.metrics )
        .join( picard_collect_wgs_metrics_bam.out.metrics )
        .join( verifybamid2_bam.out.metrics, remainder: true )
        .set { ch_bam }

// channel for samplelist input file type cram processed outputs
    Channel
        .empty()
        sample_ids
        .join( samtools_stats_cram.out.metrics )
        .join( picard_collect_multiple_metrics_cram.out.metrics )
        .join( picard_collect_wgs_metrics_cram.out.metrics )
        .join( verifybamid2_cram.out.metrics, remainder: true )
        .set { ch_cram }


// channel to mix the bam/cram process outputs and map the verifybamid2 'null' to '[]' if the verifybamid2 process output is empty
    ch_bam.mix(ch_cram)
        .map { sample, stats, quality, wgs_coverage, freemix -> [ sample, stats, quality, wgs_coverage, freemix ?: [] ] }
        .set { aln_count_in }

    count_aln ( aln_count_in )


// channel for samplelist input file type bam processed outputs
    Channel
        .empty()
        sample_ids
        .join( inputs_cond )
        .join( count_aln.out.metrics )
        .join( count_variants.out )
        .view()
        .set { ch_count }

    count_aln_vcf ( ch_count )

/*
// channel to mix the bam/cram process outputs and map the verifybamid2 'null' to '[]' if the verifybamid2 process output is empty
    ch_bam.mix(ch_cram)
        .combine(vcf_qc,by:0)
        //.view()
        .map { sample, stats, insertsize, quality, wgs_coverage, freemix, count_variants, bcftools_stats -> [ sample, stats, insertsize, quality, wgs_coverage, freemix ?: [], count_variants, bcftools_stats ] }
        //.view()
        .set { multiqc_in }
*/

/*
// channel to mix the bam/cram process outputs and map the verifybamid2 'null' to '[]' if the verifybamid2 process output is empty
    ch_bam.mix(ch_cram) // .ifEmpty([])
        //.combine(vcf_qc,by:0)
        .join(vcf_qc, remainder: true)
        //.map { it.minus(null) }
        //.map { files -> files - null }
        //.flatten()
        .map { sample, stats, insertsize, quality, wgs_coverage, freemix, count_variants, bcftools_stats -> [ sample, stats ?: [], insertsize ?: [], quality ?: [], wgs_coverage ?: [], freemix ?: [], count_variants ?: [], bcftools_stats ?: [] ] }
        .view()
        .set { multiqc_in }
*/



/*
// channel to mix the bam/cram process outputs and map the verifybamid2 'null' to '[]' if the verifybamid2 process output is empty
    Channel
        .empty()
        sample_ids
        //.combine(vcf_qc,by:0).ifEmpty([])
        //ch_bam.mix(ch_cram).ifEmpty([])
        .join(count_variants.out)
        .join(bcftools_stats.out)
        .map { sample, count_variants, bcftools_stats -> sample, count_variants, bcftools_stats }
        .mix(ch_bam,ch_cram)
        //.combine(vcf_qc,by:0).ifEmpty([])
        //.join(vcf_qc)
        //.view()
        //.mix(vcf_qc).ifEmpty([])
        //.join(vcf_qc).ifEmpty([])
        //.combine(vcf_qc,by:0).ifEmpty([])
        //.flatten()
        //.toList()
        .view()
        //.map { sample, count_variants, bcftools_stats, sample1, stats, insertsize, quality, wgs_coverage, freemix -> [ sample ?: [], count_variants ?: [], bcftools_stats ?: [], sample1 ?: [], stats ?: [], insertsize ?: [], quality ?: [], wgs_coverage ?: [], freemix ?: [] ] }
        //.view()
        .set { multiqc_in }
*/

/*
    metrics_combine = Channel.empty()
    if (params.vcf && params.aln) {
        metrics_combine =  ( count_variants.out )
    .view()
    }
*/

//    multiqc( multiqc_in )

}

//         .map { sample, stats, insertsize, quality, wgs_coverage, freemix, count_variants, bcftools_stats -> [ sample, stats, insertsize, quality, wgs_coverage, freemix ?: [], count_variants, bcftools_stats ] }

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
