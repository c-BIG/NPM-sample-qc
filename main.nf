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
SAMPLE SETTINGS
----------------------------------------------------------------------
*/

/*
----------------------------------------------------------------------
INPUT CHANNELS
----------------------------------------------------------------------
*/

Channel
    .fromPath(params.cram)
    .into { cram_ch_cram_to_bam \
          ; cram_ch_samtools_flagstat }

Channel
    .fromPath(params.vcf)
    .into { vcf_ch_bcftools_stats \
          ; vcf_ch_bcftools_gtcheck \
          ; vcf_ch_picard_collect_variant_calling_metrics }

Channel
    .fromPath(params.ref_fa)
    .into { ref_fa_ch_cram_to_bam \
          ; ref_fa_ch_picard_collect_alignment_summary_metrics \
          ; ref_fa_ch_picard_collect_wgs_metrics \
          ; ref_fa_ch_picard_collect_gc_bias_metrics \
          ; ref_fa_ch_verifybamid2 }

Channel
    .fromPath(params.dbsnp_vcf) \
    .set { dbsnp_vcf_ch_picard_collect_variant_calling_metrics }
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
    file "*" into samtools_stats_ch

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
    file "*" into bcftools_stats_ch

    script:
    """
    bcftools stats ${vcf} > ${params.sample_id}.bcftools_stats.txt
    """

}

process bcftools_gtcheck {

    publishDir "${params.publishdir}/bcftools", mode: "copy"

    input:
    file vcf from vcf_ch_bcftools_gtcheck

    output:
    file "*" into bcftools_gtcheck_ch

    script:
    """
    bcftools index --tbi ${vcf}
    bcftools gtcheck --GTs-only 1 --genotypes ${params.pst_vcf} ${vcf} > ${params.sample_id}.bcftools_gtcheck.txt
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
    picard CollectQualityYieldMetrics I=${params.sample_id}.bam O=${params.sample_id}.quality_yield_metrics.txt
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
    picard CollectAlignmentSummaryMetrics R=${ref_fa} I=${params.sample_id}.bam O=${params.sample_id}.alignment_summary_metrics.txt
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
    picard CollectWgsMetrics R=${ref_fa} I=${params.sample_id}.bam O=${params.sample_id}.wgs_metrics.txt
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
    picard CollectInsertSizeMetrics I=${params.sample_id}.bam O=${params.sample_id}.insert_size_metrics.txt H=${params.sample_id}.insert_size_histogram.pdf M=0.5
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
    picard CollectGcBiasMetrics I=${params.sample_id}.bam O=${params.sample_id}.gc_bias_metrics.txt chart=${params.sample_id}.gc_bias_metrics.pdf S=${params.sample_id}.gc_bias_summary_metrics.txt R=${ref_fa}
    """

}

//process picard_collect_variant_calling_metrics {

    //publishDir "${params.publishdir}/picard", mode: "copy"

    //input:
    //file vcf from vcf_ch_picard_collect_variant_calling_metrics
    //file dbsnp_vcf from dbsnp_vcf_ch_picard_collect_variant_calling_metrics
    
    //output:
    //file "*" into picard_collect_variant_calling_metrics_ch

    //script:
    //"""
    //picard CollectVariantCallingMetrics I=${vcf} O=${params.sample_id}.variant_calling_metrics.txt DBSNP=${dbsnp_vcf}
    //"""

//}

//process verifybamid2 {
    
    //publishDir "${params.publishdir}/verifybamid2", mode: "copy"

    //input:
    //file ref_fa from ref_fa_ch_verifybamid2
    //file "*" from cram_to_bam_ch_verifybamid2

    //output:
    //file "*" into verifybamid2_ch

    //script:
    //"""
    ///home/users/astar/gis/gonzalez/.conda/envs/mgonzalezporta-nscc/share/verifybamid2-1.0.6-0/VerifyBamID --SVDPrefix ${params.vbi2_svdprefix} --Reference ${ref_fa} --BamFile ${params.sample_id}.bam
    //"""

//}

process multiqc {

    publishDir "${params.publishdir}/multiqc", mode: "copy"

    input:
    file "samtools/*stats" from samtools_stats_ch
    file "samtools/*flagstat" from samtools_flagstat_ch
    file "bcftools/*stats" from bcftools_stats_ch
    file "bcftools/*bcftools_gtcheck" from bcftools_gtcheck_ch
    file "picard/*quality_yield_metrics.txt" from picard_collect_quality_yield_metrics_ch
    file "picard/*alignment_summary_metrics.txt" from picard_collect_alignment_summary_metrics_ch
    file "picard/*wgs_metrics.txt" from picard_collect_wgs_metrics_ch
    file "picard/*insert_size_metrics.txt" from picard_collect_insert_size_metrics_ch
    file "picard/*gc_bias_metrics.txt" from picard_collect_gc_bias_metrics_ch
    //file "picard/*variant_calling_metrics.txt" from picard_collect_variant_calling_metrics_ch
    //file "verifybamid2/*" from verifybamid2_ch

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
    file "*" into compile_metrics_ch

    script:
    """
    compile_metrics.py --multiqc_json multiqc_data/multiqc_data.json --output_json ${params.sample_id}.metrics.json
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
