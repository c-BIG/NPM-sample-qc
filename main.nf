#!/usr/bin/env nextflow

nextflow.enable.dsl=2
version = "0.6" // nf qc workflow version

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
    tag "${biosample_id}"
    publishDir "${params.publishdir}/samtools", mode: "copy"

    input:
    path cbam
    
    output:
    path "${biosample_id}.stats"

    script:
    """
      samtools stats ${cbam} > ${biosample_id}.stats
    """

}

process mosdepth {
    tag "${biosample_id}"
    publishDir "${params.publishdir}/mosdepth", mode: "copy"

    input:
    path cbam
    path cbam_idx
    path reference
    path reference_idx

    output:
    path "${biosample_id}.*"

    script:
    """
   # run mosdepth    
    mosdepth --no-per-base --by 1000 --mapq 20 --threads 4 --fasta ${reference} ${biosample_id} ${cbam}

    """

}

process mosdepth_datamash {
    tag "${biosample_id}"
    publishDir "${params.publishdir}/mosdepth", mode: "copy"

    input:
    path reference_idx
    path gap_regions
    path mosdepth 

    output:
    path "${biosample_id}.mosdepth.csv"

    script:
    """
    # filter outputs
    # focus on autosomes
    head -22 ${reference_idx} |awk '{print \$1"\t0""\t"\$2}' > autosomes.bed
    zcat "${biosample_id}.regions.bed.gz" | bedtools intersect -a stdin -b autosomes.bed | gzip -9c > "${biosample_id}.regions.autosomes.bed.gz"

    # exclude bins that overlap with N bases in ref
    zcat ${gap_regions} |cut -f2-4 -|egrep -v '_|-|X|Y'|sort -k1,1V -k2,2n > gap_regions.bed
    zcat "${biosample_id}.regions.autosomes.bed.gz" | bedtools intersect -v -a stdin -b gap_regions.bed | gzip -9c > "${biosample_id}.regions.autosomes_minus_n_bases.bed.gz"

    # calculate metrics
    BED="${biosample_id}.regions.autosomes_minus_n_bases.bed.gz";
    mean_coverage=\$(zcat \$BED | datamash --round 6 mean 4);
    sd_coverage=\$(zcat \$BED | datamash --round 6 sstdev 4);
    median_coverage=\$(zcat \$BED | datamash --round 6 median 4);
    mad_coverage=\$(zcat \$BED | datamash --round 6 madraw 4);
    total_bases=\$(zcat \$BED | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_1x_bases=\$(zcat \$BED | awk '\$4>=1' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_10x_bases=\$(zcat \$BED | awk '\$4>=10' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_15x_bases=\$(zcat \$BED | awk '\$4>=15' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_30x_bases=\$(zcat \$BED | awk '\$4>=30' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');
    ge_40x_bases=\$(zcat \$BED | awk '\$4>=40' | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}');

    # save output
    header="mean_autosome_coverage,sd_autosome_coverage,median_autosome_coverage,mad_autosome_coverage,total_autosome_bases,ge_1x_autosome_bases,ge_10x_autosome_bases,ge_15x_autosome_bases,ge_30x_autosome_bases,ge_40x_autosome_bases";
    row="\$mean_coverage,\$sd_coverage,\$median_coverage,\$mad_coverage,\$total_bases,\$ge_1x_bases,\$ge_10x_bases,\$ge_15x_bases,\$ge_30x_bases,\$ge_40x_bases";
    echo "\$header" > "${biosample_id}.mosdepth.csv";
    echo \$row >> "${biosample_id}.mosdepth.csv"
    """

}

process picard_collect_multiple_metrics {
    tag "${biosample_id}"
    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    path cbam
    path cbam_idx
    path reference
    path reference_idx

    output:
    path "${biosample_id}.*"

    script:
    """
    picard CollectMultipleMetrics  \
        I=${cbam} \
        O=${biosample_id} \
        ASSUME_SORTED=true \
        FILE_EXTENSION=".txt" \
        PROGRAM=null \
        PROGRAM=CollectQualityYieldMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=ALL_READS \
        R=${reference}
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
    path "${params.biosample_id}.metrics.json", emit: compile_metrics_out

    script:
    """
    compile_metrics.py \
        --multiqc_json multiqc_data.json \
        --output_json ${params.biosample_id}.metrics.json \
        --biosample_id ${params.biosample_id}
    """

}

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/
	
// input channels

biosample_id = params.biosample_id

reference = channel.fromPath(params.reference, checkIfExists: true)
reference_idx = channel.fromPath(params.reference + ".fai", checkIfExists: true)

input_file = file ( params.bam_cram )
index_type = input_file.getExtension()

if (index_type == "bam") {
    cbam = channel.fromPath(params.bam_cram, checkIfExists: true)
    cbam_idx = channel.fromPath(params.bam_cram + ".bai", checkIfExists: true)
}
else if (index_type == "cram") {
    cbam = channel.fromPath(params.bam_cram, checkIfExists: true)
    cbam_idx = channel.fromPath(params.bam_cram + ".crai", checkIfExists: true)
}

gap_regions = channel.fromPath(params.gap_regions, checkIfExists: true)

// main
workflow {
    samtools_stats(cbam)
    picard_collect_multiple_metrics(cbam, cbam_idx, reference, reference_idx)
    mosdepth( cbam, cbam_idx, reference, reference_idx )
    mosdepth_datamash( reference_idx, gap_regions, mosdepth.out )
    multiqc( samtools_stats.out.mix( picard_collect_multiple_metrics.out, mosdepth.out, mosdepth_datamash.out ).collect() )
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
