process picard_collect_wgs_metrics {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path autosomes_non_gap_regions
    path ref_fasta
    path ref_fastq_idx

    output:
    tuple val(sample), path("${sample}_wgs_metrics.txt"), emit: wgs_coverage

    script:
    """
    # program CollectWgsMetrics to get mean_autosome_coverage, mad_autosome_coverage and pct_autosomes_15x

    picard CollectWgsMetrics  \
        I=${bam} \
        O=${sample}_wgs_metrics.txt \
        R=${ref_fasta} \
        INTERVALS=${autosomes_non_gap_regions} \
        VALIDATION_STRINGENCY=SILENT
    """
}
