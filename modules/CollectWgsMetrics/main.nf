process picard_collect_wgs_metrics {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path autosomes_non_gap_regions
    path ref_fasta
    path ref_fastq_idx

    output:
    tuple val(sample), path("${sample}_wgs_metrics.txt"), emit: wgs_coverage
    tuple val(sample), path("${sample}_wgs_metrics.metrics"), emit: metrics

    script:
    """
    # program CollectWgsMetrics to get mean_autosome_coverage, mad_autosome_coverage and pct_autosomes_15x

    picard CollectWgsMetrics  \
        I=${bam} \
        O=${sample}_wgs_metrics.txt \
        R=${ref_fasta} \
        INTERVALS=${autosomes_non_gap_regions} \
        VALIDATION_STRINGENCY=SILENT

    #cut -f2 "${sample}_wgs_metrics.txt"  | grep "MEAN_COVERAGE" -A 1 | grep -v "MEAN_COVERAGE" | awk '{print "mean_autosome_coverage\t"\$1}' > "${sample}_wgs_metrics.metrics"
    avg=\$(cut -f2 "${sample}_wgs_metrics.txt"  | grep "MEAN_COVERAGE" -A 1 | tail -n 1)
    avgcov=\$(awk -v avgval=\$avg 'BEGIN { printf ("%.2f", avgval) }')
    echo "mean_autosome_coverage\t"\$avgcov >> "${sample}_wgs_metrics.metrics"

    cov=\$(cut -f17 NA12878-chr14-AKT1_wgs_metrics.txt | grep "PCT_15X" -A 1 | tail -n 1)
    covpct=\$(awk -v covval=\$cov 'BEGIN { printf ("%.2f", 100 * covval) }')
    echo "pct_autosomes_15x\t"\$covpct >> "${sample}_wgs_metrics.metrics"

    cut -f5 "${sample}_wgs_metrics.txt"  | grep "MAD_COVERAGE"  -A 1 | grep -v "MAD_COVERAGE"  | awk '{print "mad_autosome_coverage\t"\$1}'  >> "${sample}_wgs_metrics.metrics"
    """
}
