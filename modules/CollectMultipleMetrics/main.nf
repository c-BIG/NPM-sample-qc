process picard_collect_multiple_metrics {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fastq_idx

    output:
    tuple val(sample), path("${sample}.quality_yield_metrics.txt"), emit: quality
    tuple val(sample), path("${sample}.insert_size_metrics.txt"), emit: insert_size
    tuple val(sample), path("${sample}.quality_yield_metrics.metrics"), emit: metrics

    script:
    def reference = ref_fasta ? /R="${ref_fasta}"/ : ''
    """
    # program CollectQualityYieldMetrics to get numbers of bases that pass a base quality score 30 threshold
    # program CollectInsertSizeMetrics to get mean insert size

    picard CollectMultipleMetrics  \
        I=${bam} \
        O=${sample} \
        ${reference} \
        ASSUME_SORTED=true \
        FILE_EXTENSION=".txt" \
        PROGRAM=null \
        PROGRAM=CollectQualityYieldMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        METRIC_ACCUMULATION_LEVEL=null \
        METRIC_ACCUMULATION_LEVEL=ALL_READS 

    cut -f9 "${sample}.quality_yield_metrics.txt"  | grep "PF_Q30_BASES" -A 1 | grep -v "PF_Q30_BASES" | awk '{print "yield_bp_q30\t"\$1}' > "${sample}.quality_yield_metrics.metrics"
    """
}
