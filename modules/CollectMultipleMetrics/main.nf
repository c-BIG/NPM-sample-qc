process picard_collect_multiple_metrics {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fastq_idx

    output:
    tuple val(sample), path("${sample}.quality_yield_metrics.txt"), emit: quality
    tuple val(sample), path("${sample}.insert_size_metrics.txt"), emit: insert_size

    script:
    def reference = ref_fasta ? /REFERENCE_SEQUENCE="${ref_fasta}"/ : ''
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
    """
}
