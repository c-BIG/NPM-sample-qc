process samtools_stats {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.stats"), emit: stats
    tuple val(sample), path("${sample}.samtools.metrics"), emit: metrics

    script:
    def reference = ref_fasta ? /--reference "${ref_fasta}"/ : ''

    """
    # get the percentage of reads that have been aligned
    # get the percentage of reads that have been aligned as proper pairs

    samtools stats \\
        ${reference} \\
        "${bam}" \\
        --threads ${task.cpus} \\
        > "${sample}.stats"

    grep "insert size average:" "${sample}.stats" |cut -f3| awk '{print "mean_insert_size\t"\$1}' >"${sample}.samtools.metrics"
    grep "insert size standard deviation:" "${sample}.stats" |cut -f3| awk '{print "insert_size_std_deviation\t"\$1}' >>"${sample}.samtools.metrics"
    grep "percentage of properly paired reads" "${sample}.stats" |cut -f3| awk '{print "pct_reads_properly_paired\t"\$1}' >>"${sample}.samtools.metrics"

    seq=\$(grep "	sequences:" "${sample}.stats" |cut -f3)
    mapped_reads=\$(grep "	reads mapped:" "${sample}.stats" |cut -f3)
    div=\$(awk -v num=\$mapped_reads -v denom=\$seq 'BEGIN { printf ("%.2f", 100 * num / denom) }')
    echo "pct_reads_mapped\t"\$div >>"${sample}.samtools.metrics"
    """
}
