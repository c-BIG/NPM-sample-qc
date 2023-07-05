process samtools_stats {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.stats")

    script:
    def reference = ref_fasta ? /--reference "${ref_fasta}"/ : ''
    def threads = "${params.samtools_threads}" ? /--threads "${params.samtools_threads}"/ : ''

    """
    # get the percentage of reads that have been aligned
    # get the percentage of reads that have been aligned as proper pairs

    samtools stats \\
        ${reference} \\
        "${bam}" \\
        ${threads} \\
        > "${sample}.stats"
    """
}
