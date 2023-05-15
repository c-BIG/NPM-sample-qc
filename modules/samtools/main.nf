process samtools_stats {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.stats")

    script:
    def reference = ref_fasta ? /--reference "${ref_fasta}"/ : ''

    """
    samtools stats \\
        ${reference} \\
        "${bam}" \\
        > "${sample}.stats"
    """
}
