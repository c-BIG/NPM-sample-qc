process mosdepth {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta

    output:
    tuple val(sample), path("*.regions.bed.gz"), emit: regions
    tuple val(sample), path("*.dist.txt"), emit: dists
    tuple val(sample), path("*.summary.txt"), emit: summary

    script:
    def fasta = ref_fasta ? /--fasta "${ref_fasta}"/ : ''

    """
    mosdepth \\
        --no-per-base \\
        --by 1000 \\
        --mapq 20 \\
        --threads ${task.cpus} \\
        ${fasta} \\
        "${sample}" \\
        "${bam}"
    """
}
