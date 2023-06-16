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
   # run mosdepth    
   # do not output per-base depth and apply 1000bp window-sizes
   # use the reads with mapping quality 20 and above

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
