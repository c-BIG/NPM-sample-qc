process bcftools_stats {

    tag { sample }

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("${sample}.pass.stats")

    """
    bcftools stats \\
        -f PASS \\
        "${vcf}" \\
        > "${sample}.pass.stats"
    """
}
