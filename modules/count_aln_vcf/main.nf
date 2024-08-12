process count_aln_vcf {

    tag { sample }

    input:
    tuple val(sample), path(count_aln), path(count_variants)

    output:
    tuple val(sample), path("${sample}.aln.variant.metrics.json"), emit: metrics

    """
    cat ${count_aln} ${count_variants}  > "${sample}.aln.variant.metrics.json"
    """
}
