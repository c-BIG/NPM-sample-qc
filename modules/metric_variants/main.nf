process metric_variants {

    tag { sample }

    input:
    tuple val(sample), path(variantscount)

    output:
    tuple val(sample), path("${sample}.metrics.json"), emit: metrics

    """
    cp ${variantscount}  "${sample}.metrics.json"
    """
}
