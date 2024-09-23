process metric_aln {

    tag { sample }

    input:
    tuple val(sample), path(alncount)

    output:
    tuple val(sample), path("${sample}.metrics.json"), emit: metrics

    """
    cp ${alncount}  "${sample}.metrics.json"
    """
}
