process count_aln {

    tag { sample }

    input:
    tuple val(sample), path(stats), path(picard_quality), path(picard_wgs_coverage), path(verifybamid_freemix)

    output:
    tuple val(sample), path("${sample}.aln.metrics.json"), emit: metrics

    """
    cat ${stats} ${picard_quality} ${picard_wgs_coverage} ${verifybamid_freemix} > "${sample}.aln.metrics.json"
    """
}
