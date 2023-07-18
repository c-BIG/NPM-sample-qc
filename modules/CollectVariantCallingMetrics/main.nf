process picard_collect_variant_calling_metrics_vcf {
    tag { sample }

    input:
    tuple val(sample), path(vcf), path(tbi)
    path(ref_dbsnp) 

    output:
    tuple val(sample), path("${sample}.metrics")

    script:
    """
    picard CollectVariantCallingMetrics \
        I=${vcf} \
        O=${sample} \
        DBSNP=${ref_dbsnp}
    """

}
