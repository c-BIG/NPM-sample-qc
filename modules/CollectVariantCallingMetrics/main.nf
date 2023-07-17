process picard_collect_variant_calling_metrics_vcf {

    publishDir "${params.publishdir}/picard", mode: "copy"

    input:
    file vcf from vcf_ch_picard_collect_variant_calling_metrics_vcf
    file dbsnp_vcf from dbsnp_vcf_ch_picard_collect_variant_calling_metrics_vcf

    output:
    file "*metrics" into picard_collect_variant_calling_metrics_vcf_ch

    when:
    vcf.name =~ /.*\.vcf.gz$/

    script:
    """
    bcftools index --tbi ${vcf}
    picard CollectVariantCallingMetrics \
        I=${vcf} \
        O=${params.sample_id} \
        DBSNP=${dbsnp_vcf}
    """

}
