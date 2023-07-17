process count_variants {

    publishDir "${params.publishdir}/count_variants", mode: "copy"

    input:
    file vcf from vcf_ch_count_variants

    output:
    file "*" into count_variants_ch

    script:
    """
    count_variants.py \
        --input_vcf ${vcf} \
        --output_json ${params.sample_id}.variant_counts.json \
        --loglevel DEBUG
    """

}
