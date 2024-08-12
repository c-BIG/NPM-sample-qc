process count_variants {

    tag { sample }

    input:
    tuple val(sample), path(vcf), path(tbi)
    path(autosomes_non_gap_regions) 

    output:
    tuple val(sample), path("${sample}.variant.metrics.json")

    script:
    """
    count_variants.py \
       --input_vcf ${vcf} \
       --regions ${autosomes_non_gap_regions} \
       --output_json ${sample}.variant.metrics.json \
       --loglevel DEBUG
    """

}
