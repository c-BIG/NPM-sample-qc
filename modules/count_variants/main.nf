process count_variants {

    tag { sample }

    input:
    tuple val(sample), path(vcf), path(tbi)
    path(autosomes_non_gap_regions) 

    output:
    tuple val(sample), path("${sample}.variant_counts.json")

    script:
    """
    count_variants.py \
       --input_vcf ${vcf} \
       --regions ${autosomes_non_gap_regions} \
       --output_json ${sample}.variant_counts.json \
       --loglevel DEBUG
    """

}
