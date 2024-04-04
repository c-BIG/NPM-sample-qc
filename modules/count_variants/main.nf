process count_variants {

    tag { sample }

    input:
    tuple val(sample), path(vcf), path(tbi) 

    output:
    tuple val(sample), path("${sample}.variant_counts.json")

    script:
    """
    echo "test vcf module" ${vcf} ${tbi} >${sample}.variant_counts.json
    # count_variants.py \
    #    --input_vcf ${vcf} \
    #    --output_json ${sample}.variant_counts.json \
    #    --loglevel DEBUG
    """

}