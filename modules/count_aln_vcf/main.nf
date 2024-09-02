process count_aln_vcf {

    tag { sample }

    input:
    tuple val(sample), path(count_aln), path(count_variants)

    output:
    tuple val(sample), path("${sample}.metrics.json"), emit: metrics

    """
    # cat ${count_aln} ${count_variants}  > "${sample}.metrics.json"

    compile_aln_variants_metrics.py \
       --input_aln_metrics ${count_aln} \
       --input_variants_metrics ${count_variants} \
       --sample_id ${sample} \
       --output_json ${sample}.metrics.json
    """
}
