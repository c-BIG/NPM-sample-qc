process compile_metrics {

    tag { sample_ids }

    input:
    val(sample_ids)
    path "multiqc_data.json"

    output:
    tuple val(sample_ids), path("${sample_ids}.metrics.json"), emit: compile_metrics_out

    script:
    """
    # parse and calculate all the metrics in the multiqc output to compile

    compile_metrics.py \
        --multiqc_json multiqc_data.json \
        --output_json "${sample_ids}.metrics.json" \\
        --biosample_id "${sample_ids}" \\
    """
}
