process compile_metrics {

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)
    path multiqc

    output:
    tuple val(sample), path("${sample}.metrics.json"), emit: compile_metrics_out

    script:
    """
    # parse and calculate all the metrics in the multiqc output to compile

    compile_metrics.py \
        --multiqc_json multiqc_data.json \
        --output_json "${sample}.metrics.json" \\
        --biosample_id "${sample}" \\
    """
}
