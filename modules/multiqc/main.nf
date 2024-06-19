process multiqc {

    tag { sample }

    input:
    // tuple val(sample), path(variantscounts), path(stats)
    tuple val(sample), path(stats), path(picard_insert_size), path(picard_quality), path(picard_wgs_coverage), path(verifybamid_freemix), path(count_variants), path(bcftools_stats)
    // tuple val(sample), path(stats), path(picard_insert_size), path(picard_quality), path(verifybamid_freemix), path(count_variants), path(bcftools_stats)

    output:
    tuple val(sample), path("${sample}/multiqc_report.html"), emit: report
    tuple val(sample), path("${sample}/multiqc_data"), emit: data
    tuple val(sample), path("${sample}/multiqc_data/multiqc_data.json"), emit: json_data
    tuple val(sample), path("${sample}.metrics.json"), emit: compile_metrics_out

    """
    multiqc \\
        --data-format json \\
        --enable-npm-plugin \\
        -o ${sample} \\
        .

    compile_metrics.py \
        --multiqc_json "${sample}/multiqc_data/multiqc_data.json" \
        --output_json "${sample}.metrics.json" \\
        --biosample_id "${sample}" \\
    """
}
