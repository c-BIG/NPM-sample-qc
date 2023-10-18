process multiqc {

    tag { sample }

    input:
    tuple val(sample), path('*')

    output:
    
    tuple val(sample), path("multiqc_report.html"), emit: report
    tuple val(sample), path("multiqc_data"), emit: data
    tuple val(sample), path("${sample}/multiqc_data/multiqc_data.json"), emit: json_data

    """
    multiqc \\
        --data-format json \\
        --enable-npm-plugin \\
        .
    """
}
