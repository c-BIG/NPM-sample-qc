process multiqc {

    input:
    path 'data/*'

    output:
    path "multiqc_report.html", emit: report
    // path "multiqc_data", emit: data
    path "multiqc_data/multiqc_data.json", emit: json_data

    """
    multiqc \\
        --data-format json \\
        --enable-npm-plugin \\
        .
    """
}
