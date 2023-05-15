process multiqc {

    input:
    path 'data/*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    """
    multiqc \\
        --data-format json \\
        --enable-npm-plugin \\
        .
    """
}
