process multiqc {
    publishDir '/results/multiqc_reports', mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(raw_reports)

    output:
    path "multiqc", emit: multiqc_report

    script:
    """
    multiqc $raw_reports
    """
}