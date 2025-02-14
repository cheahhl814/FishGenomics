process multiqc {
    publishDir '/results/multiqc_reports', mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(raw_reports)

    output:
    path "multiqc", emit: multiqc_report

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate preassembly
    multiqc $raw_reports
    """
}