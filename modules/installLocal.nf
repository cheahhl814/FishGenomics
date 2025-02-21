process pipTools {
    tag "Install additional tools that cannot be installed via Conda"

    script:
    """
    pip install multiqc
    pip install nanoplot
    """
}
