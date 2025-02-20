process getMultiqc {
    tag "Install additional tools that cannot be installed via Conda"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate preassembly
    pip install multiqc
    """
}