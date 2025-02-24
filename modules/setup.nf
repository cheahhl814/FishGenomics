process pipTools {
    tag "Install additional tools that cannot be installed via Conda"

    script:
    """
    pip install multiqc
    """
}

process mkdir {
    tag "Create directory for certain tools"

    input:
    path(resultDir)

    script:
    """
    mkdir ${resultDir}/mtGenome/assembly
    mkdir ${resultDir}/mtGenome/annotation
    """
}
