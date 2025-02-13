process getPreassembly {
    tag "Install dependencies required for the preAssembly workflow"

    script:
    """
    micromamba create -y -n preassembly porechop filtlong minimap2 samtools pip seqkit
    pip install nanoplot
    """
}

process getAssembly {
    tag "Install dependencies required for the assembly workflows"

    script:
    """
    micromamba create -y -n assembly canu flye raven shasta wtdbg racon busco quast mummer4
    """
}

process getReconciliation {
    tag "Install dependencies required for the reconciliation workflows"

    script:
    """
    micromamba create -y -n reconciliation quickmerge ragtag
    """
}