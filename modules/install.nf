process getPreassembly {
    tag "Install dependencies required for the preAssembly workflow"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n preassembly -y
    micromamba activate preassembly
    micromamba install -y -c bioconda -c conda-forge porechop filtlong minimap2 samtools pip seqkit
    pip install nanoplot
    """
}

process getAssembly {
    tag "Install dependencies required for the assembly workflows"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n assembly -y
    micromamba activate assembly
    micromamba install -y -c bioconda -c conda-forge canu flye raven shasta wtdbg racon mummer4
    micromamba create -n assessment -y
    micromamba install -y -c bioconda -c conda-forge busco quast 
    """
}

process getReconciliation {
    tag "Install dependencies required for the reconciliation workflows"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n reconciliation -y
    micromamba activate reconciliation
    micromamba install -c bioconda quickmerge ragtag
    """
}