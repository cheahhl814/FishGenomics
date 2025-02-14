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
    micromamba create -n canu -y
    micromamba activate canu
    micromamba install -y -c bioconda canu
    micromamba create -n flye -y
    micromamba activate flye
    micromamba install -y -c bioconda flye
    micromamba create -n raven -y
    micromamba activate raven
    micromamba install -y -c bioconda raven
    micromamba create -n shasta -y
    micromamba activate shasta
    micromamba install -y -c bioconda shasta
    micromamba create -n wtdbg -y
    micromamba activate wtdbg
    micromamba install -y -c bioconda wtdbg
    micromamba create -n racon -y
    micromamba activate racon
    micromamba install -y -c bioconda racon
    micromamba create -n mummer4 -y
    micromamba activate mummer4
    micromamba install -y -c bioconda mummer4
    micromamba create -n assessment -y
    micromamba activate assessment
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