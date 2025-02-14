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

process getCanu {
    tag "Install Canu"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n canu -y -c bioconda canu
    """
}

process getFlye {
    tag "Install Flye"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n flye -y -c bioconda flye
    """
}

process getRaven {
    tag "Install Raven"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n raven -y -c bioconda raven
    """
}

process getShasta {
    tag "Install Shasta"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n shasta -y -c bioconda shasta
    """
}

process getwtdbg2 {
    tag "Install wtdbg2"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n wtdbg -y -c bioconda wtdbg
    """
}

process getRacon {
    tag "Install Racon"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n racon -y -c bioconda racon
    """
}

process getMummer {
    tag "Install mummer4"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n mummer4 -y -c bioconda mummer4
    """
}

process getQuast {
    tag "Install Quast"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n quast -y -c bioconda quast
    """
}

process getBusco {
    tag "Install BUSCO"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n busco -y -c bioconda -c conda-forge busco
    """
}

}

process getReconciliation {
    tag "Install dependencies required for the reconciliation workflows"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n reconciliation -y -c bioconda quickmerge ragtag
    """
}