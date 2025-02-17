process getPreassembly {
    tag "Install dependencies required for the preAssembly workflow"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n preassembly -y
    micromamba activate preassembly
    micromamba install -y -c bioconda -c conda-forge porechop filtlong minimap2 samtools pip seqkit seaborn matplotlib
    pip install nanoplot
    pip install multiqc
    """
}

process getDecon {
    tag "Install Ganon"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n decon -y
    micromamba activate decon
    micromamba install -y -c bioconda ganon seqkit
    """
}

process getCirclator {
    tag "Install Circlator"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n circlator -y
    micromamba activate circlator
    micromamba install -y -c bioconda circlator
    """
}

process getProkka {
    tag "Install Prokka"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n prokka -y
    micromamba activate prokka
    micromamba install -y -c bioconda prokka
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
    micromamba create -n racon -y -c bioconda racon minimap2
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

process getReconciliation {
    tag "Install RagTag and Quickmerge"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba create -n reconciliation -y -c bioconda quickmerge ragtag
    """
}