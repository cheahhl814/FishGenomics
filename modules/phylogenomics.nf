process setup1 {
    tag "Create gene FASTA files of closely related species for Orthofinder"
    publishDir "./results/phylogenomics/input", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    path "*_genes_filtered.fasta", emit: geneFasta

    script:
    """
    bedtools -s -name -fi ${fasta} -bed ${gff} -fo ${sample_id}_genes.fasta
    grep '^>' ${sample_id}_genes.fasta | \
        grep -v 'region' | \
        grep -v 'tRNA' | \
        grep -v 'CDS' | \
        grep -v 'exon' | \
        grep -v 'sequence' | \
        grep -v 'ncRNA_gene' | \
        sed 's/^>//' | \
        seqtk subseq ${sample_id}_genes.fasta} - > ${sample_id}_genes_filtered.fasta
    """
}

process setup2 {
    tag "Create gene FASTA files of target species for Orthofinder"
    publishDir "./results/phylogenomics/input", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(gff)
    path(asembly)

    output:
    path "*_mtGenes_filtered.fasta", emit: geneFasta

    script:
    def sample_id = mtFinal.baseName
    """
    bedtools -s -name -fi ${assembly} -bed ${gff} -fo ${sample_id}_genes.fasta
    grep '^>' ${sample_id}_genes.fasta | \
        grep -v 'region' | \
        grep -v 'tRNA' | \
        grep -v 'CDS' | \
        grep -v 'exon' | \
        grep -v 'sequence' | \
        grep -v 'ncRNA_gene' | \
        sed 's/^>//' | \
        seqtk subseq ${sample_id}_genes.fasta} - > ${sample_id}_genes_filtered.fasta
    """
}

process inferOrtho {
    tag "Orthology inference with OrthoFinder"
    publishDir "./results/phylogenomics/tree", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(inputDir)
    val(dummy1)
    val(dummy2)

    output:
    path "SpeciesTreeAlignment.fa", emit: msa

    script:
    """
    orthofinder -t ${task.cpus} -d -M msa -A mafft -oa -f ./results/phylogenomics/input
    """
}

process trimAl {
    tag "Removing gaps in alignment"
    publishDir "./results/phylogenomics/tree", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(msa)

    output:
    path "*_trimal.fasta", emit: trimal_fasta
    path "*_trimal.html", emit: trimal_html

    script:
    def sample_id = msa.baseName
    """
    trimal -in ${msa} -out ${sample_id}_trimal.fasta -htmlout ${sample_id}_trimal.html -nogaps -terminalonly
    """
}

process treeML {
    tag "Building phylogenetic tree"
    publishDir "./results/phylogenomics/tree", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(trimal_fasta)

    output:
    path "*.raxml.*", emit: tree

    script:
    def sample_id = trimal_fasta.baseName
    """
    raxml-ng --all --msa ${trimal_fasta} \
    --model GTR+G+I --bs-metric fbp,tbe --tree pars{25},rand{25} --bs-trees 1000 \
    --prefix ${sample_id} --force
    """
}