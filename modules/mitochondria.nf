process segregate {
    tag "Identify mtDNA-derived reads with Ganon"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    val(mitoDNA)
    path(fastq)

    output:
    path "*_mt.fastq", emit: mitoq
    path "*_nuclear.fastq", emit: nuclearq

    script:
    def sample_id = fastq.baseName
    """
    minimap2 -t ${task.cpus} -ax map-ont ${mitoDNA} ${fastq} | samtools view -bS - > ${sample_id}_mt.bam
    samtools fastq -F 4 ${sample_id}_mt.bam > ${sample_id}_mt.fastq
    samtools fastq -f 4 ${sample_id}_mt.bam > ${sample_id}_nuclear.fastq
    """
}

process mtAssembly {
  tag "Assemble mitogenome with Flye"
  publishDir "./results/mtGenome/assembly", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastqs)
  path(flyeDir)

  script:
  def fastq = fastqs.join(" ")
  """
  flye -t ${task.cpus} --genome-size 16k --out-dir ${flyeDir} --nano-raw ${fastq}
  """
}

process mtPolish {
  tag "Polish mitogenome with Racon"
  publishDir "./results/mtGenome/polish", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(contig)
  path(mitoq)

  output:
  path "*_racon.fasta", emit: polished_fasta

  script:
  def fastq = mitoq.join(" ")
  def sample_id = contig.baseName
  """
  cat ${mitoq} > ${sample_id}_mt.fastq
  minimap2 -t ${task.cpus} -ax map-ont ${contig} ${sample_id}_mt.fastq > ${sample_id}_mt.sam
  racon ${sample_id}_mt.fastq ${sample_id}_mt.sam ${contig} > ${sample_id}_racon.fasta
  """
}

process mtCircular {
    tag "Circularize mitogenome with Circlator"
    publishDir "./results/mtGenome/circular", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(contig)
    val(firstGene)

    output:
    path "*_circular.fasta", emit: mtFinal

    script:
    def sample_id = contig.baseName
    """
    circlator fixstart ${contig} --genes_fa ${firstGene} ${sample_id}_circular
    """
}

process mtAnnotate {
    tag "Annotate the circularized mitogenome"
    publishDir "./results/mtGenome/annotate", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(mtFinal)
    val(refseq)
    val(refseqDir)
    path(mitosDir)

    script:
    def sample_id = mtFinal.baseName
    """
    runmitos.py --input ${mtFinal} --outdir ${mitosDir} --refseq ${refseq} --code 2 --refdir ${refseqDir}
    """
}

process orthoSetup {
    tag "Create gene FASTA files of closely related species for Orthofinder"
    publishDir "./results/mtGenome/phylogenetics/input", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    path "*_mtGenes_filtered.fasta", emit: geneFasta

    script:
    """
    bedtools getfasta -s -name -fi ${fasta} -bed ${gff} -fo ${sample_id}_mtGenes.fasta
    grep '^>' ${sample_id}_mtGenes.fasta | grep -v 'region' | grep -v 'tRNA' | grep -v 'CDS' | grep -v 'exon' | grep -v 'sequence' | grep -v 'ncRNA_gene' | sed 's/^>//' | seqtk subseq ${sample_id}_mtGenes.fasta - > ${sample_id}_mtGenes_filtered.fasta
    """
}

process mtOrtho {
    tag "Create gene FASTA files of target species for Orthofinder"
    publishDir "./results/mtGenome/phylogenetics/input", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(mtGenes)
    path(mtFinal)

    output:
    path "*_mtGenes_filtered.fasta", emit: geneFasta

    script:
    def sample_id = mtFinal.baseName
    """
    bedtools getfasta -s -name -fi ${mtFinal} -bed ${mtGenes} -fo ${sample_id}_mtGenes.fasta
    grep '^>' ${sample_id}_mtGenes.fasta | grep -v 'region' | grep -v 'tRNA' | grep -v 'CDS' | grep -v 'exon' | grep -v 'sequence' | grep -v 'ncRNA_gene' | sed 's/^>//' | seqtk subseq ${sample_id}_mtGenes.fasta - > ${sample_id}_mtGenes_filtered.fasta
    """
}

process orthoFinder {
  tag "Orthology inference with OrthoFinder"
  publishDir "./results/mtGenome/phylogenetics/tree", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(inputDir)
  val(dummy1)
  val(dummy2)

  script:
  """
  orthofinder -t ${task.cpus} -d -M msa -A mafft -oa -f ${inputDir}
  """
}

process trimMSA {
  tag "Removing gaps in alignment"
  publishDir "./results/mtGenome/phylogenetics/tree", mode: 'copy', overwrite: false, pattern: '**'

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

process mtTree {
  tag "Building phylogenetic tree"
  publishDir "./results/mtGenome/phylogenetics/tree", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(trimal_fasta)

  output:
  path "*.raxml.*", emit: tree

  script:
  def sample_id = trimal_fasta.baseName
  """
  raxml-ng --all --msa ${trimal_fasta} --model GTR+G+I --bs-metric fbp,tbe --tree pars{25},rand{25} --bs-trees 1000 --prefix ${sample_id} --force
  """
}
