process identifymtDNA {
    tag "Identify mtDNA-derived reads with Ganon"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(fastq)
    path(mitoDNA)

    output:
    path "*.sam", emit: sam

    script:
    def database_id = mitoDNA.baseName
    def sample_id = fastq.baseName
    """
    minimap2 -t ${task.cpus} -ax map-ont ${mitoDNA} ${fastq} > ${sample_id}_mt.sam
    """
}

process segregateReads {
    tag "Filter out mtDNA reads"\
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(sam)
    path{fastq}

    output:
    path "*_nuclear.fastq", emit: nuclearq
    path "*_mt.fastq", emit: mitoq

    script:
    def sample_id = sam.baseName
    """
    samtools view -b -f 4 ${sam} | samtools fastq - > ${sample_id}_nuclear.fastq
    samtools view -b -F 4 ${sam} | samtools fastq - > ${sample_id}_mt.fastq
    """
}

process mtAssembly {
  tag "Assemble mitogenome with Flye"
  publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(fastqs)
  path(flyeDir)

  output:
  path "*.fasta", emit: mtContig

  script:
  def fastq = fastqs.join(" ")
  """
  flye -t ${task.cpus} --genome-size 16k --out-dir $flyeDir --nano-raw ${fastq}
  """
}

process mtPolish {
  tag "Polish mitogenome with Racon"
  publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contig)
  val(mitoq)

  output:
  path "*.paf", emit: minimap_output
  path "*_racon.fasta", emit: polished_fasta

  script:
  def fastq = mitoq.join(" ")
  def sample_id = contig.baseName
  """
  minimap2 ${contig} ${fastq} -o ${sample_id}.paf
  racon ${fastq} ${sample_id}.paf ${contig} > ${sample_id}_racon.fasta
  """
}

process mtCircular {
    tag "Circularize mitogenome with Circlator"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(contig)
    path(firstGene)

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
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(mtFinal)

    output:
    path "*.gff", emit: mtGFF
    path "*.gbk", emit: mtGenBank
    path "*.faa", emit: mtProtein
    path "*.ffn", emit: mtProteinN

    script:
    def sample_id = mtCircular.baseName
    """
    prokka --prefix ${sample_id} --compliant --addgenes --kingdom Mitochondria ${mtCircular}
    """
}

process mtOrtho {
    tag "Construct core mitogenome phylogenetic tree"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(mtProteinN)
    path(treeDir)

    output:
    path "SpeciesTreeAlignment.fa", emit: msa

    script:
    """
    cp ${mtProteinN} ${treeDir}
    orthofinder -t ${task.cpus} -d -M msa -A mafft -oa -f ${treeDir}
    """
}

process trimMSA {
  tag "Removing gaps in alignment"
  publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

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
  publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

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
