process identifymtDNA {
    tag "Identify mtDNA-derived reads with Ganon"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(sample_id), path(fastq)
    val(mitoDNA)

    output:
    path "*.sam", emit: sam

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${mitoDNA} ${fastq} > ${sample_id}_mt.sam
    """
}

process segregate {
    tag "Filter out mtDNA reads"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path('${sample_id}_nuclear.fastq'), emit: nuclearq
    tuple val(sample_id), path('${sample_id}_mt.fastq'), emit: mitoq

    script:
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

  output:
  path "./results/mtGenome"
  path "./results/mtGenome/assembly.fasta", emit: mtContig

  script:
  def fastq = fastqs.join(" ")
  """
  flye -t ${task.cpus} --genome-size 16k --out-dir ./results/mtGenome --nano-raw ${fastq}
  """
}

process mtPolish {
  tag "Polish mitogenome with Racon"
  publishDir "./results/mtGenome/polish", mode: 'copy', overwrite: false, pattern: '**'

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
    publishDir "./results/mtGenome/circular", mode: 'copy', overwrite: false, pattern: '**'

    input:
    val(contig)
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

    output:
    path "./results/mtGenome/annotate"
    path "result.{gff,bed,faa}"
    path "result.fas", emit: mtGenes

    script:
    def sample_id = mtCircular.baseName
    """
    runmitos.py --input ${mtFinal} --outdir ./results/mtGenome/annotate --refseq ${refseq} --code 2 --refdir ${refseqDir}
    """
}

process mtOrtho {
    tag "Construct core mitogenome phylogenetic tree"
    publishDir "./results/mtGenome/phylogenetics", mode: 'copy', overwrite: false, pattern: '**'

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
