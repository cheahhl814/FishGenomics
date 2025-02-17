process identifymtDNA {
    tag "Identify mtDNA-derived reads with Ganon"
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(fastq)
    path(mitoDNA)

    output:
    path "*.tre", emit: tre
    path "*.rep", emit: rep
    path "*.one", emit: one

    script:
    def database_id = mitoDNA.baseName
    """
    ganon build-custom --input-target sequence --level sequence --taxonomy skip --input ${mitoDNA} --db-prefix ${database_id} --threads ${task.cpus}
    ganon classify --db-prefix ${database_id} --single-reads ${fastq} --output-one --output-prefix mitoDNA --threads ${task.cpus}
    """
}

process segregateReads {
    tag "Filter out mtDNA reads"\
    publishDir "./results/mtGenome", mode: 'copy', overwrite: false, pattern: '**'

    input:
    val(ones)
    path{fastq}

    output:
    path "*_nuclear.fastq", emit: nuclearq
    path "*_nuclear.fastq", emit: mitoq

    script:
    def one = ones.join{" "}
    def sample_id = fastq.baseName
    """
    cat ${one} | cut -f 1 | sort | uniq | seqkit grep -v -f - ${fastq} > ${sample_id}_nuclear.fastq
    cat ${one} | cut -f 1 | sort | uniq | seqkit grep -f - ${fastq} > ${sample_id}_mt.fastq
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

process mtTree {
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