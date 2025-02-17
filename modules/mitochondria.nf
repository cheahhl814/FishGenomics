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
    output:

    script:
    """
    circlator fixstart <assembly.fasta> <output_prefix>
    """
}