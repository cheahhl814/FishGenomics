process canu {
  tag "Assemble genome with Canu"
  publishDir "./results/assembly/canu", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(fastqs)
  val(genomeSize)
  val(name)

  output:
  path "*.contigs.fasta", emit: canu_contig

  script:
  def fastq = fastqs.join(" ")
  """
  canu -p ${name} genomeSize=${genomeSize} -nanopore ${fastq}
  """
}

process wtdbg2 {
  tag "Assemble genome with wtdbg2"
  publishDir "./results/assembly/wtdbg2", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(fastqs)
  val(genomeSize)
  val(name)

  output:
  path "*.ctg.fa", emit: wtdbg2_contig

  script:
  def fastq = fastqs.join(" ")
  """
  wtdbg2 -t ${task.cpus} -x ont -g ${genomeSize} -fo ${name}_wtdbg -x ont -i ${fastq}
  wtpoa-cns -t ${task.cpus} -i ${name}_wtdbg.ctg.lay.gz -fo ${name}_wtdbg.ctg.fa
  """
}

process flye {
  tag "Assemble genome with Flye"
  publishDir "./results/assembly/flye", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(fastqs)
  val(genomeSize)

  output:
  path "./results/assembly/flye"
  path "./results/assembly/assembly.fasta", emit: flye_contig

  script:
  def fastq = fastqs.join(" ")
  """
  flye -t ${task.cpus} --genome-size ${genomeSize} --out-dir ./results/assembly/flye --nano-raw ${fastq}
  """
}

process raven {
  tag "Assemble genome with Raven"
  publishDir "./results/assembly/raven", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(fastqs)

  output:
  path "*.fasta.gz", emit: raven_contig
  path "*.gfa.gz", emit: raven_gfa

  script:
  def fastq = fastqs.join(" ")
  """
  raven -t ${task.cpus} ${fastq}
  """
}

process shasta {
  tag "Assemble genome with Shasta"
  publishDir "./results/assembly/shasta", mode: 'copy', overwrite: false

  input:
  val(fastqs)
  val(name)

  output:
  path "Assembly.fasta", emit: shasta_contig

  script:
  def fastq = fastqs.join(" ")
  """
  FastqToFasta.py ${fastq} ${name}.fasta
  shasta --input ${name}.fasta
  """
}

process racon {
  tag "Polish genome with Racon"
  publishDir "./results/assembly/racon", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contig)
  val(fastqs)

  output:
  path "*.paf", emit: minimap_output
  path "*_racon.fasta", emit: polished_fasta

  script:
  def fastq = fastqs.join(" ")
  def sample_id = contig.baseName
  """
  minimap2 ${contig} ${fastq} -o ${sample_id}.paf
  racon ${fastq} ${sample_id}.paf ${contig} > ${sample_id}_racon.fasta
  """
}
