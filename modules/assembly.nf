process canu {
  tag "Assemble genome with Canu"
  publishDir "./results/assembly/canu", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)
  val(genomeSize)

  output:
  path "*.contigs.fasta", emit: canu_contig

  script:
  def sample_id = fastq.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate canu
  canu -p $sample_id genomeSize=$genomeSize -nanopore $fastq
  """
}

process wtdbg2 {
  tag "Assemble genome with wtdbg2"
  publishDir "./results/assembly/wtdbg2", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)
  val(genomeSize)

  output:
  path "*.ctg.fa", emit: wtdbg2_contig

  script:
  def sample_id = fastq.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate wtdbg
  wtdbg2 -t ${task.cpus} -x ont -g $genomeSize -fo $baseName -x ont -i $fastq
  wtpoa-cns -t ${task.cpus} -i ${baseName}.ctg.lay.gz -fo ${baseName}.ctg.fa
  """
}

process flye {
  tag "Assemble genome with Flye"
  publishDir "./results/assembly/flye", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)
  val(genomeSize)
  path(flyeDir)

  output:
  path "*.fasta", emit: flye_contig

  script:
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate flye
  flye -t ${task.cpus} --genome-size $genomeSize --out-dir $flyeDir --nano-raw $fastq
  """
}

process raven {
  tag "Assemble genome with Raven"
  publishDir "./results/assembly/raven", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)

  output:
  path "*.fasta.gz", emit: raven_contig
  path "*.gfa.gz", emit: raven_gfa 

  script:
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate raven
  raven -t ${task.cpus} $fastq
  """
}


process shasta {
  tag "Assemble genome with Shasta"
  publishDir "./results/assembly/shasta", mode: 'copy', overwrite: false

  input:
  path(fastq)

  output:
  path "Assembly.fasta", emit: shasta_contig

  script:
  def sample_id = fastq.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate shasta
  FastqToFasta.py $fastq ${sample_id}.fasta
  shasta --input ${sample_id}.fasta
  """
}

process racon {
  tag "Polish genome with Racon"
  publishDir "./results/assembly/racon", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contig)
  path(fastq)

  output:
  path "*.paf", emit: minimap_output
  path "*_racon.fasta", emit: polished_fasta

  script:
  def sample_id = contig.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate preassembly
  minimap2 $contig $fastq -o ${sample_id}.paf
  micromamba activate racon
  racon $fastq ${sample_id}.paf $contig > ${sample_id}_racon.fasta
  """
}
