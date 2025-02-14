process quast {
  tag "Assess the genome with QUAST"
  publishDir "./results/assembly/quast", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contig)
  path(reference)

  output:
  path "report.txt", emit: summary
  path "report.tsv", emit: summary_table
  path "report.tex", emit: latex
  path "report.pdf", emit: pdf
  path "report.html", emit: html
  path "icarus.html", emit: icarus

  script:
  def sample_id = contig.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate assessment
  quast.py -t ${task.cpus} -l $sample_id $contig -r $reference
  """
}

process busco {
  tag "Assess the genome with BUSCO"
  publishDir "./results/assembly/busco", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(genome)

  output:
  path "*.gff", emit: gff
  path "full_table.tsv", emit: full_table
  path "missing_busco_list.tsv", emit: missing_busco
  path "short_summary.*.txt", emit: busco_summary_txt
  path "short_summary.*.json", emit: busco_summary_json

  script:
  def sample_id = genome.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate assessment
  busco -c ${task.cpus} -i $genome --auto-lineage-euk -o $sample_id -m genome
  """
}

process galignment {
  tag "Whole genome alignment with the reference genome"
  publishDir "./results/assembly/genomeAlignment", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(scaffold)
  path(reference)

  output:
  path "*.delta", emit:delta
  path "*.report", emit: report
  path "*.snps", emit: snps

  script:
  def sample_id = scaffold.baseName
  """
  eval "\$(micromamba shell hook --shell bash)"
  micromamba activate mummer4
  dnadiff -p $sample_id $reference $scaffold
  mummerplot -p $sample_id -t $sample_id --png $delta
  """
}
