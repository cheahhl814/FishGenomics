process nanoplot {
  tag "Create relevant plots from ONT reads using Nanoplot"
  publishDir "./results/pre-assembly/nanoplot", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)

  output:
  path "*.png", emit: png
  path "*.html", emit: html

  script:
  def fastqs = fastq.join(" ")
  """
  NanoPlot --threads ${task.cpus} --fastq ${fastqs} --maxlength 40000 --tsv_stats --plots dot --format png --info_in_report
  """
}

process porechop {
  tag "Trim adapters"
  publishDir "./results/pre-assembly/porechop", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)

  output:
  path "*_porechop.fastq", emit: porechop_fastq

  script:
  def sample_id = fastq.baseName
  """
  porechop --threads $task.cpus -i $fastq -o ${sample_id}_porechop.fastq --format fastq
  """
}

process filtlong {
    tag "Quality filtering of reads"
    publishDir "./results/pre-assembly/filtlong", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(fastq)

    output:
    path "*_filtlong.fastq", emit: filtlong_fastq

    script:
    def sample_id = fastq.baseName
    """
    filtlong --min_length 1000 --keep_percent 90 ${fastq} > ${sample_id}_filtlong.fastq
    """
}