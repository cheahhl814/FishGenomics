process sequali {
  tag "Quality checking of ONT reads using Sequali"
  publishDir "./results/pre-assembly/sequali", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)
  val(sample_id)

  output:
  path "*.html", emit: html
  path "*json", emit: json

  script:
  def fastqs = fastq.join(" ")
  """
  cat ${fastqs} > ${sample_id}.fastq
  """
}

process nanoplot {
  tag "Create relevant plots from ONT reads using Nanoplot"
  publishDir "./results/pre-assembly/nanoplot", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fastq)
  val(sample_id)

  output:
  path "${sample_id}.fastq", emit: mergedFastq
  path "*.png", emit: png
  path "*.html", emit: html

  script:
  def fastqs = fastq.join(" ")
  """
  cat ${fastqs} > ${sample_id}.fastq
  NanoPlot --threads ${task.cpus} --fastq ${sample_id}.fastq --maxlength 40000 --tsv_stats --plots dot --format png --info_in_report
  """
}

process porechop {
    tag "Trim adapters from Nanopore reads"
    publishDir "./results/pre-assembly/porechop", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(fastq)

    output:
    path "*_porechop.fastq", emit: porechop_fastq

    script:
    """
    porechop --threads ${task.cpus} -i ${fastq} -o ${baseName}_porechop.fastq --format fastq
    """
}

process filtlong {
    publishDir '.results/pre-aasembly/filtlong', mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(fastq)

    output:
    path "*.fastq", emit: filtlong_fastq

    script:
    """
    filtlong --min_length 1000 --keep_percent 90 ${fastq} > ${baseName}_filtlong.fastq
    """
}