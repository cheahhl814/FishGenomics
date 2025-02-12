process nanoplot {
  tag "Create relevant plots from ONT reads using Nanoplot"
  publishDir "./results/pre-assembly/nanoplot", mode: 'copy', overwrite: false, pattern: '**'

  input:
  tuple val(baseName), path(fastq)

  output:
  path "*.tsv", emit: tsv
  path "*.png", emit: png
  path "*.html", emit: html

  script:
  def fastq = fastq.collect { "$it" }.join(' ')
  """
  NanoPlot --threads ${task.cpus} --fastq $fastq --maxlength 40000 --tsv_stats --plots dot --format png --legacy hex --info_in_report --prefix $baseName
  """
}

process porechop {
  tag "Trim adapters"
  publishDir "./results/pre-assembly/porechop", mode: 'copy', overwrite: false, pattern: '**'

  input:
  tuple val(baseName), path(fastq)

  output:
  path "*_porechop.fastq", emit: porechop_fastq
  path "*_porechop.log", emit: porechop_log

  script:
  def fastq = fastq.collect { "$it" }.join(' ')
  """
  porechop --threads $task.cpus -i $fastq -o output_${baseName}_porechop.fastq --format fastq 2>&1 | tee ${baseName}_porechop.log
  """
}

process filtlong {
    tag "Quality filtering of reads"
    publishDir "./results/pre-aasembly/filtlong", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(baseName), path(fastq)

    output:
    path "*.fastq", emit: filtlong_fastq
    path "*.log", emit: filtlong_log

    script:
    def fastq = fastq.collect { "$it" }.join(' ')
    """
    filtlong --min_length 1000 --keep_percent 90 $fastq | gzip > ${baseName}_filtlong.fastq 2>&1 | tee ${baseName}_filtlong.log
    """
}

process removeContaminant {
  tag "Filter out contaminant reads (microbial and human DNA)"
  publishDir "./results/pre-assembly/minimap2", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(minimapIndex)
  tuple val(baseName), path(fastq)

  output:
  path "*_decontaminated.fastq", emit: deconFASTQ

  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont contaminantIndex.mmi $fastq > ${baseName}_contaminant.sam
  samtools view -b -f 4 ${baseName}_contaminant.sam > ${baseName}_decontaminated.bam
  samtools ${baseName}_decontaminated.bam > ${baseName}_decontaminated.fastq
  """
}