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
