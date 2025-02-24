process decon {
  tag "Filter out contaminant reads (microbial and human DNA)"
  publishDir "./results/decon", mode: 'copy', overwrite: false, pattern: '**'

  input:
  val(contaminants)
  path(fastq)

  output:
  path "*_decontaminated.fastq.gz", emit: deconFASTQ

  script:
  def sample_id = fastq.baseName
  """
  minimap2 -t ${task.cpus} -ax map-ont ${contaminants} ${fastq} | samtools view -bS - | samtools fastq -f 4 - | gzip > ${sample_id}_decontaminated.fastq.gz)
  """
}