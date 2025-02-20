process decon {
  tag "Filter out contaminant reads (microbial and human DNA)"
  publishDir "./results/decon", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contaminants)
  tuple val(baseName), path(fastq)

  output:
  path "*_decontaminated.fastq", emit: deconFASTQ

  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont ${contaminants} ${fastq} | samtools view -b -f 4 - | samtools fastq - > ${baseName}_decontaminated.fastq
  """
}