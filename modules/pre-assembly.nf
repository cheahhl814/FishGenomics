process nanoplot {
  tag "Create relevant plots from ONT reads using Nanoplot"
  publishDir "./results/pre-assembly/nanoplot", mode: 'copy', overwrite: false, pattern: '**'

  input:
  value(fastq)

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

process buildIndex {
    tag "Mapping long reads"
    publishDir "./results/pre-assembly/minimap2", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(genome)

    output:
    path "*.mmi", emit: mmi

    script:
    def genome_id = genome.baseName
    """
    minimap2 -d ${genome_id}.mmi $genome
    """
}

process mapReads {
    tag "Map reads to contaminant genomes (microbial and human DNA)"
    publishDir "./results/pre-assembly/minimap2", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(index)
    path(fastq)

    output:
    path "*_ids.txt", emit: contaminantsID

    script:
    def contaminant_id = index.baseName
    def sample_id = fastq.baseName
    """
    minimap2 -t ${task.cpus} -ax map-ont ${index} ${fastq} > ${sample_id}_${contaminant_id}.sam
    samtools view -F 4 ${sample_id}_${contaminant_id}.sam | awk '{print \$1}' | sort | uniq > ${sample_id}_${contaminant_id}_ids.txt
    """
}

process filterReads {
    tag "Filter out contaminant reads"
    publishDir "./results/pre-assembly/minimap2", mode: 'copy', overwrite: false, pattern: '**'

    input:
    val(contaminantID)
    path(fastq)

    output:
    path "decontaminated.fastq", emit: fastq

    script:
    def ids = contaminantID.join(" ")
    def sample_id = fastq.baseName
    """
    cat ${ids} | sort | uniq | seqkit grep -j ${task.cpus} -v -f - ${fastq} > ${sample_id}_decontaminated.fastq
    """
}