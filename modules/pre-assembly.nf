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
  micromamba activate preassembly
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

  script:
  def fastq = fastq.collect { "$it" }.join(' ')
  """
  micromamba activate preassembly
  porechop --threads $task.cpus -i $fastq -o ${baseName}_porechop.fastq --format fastq
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
    micromamba activate preassembly
    filtlong --min_length 1000 --keep_percent 90 $fastq | gzip > ${baseName}_filtlong.fastq 2>&1 | tee ${baseName}_filtlong.log
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
    micromamba activate preassembly
    minimap2 -d ${genome_id}.mmi $genome
    """
}

process mapReads {
    tag "Map reads to contaminant genomes (microbial and human DNA)"
    publishDir "./results/pre-assembly/minimap2", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path(index)
    tuple val(baseName), path(fastq)

    output:
    path "*_ids.txt", emit: contaminantsID

    script:
    def contaminant_id = index.baseName
    """
    micromamba activate preassembly
    minimap2 -t ${task.cpus} -ax map-ont $index $fastq > ${contaminant_id}.sam
    samtools view ${contaminant_id}.sam | awk '{print \$1}' | sort | uniq > ${contaminant_id}_ids.txt
    """
}

process filterReads {
    tag "Map reads to contaminant genomes (microbial and human DNA)"
    publishDir "./results/pre-assembly/minimap2", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(baseName), path(contaminantID)
    tuple val(baseName), path(fastq)

    output:
    path "decontaminated.fastq", emit: fastq

    script:
    def ids = contaminantID.collect { "$it" }.join(' ')
    """
    micromamba activate preassembly
    cat $ids | sort | uniq > contaminantID_all.txt
    seqkit grep -j ${task.cpus} -v -f contaminantID_all.txt $fastq > decontaminated.fastq
    """
}