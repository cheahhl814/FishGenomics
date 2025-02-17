process funSetup {
  tag "Set up databases needed by Funannotate"
  publishDir "./results/annotation/funannotate_db", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(fdb)
  
  output:
  path "fun_setup.log"
  
  script:
  """
  funannotate setup -d $fdb -i all -b all --force 2> fun_setup.log
  """
}

process funClean {
  tag "Clean the genome"
  publishDir "./results/annotation", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(genome)

  output:
  path "*.cleaned.fasta", emit: cleanGenome

  script:
  def sample_id = genome.baseName
  """
  funannotate clean -i ${genome} --minlen 1000 -o ${sample_id}.cleaned.fasta
  """
}

process funSort {
  tag "Sort the genome"
  publishDir "./results/annotation", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(cleanGenome)

  output:
  path "*sorted.fasta", emit: sortGenome

  script:
  def sample_id = cleanGenome.baseName
  """
  funannotate sort -i ${cleanGenome} -b scaffold -o ${sample_id}.sorted.fasta
  """
}

process funMask {
  tag "Softmask repeats"
  publishDir "./results/annotation", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(sortGenome)

  output:
  path "*.masked.fasta", emit: maskGenome

  script:
  def sample_id = sortGenome.baseName
  """
  funannotate mask --cpus ${task.cpus} -i ${sortGenome} -o ${sample_id}.masked.fasta
  """
}

process funPredict {
  tag "Predict gene models"
  publishDir "./results/annotation/predict", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(maskGenome)
  val(species)
  val(buscodb)

  output:
  path "predict"

  script:
  """
  funannotate predict --cpus ${task.cpus} -i ${maskGenome} -o predict --species "$species" --buscodb "$buscodb" --organism other --repeats2evm
  """
}

process funAnnotate {
  tag "Annotate the genome"
  publishDir "./results/annotation/annotate", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path "predict"
  val(buscodb)

  output:
  path "annotate"
  path "*.gff", emit: finalGFF

  script:
  """
  funannotate annotate --cpus ${task.cpus} --buscodb "$buscodb" --sbt -i ${predict} -o annotate
  """
}

process annotationStats {
  tag "Generate annotation reports"
  publishDir "./results/annotation/stats", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(finalGFF)

  output:
  path "stats"

  script:
  def sample_id = finalGFF.baseName
  """
  python -m jcvi.annotation.stats genestats ${finalGFF} -o $sample_id
  python -m jcvi.annotation.stats stats ${finalGFF}
  """
}

process repeatScreen {
  tag "Screen for repeats"
  publishDir "./results/annotation/repeats", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(genome)
  path(gff)

  output:
  path "*.freq"
  path "*.fasta.out"
  path "*.fasta.tbl"
  path "*.gff"
  path "*.html"
  path "*.fasta.cat.gz"
  path "*.fasta.masked"
  path "*.align", emit: repeatAlignment
  path "*_final_TE_library.fasta", emit: repeatFASTA

  script:
  def sample_id = genome.baseName
  """
  build_lmer_table -l 14 -sequence $genome -freq ${sample_id}.freq
  RepeatScout -sequence $genome -output ${sample_id}_repeats.fasta -freq ${sample_id}.freq -l 14
  cat ${sample_id}_repeats.fasta | filter-stage-1.prl > ${sample_id}_repeats.filtered1.fasta
  RepeatMasker -pa ${task.cpus} -e ncbi -s -nolow -no_is -gccalc -cutoff 200 -lib ${sample_id}_repeats.filtered1.fasta $genome
  cat ${sample_id}_repeats.filtered1.fasta | filter-stage-2.prl --cat=${sample_id}.fasta.out --thresh=10 > ${sample_id}_repeats.filtered2.fasta
  compare-out-to-gff.prl --gff=$gff --cat=${sample_id}.fasta.out --f=$genome > ${sample_id}_final_TE_library.fasta
  RepeatMasker -pa ${task.cpus} -e ncbi -s -nolow -no_is -gccalc -cutoff 200 -gff -html -a -lib ${sample_id}_final_TE_library.fasta $genome
  """
}

process repeatLandscape {
  tag "Evaluate the repeat landscape"
  publishDir "./results/annotation/repeats/landscape", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(repeatAlignment)

  output:
  path "*.divsum"
  path "*_repeatLandscape.html"

  script:
  def sample_id = repeatAlignment.baseName
  """
  calcDivergenceFromAlign.pl -s ${sample_id}.divsum $repeatAlignment
  createRepeatLandscape.pl -div ${sample_id}.divsum > ${sample_id}_repeatLandscape.html
  """
}

process classifyRepeats {
  tag "Annotate repeats"
  publishDir "./results/annotation/repeats/annotation", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(repeatFASTA)

  output:
  path "*.classified"

  script:
  def sample_id = repeatFASTA.baseName
  """
  RepeatClassifier -consensi $repeatFASTA -engine ncbi
  """
}