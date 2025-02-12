process scaffold {
  tag "Scaffolding with RagTag"
  publishDir "./results/assembly/scaffold", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contig)
  path(reference)
  path(fastq)

  output:
  path "ragtag.correct.fasta", emit: correct_fasta
  path "ragtag.correct.agp", emit: correct_agp
  path "ragtag.scaffold.agp", emit: scaffold_agp
  path "ragtag.scaffold.fasta", emit: scaffold_fasta
  path "ragtag.scaffold.stats", emit: scaffold_stats

  script:
  def sample_id = contig.baseName
  """
  ragtag.py correct $reference $contig -t ${task.cpus} -R $fastq -T ont
  ragtag.py scaffold $reference $contig -r -t ${task.cpus}
  """
}

process patch {
  tag "Patch genome assemblies"
  publishDir "./results/assembly/reconciliation", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(target)
  path(query)

  output:
  path "ragtag.patch.agp", emit: patch_agp
  path "ragtag.patch.fasta", emit: patch_fasta

  script:
  """
  ragtag.py patch $target $query -t ${task.cpus}
  """
}

process scaffold2 {
  tag "Second scaffolding after reconciliation with RagTag"
  publishDir "./results/assembly/reconciliation/ragtag", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(contig)
  path(reference)
  path(fastq)

  output:
  path "ragtag.correct.fasta", emit: correct_fasta
  path "ragtag.correct.agp", emit: correct_agp
  path "ragtag.scaffold.agp", emit: scaffold_agp
  path "ragtag.scaffold.fasta", emit: scaffold_fasta
  path "ragtag.scaffold.stats", emit: scaffold_stats

  script:
  def sample_id = contig.baseName
  """
  ragtag.py correct $reference $contig -t ${task.cpus} -R $fastq -T ont
  ragtag.py scaffold $reference $contig -r -t ${task.cpus}
  """
}

process quickmerge {
  tag "Merge genome assemblies"
  publishDir "./results/assembly/reconciliation/quickmerge", mode: 'copy', overwrite: false, pattern: '**'

  input:
  path(primary)
  path(secondary)

  output:
  path "*.fasta", emit: merge_assembly

  script:
  def primary_id = primary.baseName
  def secondary_id = secondary.baseName
  """
  merge_wrapper.py $primary $secondary --prefix ${primary_id}_${secondary_id}
  """
}
