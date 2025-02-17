#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters (Pre-assembly)
params.fastq = "${launchDir}/*.fastq" // Input WGS long read fastq files
params.conFiles = "${launchDir}/contaminants.tsv" // Contaminant reference database
params.resultDir = './results' // Directory for all results

// Parameters (Mitogenome assembly)
params.refmtDNA = "" // Mitochondria reference sequences (FASTA) of closely related species
params.firstGene = "" // FASTA sequences of genes to use as start point
params.orthoDir = "${launchDir}/orthofinder_input"

// Parameters (Assembly)
params.sample_id = ""
params.genomeSize = ""
params.reference_genome = "${launchDir}/*.{fa,fasta,fna}"

// Parameters (Reconciliation)
params.firstA = ''
params.secondA = ''
params.thirdA = ''
params.fourthA = ''
params.fifthA = ''

// Module inclusion
include { getPreassembly; getCanu; getFlye; getRaven; getShasta; getwtdbg2; getRacon; getMummer; getQuast; getBusco; getReconciliation } from './modules/installLocal.nf'
include { identifymtDNA; segregateReads; mtAssembly; mtPolish; mtCircular; mtAnnotate; mtOrtho; trimMSA; mtTree } from './modules/mitochondria.nf'
include { buildIndex; mapReads; filterReads; nanoplot as nanoplot_raw; porechop; filtlong; nanoplot as nanoplot_trimmed } from './modules/pre-assembly.nf'
include { canu; wtdbg2; flye; raven; shasta; racon } from './modules/assembly.nf'
include { scaffold; scaffold2; patch as patch1; patch as patch2; patch as patch3; patch as patch4; quickmerge as quickmerge1; quickmerge as quickmerge2; quickmerge as quickmerge3; quickmerge as quickmerge4 } from './modules/scaffolding.nf'
include { quast; quast as quast_scaffold; busco; busco as busco_scaffold; galignment } from './modules/assessment.nf'
include { multiqc as preassemblyReport; multiqc as mitoassemblyReport; multiqc as assemblyReport } from './modules/multiqc.nf'
include { ganonClassify; decon } from './modules/decon.nf'

// Workflows

workflow installLocal {
  getPreassembly()
  getDecon()
  getCirclator()
  getProkka()
  getOrthofinder()
  getTrimAl()
  getraxmlng
  getCanu()
  getFlye()
  getRaven()
  getShasta()
  getwtdbg2()
  getRacon()
  getMummer()
  getQuast()
  getBusco()
  getReconciliation()
}

workflow preAssembly {
  fastqs = Channel.fromPath("${params.fastq}").collect()
  fastq = Channel.fromPath("${params.fastq}")
  conFiles = Channel.fromPath("${params.conFiles}")

  nanoplot_raw(fastqs)
  porechop(fastq)
  filtlong(porechop.out.porechop_fastq)
  nanoplot_trimmed(filtlong.out.filtlong_fastq.collect())
  ganonClassify(fastq, conFiles)
  decon(ganonClassify.out.one.collect(), fastq)
}

workflow mitoAssembly {
  mitoDNA = Channel.fromPath("${params.refmtDNA}")
  fastq = Channel.fromPath("${params.resultDir}/pre-assembly/decon/*_decontaminated.fastq")
  asmDir = Channel.fromPath("${params.resultDir}/mtGenome/flye", type: 'dir')
  firstGene = Channel.fromPath("${params.firstGene}")
  treeDir = Channel.fromPath("${params.treeDir}", type: 'dir')

  identifymtDNA(fastq, mitoDNA)
  segregateReads(identifymtDNA.out.one.collect(), fastq)
  mtAssembly(segregateReads.out.mitoq.collect(), asmDir)
  mtPolish(mtAssembly.out.mtContig, segregateReads.out.mitoq.collect())
  mtCircular(mtPolish.out.polished_fasta, firstGene)
  mtAnnotate(mtCircular.out.mtFinal)
  mtOrtho(mtAnnotate.out.mtProteinN, orthoDir)
  trimMSA(mtOrtho.out.msa)
  mtTree(trimMSA.out.trimal_fasta)
}

workflow canuWf {
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq")
          .collect()
          .set()
  name = Channel.value(params.sample_id)
  reference_genome = Channel.fromPath(params.reference_genome)
  genomeSize = Channel.value(params.genomeSize)

  canu(fastq, genomeSize, name)
  racon(canu.out.canu_contig, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastq)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow wtdbg2Wf {
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect().set()
  name = Channel.value(params.sample_id)
  reference_genome = Channel.fromPath(params.reference_genome)
  genomeSize = Channel.value(params.genomeSize)

  wtdbg2(fastq, genomeSize, name)
  racon(wtdbg2.out.wtdbg2_contig, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastq)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow flyeWf {
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect().set()
  reference_genome = Channel.fromPath(params.reference_genome)
  genomeSize = Channel.value(params.genomeSize)
  flyeDir = Channel.fromPath("${params.resultDir}/assembly/flye", type: 'dir')

  flye(fastq, genomeSize, flyeDir)
  racon(flye.out.flye_contig, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastq)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow ravenWf {
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect().set()
  reference_genome = Channel.fromPath(params.reference_genome)

  raven(fastq)
  racon(raven.out.raven_contig, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastq)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow shastaWf {
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect().set()
  reference_genome = Channel.fromPath(params.reference_genome)
  name = Channel.value(params.sample_id)

  shasta(fastq, name)
  racon(shasta.out.shasta_contig, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastq)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow reconciliationRagTag {
  // Input channel for genome assembly FASTA files
  scaffold1 = Channel.fromPath(params.firstA)
  scaffold2 = Channel.fromPath(params.secondA)
  scaffold3 = Channel.fromPath(params.thirdA)
  scaffold4 = Channel.fromPath(params.fourthA)
  scaffold5 = Channel.fromPath(params.fifthA)
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect().set()
  reference_genome = Channel.fromPath(params.reference_genome)

  // Genome reconciliation workflow
  patch1(canuScaffold, wtdbg2Scaffold)
  patch2(patch1.out.patch_fasta, flyeScaffold)
  patch3(patch2.out.patch_fasta, ravenScaffold)
  patch4(patch3.out.patch_fasta, shastaScaffold)
  scaffold2(patch4.out.patch_fasta, reference_genome, fastq)
  racon(scaffold.out.scaffold_fasta, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  galignment(racon.out.scaffold_fasta, reference_genome)
}

workflow reconciliationQuickmerge {
  // Input channel for genome assembly FASTA files
  canuScaffold = Channel.fromPath(params.canuScaffold)
  wtdbg2Scaffold = Channel.fromPath(params.wtdbg2Scaffold)
  flyeScaffold = Channel.fromPath(params.flyeScaffold)
  ravenScaffold = Channel.fromPath(params.ravenScaffold)
  shastaScaffold = Channel.fromPath(params.shastaScaffold)
  fastq = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect().set()
  reference_genome = Channel.fromPath(params.reference_genome)

  // Genome reconciliation workflow
  quickmerge1(canuScaffold, wtdbg2Scaffold)
  quickmerge2(quickmerge1.out.merge_assembly, flyeScaffold)
  quickmerge3(quickmerge2.out.merge_assembly, ravenScaffold)
  quickmerge4(quickmerge3.out.merge_assembly, shastaScaffold)
  scaffold2(quickmerge4.out.merge_assembly, reference_genome, fastq)
  racon(scaffold.out.scaffold_fasta, fastq)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  galignment(racon.out.scaffold_fasta, reference_genome)
}

workflow generateReport {
  preassemblyReports = Channel.fromPath("${params.resultDir}/pre-assembly", type: 'dir')
  mitoassemblyReports = Channel.fromPath("${params.resultDir}/mtGenome", type: 'dir')
  assemblyReports = Channel.fromPath("${params.resultDir}/assembly", type: 'dir')

  preassemblyReport(preassemblyReports)
  mitoassemblyReport(mitoassemblyReports)
  assemblyReport(assemblyReports)
}

// Workflow error handling
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
