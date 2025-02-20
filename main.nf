#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters (Pre-assembly)
params.fastq = "${launchDir}/reads/*.fastq.gz" // Input WGS long read fastq files
params.conFasta = "${launchDir}/contaminants.fasta" // Contaminant reference database
params.resultDir = "./results" // Directory for all results

// Parameters (Mitogenome assembly)
params.refmtDNA = "${launchDir}/referenceMt/*.{fa,fasta,fna}" // Mitochondria reference sequences (FASTA) of closely related species
params.firstGene = "${launchDir}/referenceMt/firstGene.{fa,fasta,fna}" // FASTA sequences of genes to use as start point
params.orthoMt = "${launchDir}/orthofinderMt" // Input folder for Orthofinder

// Parameters (Assembly)
params.sample_id = ""
params.genomeSize = ""
params.reference_genome = "${launchDir}/referenceGenome/*.{fa,fasta,fna}" // Reference genome of a closely related species

// Parameters (Reconciliation)
params.firstA = ""
params.secondA = ""
params.thirdA = ""
params.fourthA = ""
params.fifthA = ""

// Parameters (Annotation)
params.finalAsm = "" // Final assembly after reconciliation
params.species = ""
params.buscodb = ""

// Module inclusion
include { getPreassembly; getDecon; getCirclator; getProkka; getOrthofinder; getTrimAl; getraxmlng; getCanu; getFlye; getRaven; getShasta; getwtdbg2; getRacon; getMummer; getQuast; getBusco; getReconciliation; getFunannotate } from './modules/installLocal.nf'
include { identifymtDNA; segregateReads; mtAssembly; mtPolish; mtCircular; mtAnnotate; mtOrtho; trimMSA; mtTree } from './modules/mitochondria.nf'
include { buildIndex; mapReads; filterReads; nanoplot as nanoplot_raw; porechop; filtlong; nanoplot as nanoplot_trimmed } from './modules/pre-assembly.nf'
include { canu; wtdbg2; flye; raven; shasta; racon } from './modules/assembly.nf'
include { scaffold; scaffold2; patch as patch1; patch as patch2; patch as patch3; patch as patch4; quickmerge as quickmerge1; quickmerge as quickmerge2; quickmerge as quickmerge3; quickmerge as quickmerge4 } from './modules/scaffolding.nf'
include { quast; quast as quast_scaffold; busco; busco as busco_scaffold; galignment } from './modules/assessment.nf'
include { multiqc as preassemblyReport; multiqc as mitoassemblyReport; multiqc as assemblyReport } from './modules/multiqc.nf'
include { decon } from './modules/decon.nf'
include { funClean; funSort; funMask; funPredict; funAnnotate; annotationStats } from './modules/annotate.nf'

// Workflows

workflow installLocal {
  getPreassembly()
  getDecon()
  getCirclator()
  getProkka()
  getOrthofinder()
  getTrimAl()
  getraxmlng()
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
  getFunannotate()
}

workflow deconOnly {
  fastq = Channel.fromPath("${params.fastq}")
          .map { file -> tuple(file.simpleName, file) }
          .view()
  conFasta = Channel.fromPath("${params.conFasta}").view()

  nanoplot_raw(fastq.collect())
  decon(conFasta, fastq)
}

workflow preAssembly {
  fastq_files = Channel.fromPath("${params.fastq}").collect()
  fastqs = Channel.value(fastq_files)
  fastq = Channel.fromPath("${params.fastq}")
  conFasta = Channel.fromPath("${params.conFasta}")

  nanoplot_raw(fastqs)
  porechop(fastq)
  filtlong(porechop.out.porechop_fastq)
  nanoplot_trimmed(filtlong.out.filtlong_fastq.collect())
  decon(conFasta, filtlong.out.filtlong_fastq)
}

workflow mitoAssembly {
  mitoDNA = Channel.fromPath("${params.refmtDNA}")
  fastq = Channel.fromPath("${params.resultDir}/pre-assembly/decon/*_decontaminated.fastq")
  firstGene = Channel.fromPath("${params.firstGene}")
  orthoMt = Channel.fromPath("${params.orthoMt}", type: 'dir')

  identifymtDNA(fastq, mitoDNA)
  segregateReads(identifymtDNA.out.one.collect(), fastq)
  mtAssembly(segregateReads.out.mitoq.collect(), asmDir)
  mtPolish(mtAssembly.out.mtContig, segregateReads.out.mitoq.collect())
  mtCircular(mtPolish.out.polished_fasta, firstGene)
  mtAnnotate(mtCircular.out.mtFinal)
  mtOrtho(mtAnnotate.out.mtProteinN, orthoMt)
  trimMSA(mtOrtho.out.msa)
  mtTree(trimMSA.out.trimal_fasta)
}

workflow canuWf {
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  name = Channel.value("${params.sample_id}")
  reference_genome = Channel.fromPath("${params.reference_genome}")
  genomeSize = Channel.value("${params.genomeSize}")

  canu(fastqs, genomeSize, name)
  racon(canu.out.canu_contig, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastqs)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow wtdbg2Wf {
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  name = Channel.value("${params.sample_id}")
  reference_genome = Channel.fromPath("${params.reference_genome}")
  genomeSize = Channel.value("${params.genomeSize}")

  wtdbg2(fastqs, genomeSize, name)
  racon(wtdbg2.out.wtdbg2_contig, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastqs)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow flyeWf {
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  reference_genome = Channel.fromPath("${params.reference_genome}")
  genomeSize = Channel.value("${params.genomeSize}")
  flyeDir = Channel.fromPath("${params.resultDir}/assembly/flye", type: 'dir')

  flye(fastq, genomeSize, flyeDir)
  racon(flye.out.flye_contig, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastqs)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow ravenWf {
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  reference_genome = Channel.fromPath("${params.reference_genome}")

  raven(fastq)
  racon(raven.out.raven_contig, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastqs)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow shastaWf {
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  reference_genome = Channel.fromPath("${params.reference_genome}")
  name = Channel.value("${params.sample_id}")

  shasta(fastqs, name)
  racon(shasta.out.shasta_contig, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, fastqs)
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow reconciliationRagTag {
  // Input channel for genome assembly FASTA files
  scaffold1 = Channel.fromPath("${params.firstA}")
  scaffold2 = Channel.fromPath("${params.secondA}")
  scaffold3 = Channel.fromPath("${params.thirdA}")
  scaffold4 = Channel.fromPath("${params.fourthA}")
  scaffold5 = Channel.fromPath("${params.fifthA}")
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  reference_genome = Channel.fromPath("${params.reference_genome}")

  // Genome reconciliation workflow
  patch1(canuScaffold, wtdbg2Scaffold)
  patch2(patch1.out.patch_fasta, flyeScaffold)
  patch3(patch2.out.patch_fasta, ravenScaffold)
  patch4(patch3.out.patch_fasta, shastaScaffold)
  scaffold2(patch4.out.patch_fasta, reference_genome, fastqs)
  racon(scaffold.out.scaffold_fasta, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  galignment(racon.out.scaffold_fasta, reference_genome)
}

workflow reconciliationQuickmerge {
  // Input channel for genome assembly FASTA files
  scaffold1 = Channel.fromPath("${params.firstA}")
  scaffold2 = Channel.fromPath("${params.secondA}")
  scaffold3 = Channel.fromPath("${params.thirdA}")
  scaffold4 = Channel.fromPath("${params.fourthA}")
  scaffold5 = Channel.fromPath("${params.fifthA}")
  fastq_files = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq").collect()
  fastqs = Channel.value(fastq_files)
  reference_genome = Channel.fromPath("${params.reference_genome}")

  // Genome reconciliation workflow
  quickmerge1(canuScaffold, wtdbg2Scaffold)
  quickmerge2(quickmerge1.out.merge_assembly, flyeScaffold)
  quickmerge3(quickmerge2.out.merge_assembly, ravenScaffold)
  quickmerge4(quickmerge3.out.merge_assembly, shastaScaffold)
  scaffold2(quickmerge4.out.merge_assembly, reference_genome, fastqs)
  racon(scaffold.out.scaffold_fasta, fastqs)
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  galignment(racon.out.scaffold_fasta, reference_genome)
}

workflow annotation {
  // Input channel for the final genome assembly
  input_assembly = Channel.fromPath("${params.finalAsm}")
  species = Channel.value("${params.species}")
  buscodb = Channel.value("${params.buscodb}")

  // Annotation workflow
  funClean(input_assembly)
  funSort(funClean.out.cleanGenome)
  funMask(funSort.out.sortGenome)
  funPredict(funMask.out.maskGenome, species, buscodb)
  funAnnotate(funPredict.out, buscodb)
  annotationStats(annotateGenes.out.finalGFF)
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
