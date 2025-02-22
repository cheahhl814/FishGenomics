#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters (Pre-assembly)
params.fastq = "${launchDir}/reads/*.fastq.gz" // Input WGS long read fastq files
params.conFasta = "${launchDir}/contaminants.fasta" // Contaminant reference database
params.resultDir = "./results" // Directory for all results

// Parameters (Mitogenome assembly)
params.refmtDNA = "${launchDir}/referenceMt" // Mitochondria reference sequences (FASTA), annotation, and MITOS reference data of closely related species
params.orthoMt = "${launchDir}/orthofinderMt" // Input folder for Orthofinder

// Parameters (Assembly)
params.sample_id = "" // Some tools require prefix
params.genomeSize = "" // Most assembly tools require genome size of the target species
params.reference_genome = "${launchDir}/referenceGenome/*.{fa,fasta,fna}" // Reference genome of a closely related species

// Parameters (Reconciliation)
params.firstA = "" // The best assembly on all parameters
params.secondA = "" // Second best
params.thirdA = "" // Third
params.fourthA = "" // Fourth
params.fifthA = "" // Last

// Parameters (Annotation)
params.finalAsm = "" // Final assembly after reconciliation
params.species = ""
params.buscodb = ""

// Module inclusion
include { pipTools } from './modules/installLocal.nf'
include { nanoplot } from './modules/pre-assembly.nf'
include { segregate; mtAssembly; mtPolish; mtCircular; mtAnnotate; orthoSetup; mtOrtho; orthoFinder; trimMSA; mtTree } from './modules/mitochondria.nf'
include { canu; wtdbg2; flye; raven; shasta; racon } from './modules/assembly.nf'
include { scaffold; scaffold2; patch as patch1; patch as patch2; patch as patch3; patch as patch4; quickmerge as quickmerge1; quickmerge as quickmerge2; quickmerge as quickmerge3; quickmerge as quickmerge4 } from './modules/scaffolding.nf'
include { quast; quast as quast_scaffold; busco; busco as busco_scaffold; galignment } from './modules/assessment.nf'
include { multiqc as preassemblyReport; multiqc as mitoassemblyReport; multiqc as assemblyReport } from './modules/multiqc.nf'
include { decon } from './modules/decon.nf'
include { funClean; funSort; funMask; funPredict; funAnnotate; annotationStats } from './modules/annotate.nf'

// Workflows

workflow extraTools {
  pipTools()
}

workflow plotNdecon {
  reads = Channel.fromPath("${params.fastq}")
  conFasta = Channel.value("${params.conFasta}")

  nanoplot(reads.collect())
  decon(conFasta, reads)
}

workflow mitoAssembly {
  mitoDNA = Channel.value("${params.refmtDNA}/*_mt.{fa,fasta,fna}") // Mitochondrial DNAs of related species
  reads = Channel.fromPath("${params.resultDir}/decon/*_decontaminated.fastq.gz") // Decontaminated reads
  firstGene = Channel.value("${params.refmtDNA}/firstGene.{fna,fa,fasta}") // Desired first gene in mitochondrial DNA circularization, e.g., COI.
  refseq = Channel.value("refseq63f") // Reference folder name for MITOS2
  refseqDir = Channel.value("${params.refmtDNA}", type: 'dir') // Parent directory for refseq63f
  orthoFasta = Channel.fromPath("${params.orthoMt}/*.fasta") // FASTA files of mitochondrial DNAs of related species
  orthoGFF = Channel.fromPath("${params.orthoMt}/*.gff3") // GFF files of mitochondrial genes of related species
  match = orthoFasta.map { file -> tuple(file.baseName, file) }
    .combine(orthoGFF.map { file -> tuple(file.baseName, file) })
    .filter { it[0] == it[2] }
    .map { it -> tuple(it[0], it[1], it[3]) }
  orthofinderInput = Channel.fromPath("${params.resultDir}/mtGenome/phylogenetics/input", type: 'dir') // Input folder for OrthoFinder

  segregate(reads, mitoDNA)
  mtAssembly(segregate.out.mitoq.collect())
  mtPolish(mtAssembly.out.mtContig, segregate.out.mitoq.collect())
  mtCircular(mtPolish.out.polished_fasta, firstGene)
  mtAnnotate(mtCircular.out.mtFinal, refseq, refseqDir)
  orthoSetup(match)
  mtOrtho(mtAnnotate.out.mtGenes, mtCircular.out.mtFinal)
  orthoFinder(orthofinderInput, orthoSetup.out.geneFasta, mtOrtho.out.geneFasta)
  trimMSA(orthoFinder.out.msa)
  mtTree(trimMSA.out.trimal_fasta)
}

workflow canuWf {
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  name = Channel.value("${params.sample_id}")
  reference_genome = Channel.value("${params.reference_genome}")
  genomeSize = Channel.value("${params.genomeSize}")

  canu(reads.collect(), genomeSize, name)
  racon(canu.out.canu_contig, reads.collect())
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, reads.collect())
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow wtdbg2Wf {
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  name = Channel.value("${params.sample_id}")
  reference_genome = Channel.value("${params.reference_genome}")
  genomeSize = Channel.value("${params.genomeSize}")

  wtdbg2(reads.collect(), genomeSize, name)
  racon(wtdbg2.out.wtdbg2_contig, reads.collect())
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, reads.collect())
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow flyeWf {
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  reference_genome = Channel.value("${params.reference_genome}")
  genomeSize = Channel.value("${params.genomeSize}")
  flyeDir = Channel.fromPath("${params.resultDir}/assembly/flye", type: 'dir')

  flye(reads.collect(), genomeSize, flyeDir)
  racon(flye.out.flye_contig, reads.collect())
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, reads.collect())
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow ravenWf {
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  reference_genome = Channel.fromPath("${params.reference_genome}")

  raven(reads.collect())
  racon(raven.out.raven_contig, reads.collect())
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, reads.collect())
  quast_scaffold(scaffold.out.scaffold_fasta, reference_genome)
  busco_scaffold(scaffold.out.scaffold_fasta)
  galignment(scaffold.out.scaffold_fasta, reference_genome)
}

workflow shastaWf {
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  reference_genome = Channel.fromPath("${params.reference_genome}")
  name = Channel.value("${params.sample_id}")

  shasta(reads.collect(), name)
  racon(shasta.out.shasta_contig, reads.collect())
  quast(racon.out.polished_fasta, reference_genome)
  busco(racon.out.polished_fasta)
  scaffold(racon.out.polished_fasta, reference_genome, reads.collect())
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
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  reference_genome = Channel.fromPath("${params.reference_genome}")

  // Genome reconciliation workflow
  patch1(scaffold1, scaffold2)
  patch2(patch1.out.patch_fasta, scaffold3)
  patch3(patch2.out.patch_fasta, scaffold4)
  patch4(patch3.out.patch_fasta, scaffold5)
  scaffold2(patch4.out.patch_fasta, reference_genome, reads.collect())
  racon(scaffold2.out.scaffold_fasta, reads.collect())
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
  reads = Channel.fromPath("${params.resultDir}/mtGenome/*_nuclear.fastq.gz")
  reference_genome = Channel.value("${params.reference_genome}")

  // Genome reconciliation workflow
  quickmerge1(canuScaffold, wtdbg2Scaffold)
  quickmerge2(quickmerge1.out.merge_assembly, flyeScaffold)
  quickmerge3(quickmerge2.out.merge_assembly, ravenScaffold)
  quickmerge4(quickmerge3.out.merge_assembly, shastaScaffold)
  scaffold2(quickmerge4.out.merge_assembly, reference_genome, reads.collect())
  racon(scaffold2.out.scaffold_fasta, reads.collect())
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
