# Nextflow Genomic Pipeline

A Nextflow pipeline for whole genome assembly (using Oxford Nanopore Technologies (ONT) long reads, implementing multiple assemblers and reconciliation approaches).

## Features
- Pre-assembly QC and read processing
- Removal of contaminant reads
- Mitochondrial genome assembly, annotation, and phylogenetic analysis
- Multiple WGA methods (Canu, wtdbg2, Flye, Raven, Shasta)
- Polishing with Racon
- Genome reconciliation using RagTag or Quickmerge
- Assembly evaluation (QUAST, BUSCO)
- Genome annotation (Funannotate)
- MultiQC report generation

## General Requirements
- Nextflow (DSL2)
    ```bash
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow $HOME/.local/bin/
    ```
- Micromamba (Recommended) (for tool dependencies). To install micromamba, run:
    ```bash
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
    ```
- Tool dependencies. To install dependencies, run:
    ```bash
    nextflow run cheahhl814/FishGenomics -entry installLocal -profile local
    ```

## Workflow-specific Requirements
### deconOnly and preAssembly Workflow
- Input ONT reads (`.fastq`, `fq`, `fasta`, or `fa`) in current directory.
- Contaminant genomes (`./contaminants/*.fasta`).

### mitoAssembly
- Reference (complete) mitochondrial genome (`./referenceMt/*.{fa,fasta,fna}`)
- DNA sequence of the first gene in mitochondrial genomes (`/referenceMt/firstGene.{fa,fasta,fna}`)
- Input folder for Orthofinder (mitochondrial genome) (`./orthofinderMt`). This folder should contain nucleotide FASTA files of core genes (13 protein-coding genes + 2 ribosomal RNA genes) from the target species, closely related species (same genus), and at least one outgroup species.

### canuWf, wtdbg2Wf, flyeWf, ravenWf, and shastaWf Workflows
- Reference genome (`./referenceGenome/*.fasta`, for scaffolding).
- `--sample_id` "prefix" (Prefix of genome assembly)
- `--genomeSize` "size" (Expected genome size)

### reconciliationRagTag and reconciliationQuickmerge Workflows
- `--firstA` "First assembly FASTA"
- `--secondA` = "Second assembly FASTA"
- `--thirdA` "Third assembly FASTA"
- `--fourthA` "Fourth assembly FASTA"
- `--fifthA` "Fifth assembly FASTA"

## Usage

Basic usage:
```bash
nextflow run main.nf -entry [workflow_name] [parameters]
```
You can also run it directly without cloning/downloading the script using:
```bash
nextflow run cheahhl814/FishGenomics -entry [workflow_name] [parameters]
```
You can select whether to run it on local desktop (using Micromamba) or HPC system (SLURM):
```bash
nextflow run cheahhl814/FishGenomics -entry [workflow_name] -profile [local,hpc] [parameters]
```

### Parameters (Pre-assembly)
```bash
nextflow run main.nf -entry preAssembly [parameters]
```
### Parameters (Mitochondrial genome assembly)
```bash
nextflow run main.nf -entry mitoAssembly [parameters]
```
### Parameters (Assembly)
```bash
nextflow run main.nf -entry [canuWf,wtdbg2Wf,flyeWf,ravenWf,shastaWf] --genomeSize "estimated_size"
```
- `--genomeSize`: Estimated genome size (required for some assemblers)

### Parameters (Genome reconciliation)
```bash
nextflow run main.nf -entry [reconciliationRagTag,reconciliationQuickmerge]
```
- `--canuScaffold`: Path to scaffold fasta file from the Canu pipeline
- `--wtdbg2Scaffold`: Path to scaffold fasta file from the wtdbg2 pipeline
- `--flyeScaffold`: Path to scaffold fasta file from the Flye pipeline
- `--ravenScaffold`: Path to scaffold fasta file from the Raven pipeline
- `--shastaScaffold`: Path to scaffold fasta file from the Shasta pipeline

### Output Directories
- `./results/`

### Custom Configuration
We understand that some tools required in this pipeline may be not available on your HPC system. In that case, you can run the pipeline with your custom configuration file (template: https://github.com/cheahhl814/FishGenomics/blob/main/customNextflow.config):
```bash
nextflow run cheahhl814/FishGenomics -entry [workflow_name] [parameters] -c [customConfigFile]
```

## Workflows
1. Install dependencies: `installLocal`
2. Read QC plot and decontamination only: `deconOnly`
3. Full pre-assembly workflow (Read QC, filtering, and decontamination): `preAssembly`
4. Mitochondrial genome assembly, annotation, and phylogenetic analysis: `mitoAssembly`
5. Genome assembly: `canuWf`, `wtdbg2Wf`, `flyeWf`, `ravenWf`, `shastaWf`
6. Genome reconciliation: `reconciliationRagTag` or `reconciliationQuickmerge`
7. Genome annotation: `annotation`
