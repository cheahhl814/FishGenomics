# Nextflow Genomic Pipeline

A Nextflow pipeline for genome assembly (using Oxford Nanopore Technologies (ONT) long reads, implementing multiple assemblers and reconciliation approaches), genome annotation, and comparative genomic analyses.

## Features
- Pre-assembly QC and read processing
- Multiple assembly methods (Canu, wtdbg2, Flye, Raven, Shasta)
- Assembly polishing with Racon
- Assembly evaluation (QUAST, BUSCO)
- Genome reconciliation using RagTag and Quickmerge
- Genome annotation using Funannotate
- Comparative genomic analyses using JCVI
- MultiQC report generation

## Requirements
- Nextflow (DSL2)
- Conda/Mamba/Micromamba (Recommended) (for tool dependencies)
- Input ONT reads (.fastq/.fq.gz)
- Reference genome (fasta, for scaffolding)

## Usage

Basic usage:
```bash
nextflow run main.nf -entry [workflow_name] [parameters]
```
You can also run it directly without cloning/downloading the script using:
```bash
nextflow run cheahhl814/NextGenomics -entry [workflow_name] [parameters]
```

### Parameters (Pre-assembly)
```bash
nextflow run main.nf -entry preAssembly [parameters]
```

### Parameters (Assembly)
```bash
nextflow run main.nf -entry [canu,wtdbg2,flye,raven,shasta] --genomeSize "estimated_size"
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

## Workflows
1. Read QC and filtering: `preAssembly`
2. Genome assembly: `canu`, `wtdbg2`, `flye`, `raven`, `shasta`
3. Genome reconciliation: `reconciliationRagTag` or `reconciliationQuickmerge`

## 
