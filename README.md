# Nextflow Genomic Pipeline

A Nextflow pipeline for whole genome assembly (using Oxford Nanopore Technologies (ONT) long reads, implementing multiple assemblers and reconciliation approaches).

## Features
- Pre-assembly QC and read processing
- Multiple assembly methods (Canu, wtdbg2, Flye, Raven, Shasta)
- Assembly polishing with Racon
- Assembly evaluation (QUAST, BUSCO)
- Genome reconciliation using RagTag and Quickmerge
- MultiQC report generation

## Requirements
- Nextflow (DSL2)
    ```bash
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow $HOME/.local/bin/
    ```
- Micromamba (Recommended) (for tool dependencies)
    ```bash
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
    ```
- Tool dependencies
    ```bash
    nextflow run cheahhl814/FishGenomics -entry install -profile local
    ```
Input ONT reads (.fastq) in current directory
- Contaminant genomes (.fasta) in the `contaminants/` folder of current directory
- Reference genome (.fasta, for scaffolding) in current directory

## Usage

Basic usage:
```bash
nextflow run main.nf -entry [workflow_name] [parameters]
```
You can also run it directly without cloning/downloading the script using:
```bash
nextflow run cheahhl814/FishGenomics -entry [workflow_name] [parameters]
```

### Parameters (Pre-assembly)
```bash
nextflow run main.nf -entry preAssembly [parameters]
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

## Workflows
1. Read QC and filtering: `preAssembly`
2. Genome assembly: `canuWf`, `wtdbg2Wf`, `flyeWf`, `ravenWf`, `shastaWf`
3. Genome reconciliation: `reconciliationRagTag` or `reconciliationQuickmerge`

## 
