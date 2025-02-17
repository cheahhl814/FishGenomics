#!/bin/bash
set -euo pipefail

# Set launch directory (if not defined elsewhere)
LAUNCH_DIR=$(pwd)

# Print general usage information.
usage() {
    cat <<EOF
Usage: $0 <workflow> [options]

Available workflows:
  installLocal            Run the installLocal workflow
  preAssembly             Run the preAssembly workflow
  mitoAssembly            Run the mitoAssembly workflow
  canuWf                  Run the canuWf workflow
  wtdbg2Wf                Run the wtdbg2Wf workflow
  flyeWf                  Run the flyeWf workflow
  ravenWf                 Run the ravenWf workflow
  shastaWf                Run the shastaWf workflow
  reconciliationRagTag    Run the reconciliationRagTag workflow
  reconciliationQuickmerge Run the reconciliationQuickmerge workflow
  generateReport          Run the generateReport workflow

General Options:
  --help                  Show this help message.

For workflow-specific help, run:
  $0 <workflow> --help

EOF
    exit 1
}

# Ensure at least one argument is provided.
if [ "$#" -lt 1 ]; then
    usage
fi

# First argument is the workflow name.
WORKFLOW="$1"
shift

# If the next argument is '--help', show workflow-specific usage.
if [ "${1:-}" == "--help" ]; then
    case "$WORKFLOW" in
        installLocal)
            cat <<EOF
Usage: $0 installLocal [nextflow options]
This workflow installs all required modules and tools locally.
(No additional parameters required.)
EOF
            exit 0
            ;;
        preAssembly)
            cat <<EOF
Usage: $0 preAssembly [--fastq <pattern>] [--conFiles <path>] [--resultDir <dir>]
This workflow performs pre-assembly operations.
Parameters:
  --fastq     Input FASTQ file pattern (default: "${LAUNCH_DIR}/*.fastq")
  --conFiles  Path to contaminants file (default: "${LAUNCH_DIR}/contaminants.tsv")
  --resultDir Directory to store results (default: "./results")
EOF
            exit 0
            ;;
        mitoAssembly)
            cat <<EOF
Usage: $0 mitoAssembly [--refmtDNA <path>] [--firstGene <path>] [--orthoDir <dir>]
                 [--resultDir <dir>] [--treeDir <dir>]
This workflow performs mitogenome assembly.
Parameters:
  --refmtDNA   Reference mitochondria FASTA (default: "")
  --firstGene  FASTA for gene start (default: "")
  --orthoDir   Orthofinder input directory (default: "${LAUNCH_DIR}/orthofinder_input")
  --resultDir  Results directory (default: "./results")
  --treeDir    Directory for tree outputs (if used)
EOF
            exit 0
            ;;
        canuWf)
            cat <<EOF
Usage: $0 canuWf [--sample_id <id>] [--genomeSize <size>] [--reference_genome <pattern>]
This workflow performs nuclear genome assembly using canu.
Parameters:
  --sample_id         Sample identifier (default: "")
  --genomeSize        Genome size (default: "")
  --reference_genome  Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
EOF
            exit 0
            ;;
        wtdbg2Wf)
            cat <<EOF
Usage: $0 wtdbg2Wf [--sample_id <id>] [--genomeSize <size>] [--reference_genome <pattern>]
This workflow performs nuclear genome assembly using wtdbg2.
Parameters:
  --sample_id         Sample identifier
  --genomeSize        Genome size
  --reference_genome  Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
EOF
            exit 0
            ;;
        flyeWf)
            cat <<EOF
Usage: $0 flyeWf [--genomeSize <size>] [--reference_genome <pattern>] [--resultDir <dir>]
This workflow performs nuclear genome assembly using flye.
Parameters:
  --genomeSize        Genome size
  --reference_genome  Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
  --resultDir         Results directory (default: "./results")
EOF
            exit 0
            ;;
        ravenWf)
            cat <<EOF
Usage: $0 ravenWf [--reference_genome <pattern>]
This workflow performs nuclear genome assembly using raven.
Parameters:
  --reference_genome  Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
EOF
            exit 0
            ;;
        shastaWf)
            cat <<EOF
Usage: $0 shastaWf [--sample_id <id>] [--reference_genome <pattern>]
This workflow performs nuclear genome assembly using shasta.
Parameters:
  --sample_id         Sample identifier
  --reference_genome  Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
EOF
            exit 0
            ;;
        reconciliationRagTag)
            cat <<EOF
Usage: $0 reconciliationRagTag [--firstA <path>] [--secondA <path>] [--thirdA <path>]
                           [--fourthA <path>] [--fifthA <path>] [--reference_genome <pattern>]
This workflow performs genome reconciliation using RagTag.
Parameters:
  --firstA          First assembly scaffold file
  --secondA         Second assembly scaffold file
  --thirdA          Third assembly scaffold file
  --fourthA         Fourth assembly scaffold file
  --fifthA          Fifth assembly scaffold file
  --reference_genome Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
EOF
            exit 0
            ;;
        reconciliationQuickmerge)
            cat <<EOF
Usage: $0 reconciliationQuickmerge [--canuScaffold <path>] [--wtdbg2Scaffold <path>]
                                  [--flyeScaffold <path>] [--ravenScaffold <path>]
                                  [--shastaScaffold <path>] [--reference_genome <pattern>]
This workflow performs genome reconciliation using Quickmerge.
Parameters:
  --canuScaffold     Assembly scaffold from canu
  --wtdbg2Scaffold   Assembly scaffold from wtdbg2
  --flyeScaffold     Assembly scaffold from flye
  --ravenScaffold    Assembly scaffold from raven
  --shastaScaffold   Assembly scaffold from shasta
  --reference_genome Reference genome file pattern (default: "${LAUNCH_DIR}/*.{fa,fasta,fna}")
EOF
            exit 0
            ;;
        generateReport)
            cat <<EOF
Usage: $0 generateReport [--resultDir <dir>]
This workflow generates various reports.
Parameters:
  --resultDir  Base results directory (default: "./results")
EOF
            exit 0
            ;;
        *)
            echo "Unknown workflow: $WORKFLOW"
            usage
            ;;
    esac
fi

# Set default values (using LAUNCH_DIR for defaults)
FASTQ_PATTERN="${LAUNCH_DIR}/*.fastq"
CON_FILES="${LAUNCH_DIR}/contaminants.tsv"
RESULT_DIR="./results"
REF_MTDNA=""
FIRST_GENE=""
ORTHO_DIR="${LAUNCH_DIR}/orthofinder_input"
SAMPLE_ID=""
GENOME_SIZE=""
REFERENCE_GENOME="${LAUNCH_DIR}/*.{fa,fasta,fna}"
FIRST_A=""
SECOND_A=""
THIRD_A=""
FOURTH_A=""
FIFTH_A=""
CANU_SCAFFOLD=""
WTDBG2_SCAFFOLD=""
FLYE_SCAFFOLD=""
RAVEN_SCAFFOLD=""
SHASTA_SCAFFOLD=""

# Parse additional command-line parameters.
# We expect options in the form: --key value
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fastq) FASTQ_PATTERN="$2"; shift ;;
        --conFiles) CON_FILES="$2"; shift ;;
        --resultDir) RESULT_DIR="$2"; shift ;;
        --refmtDNA) REF_MTDNA="$2"; shift ;;
        --firstGene) FIRST_GENE="$2"; shift ;;
        --orthoDir) ORTHO_DIR="$2"; shift ;;
        --sample_id) SAMPLE_ID="$2"; shift ;;
        --genomeSize) GENOME_SIZE="$2"; shift ;;
        --reference_genome) REFERENCE_GENOME="$2"; shift ;;
        --firstA) FIRST_A="$2"; shift ;;
        --secondA) SECOND_A="$2"; shift ;;
        --thirdA) THIRD_A="$2"; shift ;;
        --fourthA) FOURTH_A="$2"; shift ;;
        --fifthA) FIFTH_A="$2"; shift ;;
        --canuScaffold) CANU_SCAFFOLD="$2"; shift ;;
        --wtdbg2Scaffold) WTDBG2_SCAFFOLD="$2"; shift ;;
        --flyeScaffold) FLYE_SCAFFOLD="$2"; shift ;;
        --ravenScaffold) RAVEN_SCAFFOLD="$2"; shift ;;
        --shastaScaffold) SHASTA_SCAFFOLD="$2"; shift ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
    shift
done

# Build the Nextflow run command.
# The "-entry" option selects the workflow (sub-workflow) to run.
CMD="nextflow run pipeline.nf -resume -entry ${WORKFLOW}"

# Append common parameters.
CMD+=" --fastq \"$FASTQ_PATTERN\" --conFiles \"$CON_FILES\" --resultDir \"$RESULT_DIR\""
CMD+=" --refmtDNA \"$REF_MTDNA\" --firstGene \"$FIRST_GENE\" --orthoDir \"$ORTHO_DIR\""
CMD+=" --sample_id \"$SAMPLE_ID\" --genomeSize \"$GENOME_SIZE\" --reference_genome \"$REFERENCE_GENOME\""
CMD+=" --firstA \"$FIRST_A\" --secondA \"$SECOND_A\" --thirdA \"$THIRD_A\" --fourthA \"$FOURTH_A\" --fifthA \"$FIFTH_A\""
CMD+=" --canuScaffold \"$CANU_SCAFFOLD\" --wtdbg2Scaffold \"$WTDBG2_SCAFFOLD\" --flyeScaffold \"$FLYE_SCAFFOLD\" --ravenScaffold \"$RAVEN_SCAFFOLD\" --shastaScaffold \"$SHASTA_SCAFFOLD\""

echo "Running: $CMD"
eval $CMD