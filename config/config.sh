#!/usr/bin/env bash

# Runtime configuration.
# Nothing here depends on a particular personal/server filesystem.
#
# Example:
#   export DATA_ROOT=/path/to/project_data
#   export REFERENCE_FASTA=/path/to/PlasmoDB-3D7.fasta
#   export REFERENCE_GFF=/path/to/PlasmoDB-3D7.gff
#   export ILRA_SH=ILRA.sh
#   export ILRA_PATH=/path/to/ILRA/path_to_source
#   export PILON_JAR=pilon.jar
#
# When DATA_ROOT is not provided, data are expected under ./data.
# Large sequencing/annotation files should normally live outside Git.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}/data}"
export WORK_ROOT="${WORK_ROOT:-${REPO_ROOT}/results}"

export REFERENCE_FASTA="${REFERENCE_FASTA:-${DATA_ROOT}/reference/PlasmoDB63-3D7.fasta}"
export REFERENCE_GFF="${REFERENCE_GFF:-${DATA_ROOT}/reference/PlasmoDB63-3D7.gff}"

export THREADS="${THREADS:-30}"
export CANU_GENOME_SIZE="${CANU_GENOME_SIZE:-24000000}"
export CANU_READ_COVERAGE="${CANU_READ_COVERAGE:-200}"
export CANU_COR_MAX_EVIDENCE_ERATE="${CANU_COR_MAX_EVIDENCE_ERATE:-0.15}"
export CANU_USE_GRID="${CANU_USE_GRID:-0}"

export ASSEMBLY_ENV="${ASSEMBLY_ENV:-assembly}"

export ILRA_SH="${ILRA_SH:-ILRA.sh}"
export ILRA_PATH="${ILRA_PATH:-}"
export ILRA_ENV="${ILRA_ENV:-ILRA_env}"

export PILON_JAR="${PILON_JAR:-pilon.jar}"
export PILON_MEMORY="${PILON_MEMORY:-100G}"

export THREE_D7_PACBIO_FASTQ="${THREE_D7_PACBIO_FASTQ:-${DATA_ROOT}/3D7/SRR31948_merged.fastq}"
export THREE_D7_ILLUMINA_PREFIX="${THREE_D7_ILLUMINA_PREFIX:-${DATA_ROOT}/3D7/ERR04_merged}"

export CANU_SPEC="${REPO_ROOT}/config/canu_hifi.spec"
export SAMPLE_TABLE="${REPO_ROOT}/config/samples.tsv"
export DOWNSAMPLING_TABLE="${REPO_ROOT}/config/downsampling.tsv"
export FRAMESHIFT_TABLE="${REPO_ROOT}/config/frameshift_samples.tsv"

require_file() {
    local f="$1"
    [[ -f "$f" ]] || {
        echo "ERROR: file not found: $f" >&2
        exit 1
    }
}

require_dir() {
    local d="$1"
    [[ -d "$d" ]] || {
        echo "ERROR: directory not found: $d" >&2
        exit 1
    }
}

require_prefix() {
    local prefix="$1"
    compgen -G "${prefix}*" >/dev/null || {
        echo "ERROR: no files found for prefix: ${prefix}*" >&2
        exit 1
    }
}

activate_conda_env() {
    local env="$1"

    command -v conda >/dev/null 2>&1 || {
        echo "ERROR: conda is not available in PATH." >&2
        exit 1
    }

    eval "$(conda shell.bash hook)"
    conda activate "$env"
}

sample_fastq_path() {
    local sample="$1"

    awk -F '\t' -v target="$sample" '
        $0 !~ /^#/ && $1 == target {print $3; exit}
    ' "$SAMPLE_TABLE"
}

sample_shipment() {
    local sample="$1"

    awk -F '\t' -v target="$sample" '
        $0 !~ /^#/ && $1 == target {print $2; exit}
    ' "$SAMPLE_TABLE"
}

sample_illumina_prefix() {
    local sample="$1"

    awk -F '\t' -v target="$sample" '
        $0 !~ /^#/ && $1 == target {print $4; exit}
    ' "$SAMPLE_TABLE"
}

log_dir() {
    mkdir -p "${WORK_ROOT}/logs"
}
