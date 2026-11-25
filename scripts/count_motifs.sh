#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash scripts/count_motifs.sh \
      --genome path/to/assembly.fasta \
      --output path/to/results.tsv

Options:
  --motif-dir DIR   Directory containing the motif FASTA files (default: data/motifs)
  --threads N       Number of BLAST threads (default: 1)
  --keep-db         Keep the temporary BLAST database
USAGE
}

GENOME=""
OUTPUT=""
MOTIF_DIR="data/motifs"
THREADS=1
KEEP_DB=0

require_option_value() {
    if (($# < 2)) || [[ -z "${2:-}" ]]; then
        echo "ERROR: option $1 requires a value." >&2
        usage >&2
        exit 1
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)
            require_option_value "$@"
            GENOME="$2"
            shift 2
            ;;
        --output)
            require_option_value "$@"
            OUTPUT="$2"
            shift 2
            ;;
        --motif-dir)
            require_option_value "$@"
            MOTIF_DIR="$2"
            shift 2
            ;;
        --threads)
            require_option_value "$@"
            THREADS="$2"
            shift 2
            ;;
        --keep-db)
            KEEP_DB=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

[[ -n "$GENOME" && -n "$OUTPUT" ]] || {
    echo "ERROR: --genome and --output are required." >&2
    usage >&2
    exit 1
}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: --threads must be a positive integer." >&2
    exit 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}/data}"

# Absolute FASTA paths are expected when the caller supplies assemblies in
# WORK_ROOT. Relative paths are resolved from the current directory first,
# then DATA_ROOT, then the repository root.
resolve_input_file() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        [[ -s "$path" ]] || return 1
        printf '%s\n' "$path"
    elif [[ -s "$path" ]]; then
        (cd "$(dirname "$path")" && printf '%s/%s\n' "$PWD" "$(basename "$path")")
    elif [[ -s "${DATA_ROOT}/${path}" ]]; then
        printf '%s\n' "${DATA_ROOT}/${path}"
    elif [[ -s "${REPO_ROOT}/${path}" ]]; then
        printf '%s\n' "${REPO_ROOT}/${path}"
    else
        return 1
    fi
}

resolve_directory() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        [[ -d "$path" ]] || return 1
        printf '%s\n' "$path"
    elif [[ -d "$path" ]]; then
        (cd "$path" && printf '%s\n' "$PWD")
    elif [[ -d "${REPO_ROOT}/${path}" ]]; then
        printf '%s\n' "${REPO_ROOT}/${path}"
    elif [[ -d "${DATA_ROOT}/${path}" ]]; then
        printf '%s\n' "${DATA_ROOT}/${path}"
    else
        return 1
    fi
}

GENOME_PATH="$(resolve_input_file "$GENOME")" || {
    echo "ERROR: genome FASTA not found: ${GENOME}" >&2
    exit 1
}
MOTIF_PATH="$(resolve_directory "$MOTIF_DIR")" || {
    echo "ERROR: motif directory not found: ${MOTIF_DIR}" >&2
    exit 1
}

# Relative output paths are always placed under the repository root; absolute
# paths (such as WORK_ROOT from config.sh) are preserved as provided.
if [[ "$OUTPUT" = /* ]]; then
    OUTPUT_PATH="$OUTPUT"
else
    OUTPUT_PATH="${REPO_ROOT}/${OUTPUT}"
fi

for motif_file in VAR.patter.LARS.fasta Stevor.motif.fasta Rifin.motif.fasta; do
    [[ -s "${MOTIF_PATH}/${motif_file}" ]] || {
        echo "ERROR: required motif FASTA is missing or empty: ${MOTIF_PATH}/${motif_file}" >&2
        exit 1
    }
done

command -v makeblastdb >/dev/null 2>&1 || {
    echo "ERROR: makeblastdb not found in PATH." >&2
    exit 1
}
command -v tblastn >/dev/null 2>&1 || {
    echo "ERROR: tblastn not found in PATH." >&2
    exit 1
}

mkdir -p "$(dirname "$OUTPUT_PATH")"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/motif_scan.XXXXXX")"
cleanup() {
    if [[ "$KEEP_DB" -eq 0 ]]; then
        rm -rf "$TMP_DIR"
    else
        echo "BLAST database retained at: $TMP_DIR"
    fi
}
trap cleanup EXIT

DB_PREFIX="${TMP_DIR}/genome"
LARS_HITS="${TMP_DIR}/lars.tsv"
STEVOR_HITS="${TMP_DIR}/stevor.tsv"
RIFIN_HITS="${TMP_DIR}/rifin.tsv"

printf 'Building BLAST database from: %s\n' "$GENOME_PATH"
makeblastdb -in "$GENOME_PATH" -dbtype nucl -out "$DB_PREFIX" >/dev/null

tblastn \
    -db "$DB_PREFIX" \
    -query "${MOTIF_PATH}/VAR.patter.LARS.fasta" \
    -outfmt '6 qseqid sseqid pident length' \
    -num_threads "$THREADS" \
    > "$LARS_HITS"

tblastn \
    -db "$DB_PREFIX" \
    -query "${MOTIF_PATH}/Stevor.motif.fasta" \
    -outfmt '6 qseqid sseqid pident length' \
    -num_threads "$THREADS" \
    > "$STEVOR_HITS"

tblastn \
    -db "$DB_PREFIX" \
    -query "${MOTIF_PATH}/Rifin.motif.fasta" \
    -outfmt '6 qseqid sseqid pident length' \
    -num_threads "$THREADS" \
    > "$RIFIN_HITS"

lars_hits="$(awk '$3 > 90 && $4 > 10 {n++} END {print n+0}' "$LARS_HITS")"
lars_contigs="$(awk '$3 > 90 && $4 > 10 {print $2}' "$LARS_HITS" | sort -u | awk 'END {print NR+0}')"

stevor_hits="$(awk '$3 > 90 {n++} END {print n+0}' "$STEVOR_HITS")"
stevor_contigs="$(awk '$3 > 90 {print $2}' "$STEVOR_HITS" | sort -u | awk 'END {print NR+0}')"

rifin_all="$(awk 'END {print NR+0}' "$RIFIN_HITS")"
rifin_contigs="$(cut -f2 "$RIFIN_HITS" | sort -u | awk 'END {print NR+0}')"
rifin_motif2="$(awk '$1 == "rifin.motif2" {n++} END {print n+0}' "$RIFIN_HITS")"

{
    printf 'metric\tvalue\n'
    printf 'genome\t%s\n' "$GENOME_PATH"
    printf 'LARS_motifs\t%s\n' "$lars_hits"
    printf 'LARS_motifs_on_different_contigs\t%s\n' "$lars_contigs"
    printf 'Stevor_motifs\t%s\n' "$stevor_hits"
    printf 'Stevor_on_different_contigs\t%s\n' "$stevor_contigs"
    printf 'Rifin_motif2\t%s\n' "$rifin_motif2"
    printf 'Rifin_all_hits\t%s\n' "$rifin_all"
    printf 'Rifin_on_different_contigs\t%s\n' "$rifin_contigs"
} > "$OUTPUT_PATH"

printf 'Results written to: %s\n' "$OUTPUT_PATH"
cat "$OUTPUT_PATH"
