#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash scripts/count_motifs.sh \
      --genome path/to/assembly.fasta \
      --output path/to/results.tsv

Optional:
      --motif-dir data/motifs
      --threads 8
      --keep-db
EOF
}

GENOME=""
OUTPUT=""
MOTIF_DIR="data/motifs"
THREADS=1
KEEP_DB=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome) GENOME="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --motif-dir) MOTIF_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --keep-db) KEEP_DB=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

[[ -n "$GENOME" && -n "$OUTPUT" ]] || {
    echo "Genome and output are required." >&2
    usage >&2
    exit 1
}

# Resolve paths relative to the repository, but never require the user's
# original machine paths.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

case "$GENOME" in
    /*) echo "ERROR: use a path relative to the repository or DATA_ROOT." >&2; exit 1 ;;
esac
case "$MOTIF_DIR" in
    /*) echo "ERROR: use a path relative to the repository or DATA_ROOT." >&2; exit 1 ;;
esac
case "$OUTPUT" in
    /*) echo "ERROR: use a path relative to the repository." >&2; exit 1 ;;
esac

DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}}"
GENOME_PATH="${DATA_ROOT}/${GENOME}"
MOTIF_PATH="${REPO_ROOT}/${MOTIF_DIR}"
OUTPUT_PATH="${REPO_ROOT}/${OUTPUT}"

[[ -s "$GENOME_PATH" ]] || { echo "ERROR: genome not found: $GENOME_PATH" >&2; exit 1; }
[[ -s "${MOTIF_PATH}/VAR.patter.LARS.fasta" ]] || { echo "ERROR: missing LARS motif FASTA" >&2; exit 1; }
[[ -s "${MOTIF_PATH}/Stevor.motif.fasta" ]] || { echo "ERROR: missing Stevor motif FASTA" >&2; exit 1; }
[[ -s "${MOTIF_PATH}/Rifin.motif.fasta" ]] || { echo "ERROR: missing Rifin motif FASTA" >&2; exit 1; }

command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found." >&2; exit 1; }
command -v tblastn >/dev/null 2>&1 || { echo "ERROR: tblastn not found." >&2; exit 1; }

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

echo "Building BLAST database from: $GENOME_PATH"
makeblastdb -in "$GENOME_PATH" -dbtype nucl -out "$DB_PREFIX" >/dev/null

LARS_HITS="${TMP_DIR}/lars.tsv"
STEVOR_HITS="${TMP_DIR}/stevor.tsv"
RIFIN_HITS="${TMP_DIR}/rifin.tsv"

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
    printf "metric\tvalue\n"
    printf "genome\t%s\n" "$GENOME"
    printf "LARS_motifs\t%s\n" "$lars_hits"
    printf "LARS_motifs_on_different_contigs\t%s\n" "$lars_contigs"
    printf "Stevor_motifs\t%s\n" "$stevor_hits"
    printf "Stevor_on_different_contigs\t%s\n" "$stevor_contigs"
    printf "Rifin_motif2\t%s\n" "$rifin_motif2"
    printf "Rifin_all_hits\t%s\n" "$rifin_all"
    printf "Rifin_on_different_contigs\t%s\n" "$rifin_contigs"
} > "$OUTPUT_PATH"

echo "Results written to: $OUTPUT_PATH"
cat "$OUTPUT_PATH"
