#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${REPO_ROOT}/config/config.sh"

require_file "$SAMPLE_TABLE"
MOTIF_OUTPUT_DIR="${WORK_ROOT}/motifs"
mkdir -p "$MOTIF_OUTPUT_DIR"

while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    isolate="${isolate//$'\r'/}"
    [[ -z "$isolate" || "$isolate" == \#* || "$isolate" == "isolate" ]] && continue

    # The Pilon script writes <WORK_ROOT>/pilon/<isolate>/<isolate>_pilon.fasta.
    assembly="${WORK_ROOT}/pilon/${isolate}/${isolate}_pilon.fasta"
    output="${MOTIF_OUTPUT_DIR}/${isolate}.tsv"

    if [[ ! -s "$assembly" ]]; then
        echo "WARNING: final Pilon assembly not found for ${isolate}: ${assembly}" >&2
        continue
    fi

    bash "${SCRIPT_DIR}/count_motifs.sh" \
        --genome "$assembly" \
        --output "$output" \
        --motif-dir "${REPO_ROOT}/data/motifs" \
        --threads "$THREADS"
done < "$SAMPLE_TABLE"
