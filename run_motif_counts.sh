#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"

while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    [[ -z "$isolate" || "$isolate" == \#* ]] && continue

    assembly="${WORK_ROOT}/pilon/${isolate}/${isolate}_pilon.fasta"

    if [[ ! -s "$assembly" ]]; then
        echo "WARNING: final Pilon assembly not found for ${isolate}: ${assembly}" >&2
        continue
    fi

    bash "${REPO_ROOT}/scripts/count_motifs.sh" \
        --genome "$assembly" \
        --output "results/motifs/${isolate}.tsv"
done < "$SAMPLE_TABLE"
