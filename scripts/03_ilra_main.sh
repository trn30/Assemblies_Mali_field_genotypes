#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"
log_dir

require_file "$REFERENCE_FASTA"
require_file "$REFERENCE_GFF"
require_file "$ILRA_SH"

activate_conda_env "$ILRA_ENV"
[[ -f "$ILRA_PATH" ]] && source "$ILRA_PATH"

while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    [[ -z "$isolate" || "$isolate" == \#* ]] && continue

    assembly="${WORK_ROOT}/canu/${isolate}/${isolate}.contigs.fasta"
    out_dir="${WORK_ROOT}/ilra/${isolate}"
    mkdir -p "$out_dir"

    require_file "$assembly"

    if find "$out_dir" -type f -name "*ILRA.fasta" | grep -q .; then
        echo "==> already processed by ILRA: ${isolate}"
        continue
    fi

    echo "==> ILRA (PacBio only): ${isolate}"

    "$ILRA_SH" \
        -a "$assembly" \
        -o "$out_dir" \
        -n "$isolate" \
        -C no \
        -r "$REFERENCE_FASTA" \
        -t "$THREADS" \
        -g "$REFERENCE_GFF" \
        -L pb \
        -K yes \
        -M 200g \
        -p yes \
        -l \
        > "${WORK_ROOT}/logs/${isolate}_ilra.log" 2>&1
done < "$SAMPLE_TABLE"
