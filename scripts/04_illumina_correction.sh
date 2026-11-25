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
    [[ "$illumina_prefix" == "-" ]] && continue

    assembly="${WORK_ROOT}/canu/${isolate}/${isolate}.contigs.fasta"
    out_root="${WORK_ROOT}/illumina_correction/${isolate}"

    require_file "$assembly"
    require_prefix "${DATA_ROOT}/${illumina_prefix}"

    mkdir -p "$out_root"

    echo "==> ${isolate}: ILRA + Illumina correction, single pass"

    "$ILRA_SH" \
        -a "$assembly" \
        -o "${out_root}/x3" \
        -n "${isolate}_ILRA_pilon" \
        -C yes \
        -I "${DATA_ROOT}/${illumina_prefix}" \
        -r "$REFERENCE_FASTA" \
        -t "$THREADS" \
        -g "$REFERENCE_GFF" \
        -L pb \
        -K yes \
        -Mj 200g \
        -p yes \
        > "${WORK_ROOT}/logs/${isolate}_ilra_illumina_x3.log" 2>&1

    echo "==> ${isolate}: ILRA + Illumina correction, five iterations"

    "$ILRA_SH" \
        -a "$assembly" \
        -o "${out_root}/x5" \
        -n "${isolate}_ILRA_pilon_x5" \
        -C yes \
        -I "${DATA_ROOT}/${illumina_prefix}" \
        -r "$REFERENCE_FASTA" \
        -t "$THREADS" \
        -g "$REFERENCE_GFF" \
        -L pb \
        -K yes \
        -Mj 200g \
        -p yes \
        -i 5 \
        > "${WORK_ROOT}/logs/${isolate}_ilra_illumina_x5.log" 2>&1
done < "$SAMPLE_TABLE"
