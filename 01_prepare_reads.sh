#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"
log_dir

for shipment in Clones_ship1 clones_ship2 Clones_ship3; do
    raw_dir="${DATA_ROOT}/${shipment}/RAW"
    fastqc_dir="${raw_dir}/Fastqc"

    require_dir "$raw_dir"
    mkdir -p "$fastqc_dir"

    shopt -s nullglob
    bams=( "$raw_dir"/*.bam )
    shopt -u nullglob

    for bam in "${bams[@]}"; do
        output_prefix="${bam}.fastq"

        if [[ ! -e "$output_prefix" && ! -e "${output_prefix}.fastq" ]]; then
            echo "==> BAM -> FASTQ: $(basename "$bam")"
            bam2fastq -o "$output_prefix" -j "$THREADS" -u "$bam"
        else
            echo "==> FASTQ already present for $(basename "$bam")"
        fi
    done

    shopt -s nullglob
    fastqs=( "$raw_dir"/*.fastq )
    shopt -u nullglob

    if ((${#fastqs[@]} > 0)); then
        fastqc -t "$THREADS" "${fastqs[@]}" -o "$fastqc_dir" \
            > "${WORK_ROOT}/logs/${shipment}_fastqc.log" 2>&1
    fi
done
