#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"
log_dir
activate_conda_env "$ASSEMBLY_ENV"

mkdir -p "${WORK_ROOT}/canu"

while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    [[ -z "$isolate" || "$isolate" == \#* ]] && continue

    fq="${DATA_ROOT}/${raw_fastq}"
    out_dir="${WORK_ROOT}/canu/${isolate}"
    result_fasta="${out_dir}/${isolate}.contigs.fasta"

    require_file "$fq"
    mkdir -p "$out_dir"

    if [[ -s "$result_fasta" ]]; then
        echo "==> already assembled: ${isolate}"
        continue
    fi

    echo "==> Canu HiFi assembly: ${isolate}"

    canu \
        corMaxEvidenceErate="${CANU_COR_MAX_EVIDENCE_ERATE}" \
        maxThreads="${THREADS}" \
        genomeSize="${CANU_GENOME_SIZE}" \
        useGrid="${CANU_USE_GRID}" \
        readSamplingCoverage="${CANU_READ_COVERAGE}" \
        -p "${isolate}" \
        -d "$out_dir" \
        -s "$CANU_SPEC" \
        -pacbio-hifi "$fq" \
        > "${WORK_ROOT}/logs/${isolate}_canu.log" 2>&1

    generated="$(find "$out_dir" -maxdepth 1 -type f -name "*.contigs.fasta" -print -quit)"
    [[ -n "$generated" ]] || {
        echo "ERROR: Canu did not produce a contigs FASTA for ${isolate}" >&2
        exit 1
    }

    if [[ "$generated" != "$result_fasta" ]]; then
        mv "$generated" "$result_fasta"
    fi
done < "$SAMPLE_TABLE"
