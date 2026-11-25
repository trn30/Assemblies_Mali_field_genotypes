#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"
log_dir
activate_conda_env "$ASSEMBLY_ENV"

mkdir -p "${WORK_ROOT}/coverage_experiment/reads"
mkdir -p "${WORK_ROOT}/coverage_experiment/canu"

while IFS=$'\t' read -r isolate coverage fraction seed; do
    [[ -z "$isolate" || "$isolate" == \#* ]] && continue

    raw_fastq="$(sample_fastq_path "$isolate")"
    shipment="$(sample_shipment "$isolate")"

    require_file "${DATA_ROOT}/${raw_fastq}"

    sampled_fq="${WORK_ROOT}/coverage_experiment/reads/${isolate}_${coverage}X.fastq"
    out_dir="${WORK_ROOT}/coverage_experiment/canu/${isolate}_${coverage}X"

    mkdir -p "$out_dir"

    if [[ ! -s "$sampled_fq" ]]; then
        echo "==> Downsampling ${isolate} to ${coverage}X"
        seqtk sample -s "$seed" "${DATA_ROOT}/${raw_fastq}" "$fraction" > "$sampled_fq"
    fi

    if [[ ! -s "${out_dir}/${isolate}_${coverage}X_assembly.contigs.fasta" ]]; then
        canu \
            corMaxEvidenceErate="${CANU_COR_MAX_EVIDENCE_ERATE}" \
            maxThreads="${THREADS}" \
            genomeSize="${CANU_GENOME_SIZE}" \
            useGrid="${CANU_USE_GRID}" \
            readSamplingCoverage="${coverage}" \
            -p "${isolate}_${coverage}X_assembly" \
            -d "$out_dir" \
            -s "$CANU_SPEC" \
            -pacbio-hifi "$sampled_fq" \
            > "${WORK_ROOT}/logs/${isolate}_${coverage}X_canu.log" 2>&1
    fi
done < "$DOWNSAMPLING_TABLE"

# The original analysis used the full read sets for 50X and 100X, letting
# Canu control the working coverage with readSamplingCoverage.
for isolate in BG1905 BG0801 KK2504 KK2202 BG1701 KK2403; do
    raw_fastq="$(sample_fastq_path "$isolate")"
    require_file "${DATA_ROOT}/${raw_fastq}"

    for coverage in 50 100; do
        out_dir="${WORK_ROOT}/coverage_experiment/canu/${isolate}_${coverage}X"
        mkdir -p "$out_dir"

        if [[ ! -s "${out_dir}/${isolate}_${coverage}X_assembly.contigs.fasta" ]]; then
            echo "==> Canu ${isolate} at ${coverage}X"

            canu \
                corMaxEvidenceErate="${CANU_COR_MAX_EVIDENCE_ERATE}" \
                maxThreads="${THREADS}" \
                genomeSize="${CANU_GENOME_SIZE}" \
                useGrid="${CANU_USE_GRID}" \
                readSamplingCoverage="${coverage}" \
                -p "${isolate}_${coverage}X_assembly" \
                -d "$out_dir" \
                -s "$CANU_SPEC" \
                -pacbio-hifi "${DATA_ROOT}/${raw_fastq}" \
                > "${WORK_ROOT}/logs/${isolate}_${coverage}X_canu.log" 2>&1
        fi
    done
done
