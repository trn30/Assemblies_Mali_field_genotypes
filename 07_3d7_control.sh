#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"
log_dir

require_file "$REFERENCE_FASTA"
require_file "$REFERENCE_GFF"
require_file "$THREE_D7_PACBIO_FASTQ"
require_prefix "$THREE_D7_ILLUMINA_PREFIX"
require_file "$CANU_SPEC"

activate_conda_env "$ASSEMBLY_ENV"

control_root="${WORK_ROOT}/3D7_control"
canu_root="${control_root}/canu"
ilra_root="${control_root}/ilra"

mkdir -p "$canu_root" "$ilra_root"

for coverage in 25 50 100 200 450; do
    out_dir="${canu_root}/${coverage}X"
    mkdir -p "$out_dir"

    if [[ ! -s "${out_dir}/3D7_${coverage}X.contigs.fasta" ]]; then
        canu \
            corMaxEvidenceErate=0.15 \
            maxThreads="${THREADS}" \
            genomeSize="${CANU_GENOME_SIZE}" \
            useGrid=0 \
            readSamplingCoverage="${coverage}" \
            stopOnLowCoverage=8 \
            -p "3D7_${coverage}X" \
            -d "$out_dir" \
            -s "$CANU_SPEC" \
            -pacbio "$THREE_D7_PACBIO_FASTQ" \
            > "${WORK_ROOT}/logs/3D7_${coverage}X_canu.log" 2>&1
    fi
done

activate_conda_env "$ILRA_ENV"
[[ -f "$ILRA_PATH" ]] && source "$ILRA_PATH"

for coverage in 200 450; do
    assembly="${canu_root}/${coverage}X/3D7_${coverage}X.contigs.fasta"
    out_dir="${ilra_root}/${coverage}X"

    require_file "$assembly"
    mkdir -p "$out_dir"

    if ! find "$out_dir" -type f -name "*ILRA.fasta" | grep -q .; then
        "$ILRA_SH" \
            -a "$assembly" \
            -o "$out_dir" \
            -n "3D7_${coverage}X" \
            -C yes \
            -I "$THREE_D7_ILLUMINA_PREFIX" \
            -r "$REFERENCE_FASTA" \
            -t "$THREADS" \
            -g "$REFERENCE_GFF" \
            -L pb \
            -K yes \
            -Mj 200g \
            -p yes \
            > "${WORK_ROOT}/logs/3D7_${coverage}X_ilra.log" 2>&1
    fi
done
