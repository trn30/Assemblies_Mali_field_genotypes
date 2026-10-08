#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"
log_dir

require_file "$REFERENCE_FASTA"
require_file "$PILON_JAR"
activate_conda_env "$ASSEMBLY_ENV"

while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    [[ -z "$isolate" || "$isolate" == \#* ]] && continue

    fq="${DATA_ROOT}/${raw_fastq}"
    ilra_fasta="$(find "${WORK_ROOT}/ilra/${isolate}" -type f -name "*ILRA.fasta" -print -quit)"

    require_file "$fq"

    if [[ -z "$ilra_fasta" ]]; then
        echo "ERROR: no ILRA FASTA found for ${isolate}" >&2
        exit 1
    fi

    out_dir="${WORK_ROOT}/pilon/${isolate}"
    align_dir="${out_dir}/alignment"
    bam="${align_dir}/${isolate}.sorted.bam"

    mkdir -p "$align_dir"

    if [[ ! -e "${ilra_fasta}.bwt" ]]; then
        bwa index "$ilra_fasta"
    fi

    if [[ ! -s "$bam" ]]; then
        echo "==> Mapping reads for Pilon: ${isolate}"

        bwa mem \
            -x pacbio \
            -t "$THREADS" \
            -B 100 \
            -M \
            -T 100 \
            "$ilra_fasta" \
            "$fq" \
            | samtools view -@ "$THREADS" -bh - \
            | samtools sort -@ "$THREADS" -o "$bam" -

        samtools index "$bam"
    fi

    output_prefix="${isolate}_pilon"

    if [[ ! -s "${out_dir}/${output_prefix}.fasta" ]]; then
        echo "==> Pilon polishing: ${isolate}"

        java -Xmx"${PILON_MEMORY}" -jar "$PILON_JAR" \
            --genome "$ilra_fasta" \
            --bam "$bam" \
            --output "$output_prefix" \
            --outdir "$out_dir" \
            --changes \
            --vcf \
            --tracks \
            --threads "$THREADS" \
            > "${WORK_ROOT}/logs/${isolate}_pilon.log" 2>&1
    fi
done < "$SAMPLE_TABLE"
