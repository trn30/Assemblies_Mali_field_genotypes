#!/usr/bin/env bash
set -euo pipefail

# Run TIDK for 3D7 and each of the eight isolates, then run the R plotting and
# counting script. This script is called at the end of run_main.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${REPO_ROOT}/config/config.sh"

for tool in tidk samtools Rscript; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: required program '$tool' is not in PATH." >&2
        exit 1
    }
done

require_file "$REFERENCE_FASTA"
require_file "$SAMPLE_TABLE"

TEL_ROOT="${WORK_ROOT}/telomeres"
TIDK_DIR="${TEL_ROOT}/tidk"
FIGURE_DIR="${TEL_ROOT}/figures"
INPUT_TABLE="${TEL_ROOT}/telomere_inputs.tsv"

mkdir -p "$TIDK_DIR" "$FIGURE_DIR"
printf "Sample\tFasta\tTidkTSV\n" > "$INPUT_TABLE"

run_tidk() {
    local sample="$1"
    local fasta="$2"
    local fai="${fasta}.fai"
    local tsv="${TIDK_DIR}/${sample}_telomeric_repeat_windows.tsv"

    [[ -s "$fasta" ]] || {
        echo "ERROR: FASTA for ${sample} not found: ${fasta}" >&2
        exit 1
    }

    # Build/update the FASTA index required by the R plotting script.
    if [[ ! -s "$fai" || "$fasta" -nt "$fai" ]]; then
        echo "==> samtools faidx: ${sample}"
        samtools faidx "$fasta"
    fi

    echo "==> TIDK search: ${sample}"
    tidk search \
        -s TTTAGGG \
        -w 10000 \
        -d "$TIDK_DIR" \
        -o "$sample" \
        "$fasta"

    [[ -s "$tsv" ]] || {
        echo "ERROR: expected TIDK output not found: ${tsv}" >&2
        echo "Check the TIDK version and the output filename it generates." >&2
        exit 1
    }

    printf "%s\t%s\t%s\n" "$sample" "$fasta" "$tsv" >> "$INPUT_TABLE"
}

# The reference is processed first.
run_tidk "3D7" "$REFERENCE_FASTA"

# Process all isolates defined in config/samples.tsv. Assembly FASTAs are
# expected under $WORK_ROOT/ilra/<isolate>/ and have the suffix *ILRA.fasta.
while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    [[ -z "$isolate" || "$isolate" == \#* ]] && continue

    sample_dir="${WORK_ROOT}/ilra/${shipment}/${isolate}"
    if [[ ! -d "$sample_dir" ]]; then
        # Fallback for installations with a flat ILRA output layout.
        sample_dir="${WORK_ROOT}/ilra/${isolate}"
    fi

    [[ -d "$sample_dir" ]] || {
        echo "ERROR: ILRA output directory not found for ${isolate}: ${sample_dir}" >&2
        echo "Run the assembly/ILRA stages first." >&2
        exit 1
    }

    mapfile -t candidates < <(find "$sample_dir" -type f -name "*ILRA.fasta" -print | sort)
    if (("${#candidates[@]}" == 0)); then
        echo "ERROR: no *ILRA.fasta assembly found for ${isolate} under ${sample_dir}" >&2
        exit 1
    elif (("${#candidates[@]}" > 1)); then
        echo "ERROR: multiple ILRA FASTAs found for ${isolate}:" >&2
        printf '  %s\n' "${candidates[@]}" >&2
        echo "Choose the intended final assembly or narrow the find pattern in this script." >&2
        exit 1
    fi

    run_tidk "$isolate" "${candidates[0]}"
done < "$SAMPLE_TABLE"

echo "==> Plotting and counting telomeric repeat windows in R"
Rscript "${SCRIPT_DIR}/plot_telomeres_repetitions.R" \
    --input-table "$INPUT_TABLE" \
    --output-dir "$FIGURE_DIR"

echo "Telomere analysis finished."
echo "TIDK output: ${TIDK_DIR}"
echo "Figures and counts: ${FIGURE_DIR}"
