#!/usr/bin/env bash
set -euo pipefail

# Run TIDK on the 3D7 reference and each ILRA assembly, then pass the
# generated FASTA/FAI and TIDK TSV paths to the plotting script.

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
printf 'Sample\tFasta\tTidkTSV\n' > "$INPUT_TABLE"

run_tidk() {
    local sample="$1"
    local fasta="$2"
    local fai="${fasta}.fai"
    local tsv="${TIDK_DIR}/${sample}_telomeric_repeat_windows.tsv"

    [[ -s "$fasta" ]] || {
        echo "ERROR: FASTA for ${sample} not found or empty: ${fasta}" >&2
        exit 1
    }

    # TIDK's R plot uses the FASTA index to draw chromosome backgrounds.
    if [[ ! -s "$fai" || "$fasta" -nt "$fai" ]]; then
        echo "==> samtools faidx: ${sample}"
        samtools faidx "$fasta"
    fi

    echo "==> TIDK search: ${sample}"
    tidk search \
        --string TTTAGGG \
        --window 10000 \
        --dir "$TIDK_DIR" \
        --output "$sample" \
        "$fasta"

    [[ -s "$tsv" ]] || {
        echo "ERROR: expected TIDK output not found or empty: ${tsv}" >&2
        echo "Expected filename: ${sample}_telomeric_repeat_windows.tsv" >&2
        exit 1
    }

    printf '%s\t%s\t%s\n' "$sample" "$fasta" "$tsv" >> "$INPUT_TABLE"
}

# Process the reference first.
run_tidk "3D7" "$REFERENCE_FASTA"

# Process all isolates in config/samples.tsv.
while IFS=$'\t' read -r isolate shipment raw_fastq illumina_prefix notes; do
    # Remove a possible Windows carriage return and skip headers/comments.
    isolate="${isolate//$'\r'/}"
    [[ -z "$isolate" || "$isolate" == \#* || "$isolate" == "isolate" ]] && continue

    sample_dir="${WORK_ROOT}/ilra/${shipment}/${isolate}"
    if [[ ! -d "$sample_dir" ]]; then
        # The main ILRA script currently writes to results/ilra/<isolate>/.
        sample_dir="${WORK_ROOT}/ilra/${isolate}"
    fi

    [[ -d "$sample_dir" ]] || {
        echo "ERROR: ILRA output directory not found for ${isolate}: ${sample_dir}" >&2
        echo "Run the assembly/ILRA stages first." >&2
        exit 1
    }

    mapfile -t candidates < <(find "$sample_dir" -type f -name '*ILRA.fasta' -print | sort)
    if ((${#candidates[@]} == 0)); then
        echo "ERROR: no *ILRA.fasta assembly found for ${isolate} under ${sample_dir}" >&2
        exit 1
    elif ((${#candidates[@]} > 1)); then
        echo "ERROR: multiple ILRA FASTAs found for ${isolate}:" >&2
        printf '  %s\n' "${candidates[@]}" >&2
        echo "Choose the intended final assembly or narrow the find pattern in this script." >&2
        exit 1
    fi

    run_tidk "$isolate" "${candidates[0]}"
done < "$SAMPLE_TABLE"

echo "==> Plotting telomeric repeat windows in R"
Rscript "${SCRIPT_DIR}/plot_telomeres_repetitions.R" \
    --input-table "$INPUT_TABLE" \
    --output-dir "$FIGURE_DIR"

echo "Telomere analysis finished."
echo "TIDK output: ${TIDK_DIR}"
echo "Figures: ${FIGURE_DIR}"
