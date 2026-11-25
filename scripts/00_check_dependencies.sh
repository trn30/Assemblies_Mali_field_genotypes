#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${REPO_ROOT}/config/config.sh"

usage() {
    cat <<'USAGE'
Usage:
  bash scripts/00_check_dependencies.sh [--stage STAGE]

Stages:
  main        Default workflow: read preparation, Canu and PacBio-only ILRA
  pilon       Optional standalone Pilon run in 04_pilon_main.sh
  illumina    ILRA with Illumina correction (ILRA -p yes)
  coverage    Coverage/downsampling experiment (Canu + seqtk)
  3d7         3D7 control experiment (Canu + ILRA/Illumina)
  telomeres   TIDK search and R/ggplot2 plots
  motifs      BLAST+ motif-counting analysis
  frameshifts R frameshift analysis
  resistance R known-resistance-mutation analysis
  all         Check dependencies for every stage above

The default stage is "main", matching scripts/run_main.sh.
USAGE
}

STAGE="main"
while (($#)); do
    case "$1" in
        --stage)
            if (($# < 2)); then
                echo "ERROR: --stage requires a value." >&2
                usage >&2
                exit 2
            fi
            STAGE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$STAGE" in
    assembly) STAGE="main" ;;
    hifi-pilon) STAGE="pilon" ;;
    resistance-mutations) STAGE="resistance" ;;
esac

VALID_STAGES=(main pilon illumina coverage 3d7 telomeres motifs frameshifts resistance all)
if [[ ! " ${VALID_STAGES[*]} " =~ " ${STAGE} " ]]; then
    echo "ERROR: unknown stage '${STAGE}'." >&2
    usage >&2
    exit 2
fi

MISSING=()
WARNINGS=()
OK_COUNT=0

ok() {
    printf '  [OK] %s\n' "$1"
    ((OK_COUNT += 1))
}

missing() {
    printf '  [MISSING] %s\n' "$1" >&2
    MISSING+=("$1")
}

warning() {
    printf '  [WARNING] %s\n' "$1" >&2
    WARNINGS+=("$1")
}

check_command() {
    local tool="$1"
    local context="${2:-PATH}"
    if command -v "$tool" >/dev/null 2>&1; then
        ok "${tool} (${context})"
    else
        missing "${tool} (${context})"
    fi
}

check_file() {
    local path="$1"
    local label="${2:-$1}"
    if [[ -f "$path" && -r "$path" ]]; then
        ok "${label}: ${path}"
    else
        missing "${label}: ${path}"
    fi
}

check_executable_file() {
    local path="$1"
    local label="${2:-executable script}"
    if [[ -f "$path" && -x "$path" ]]; then
        ok "${label}: ${path}"
    else
        missing "${label} must be an executable file (set an executable path in config/config.sh): ${path}"
    fi
}


check_conda_env() {
    local env_name="$1"

    if ! command -v conda >/dev/null 2>&1; then
        missing "Conda environment '${env_name}' (conda is not in PATH)"
        return 1
    fi

    if conda run -n "$env_name" bash -c 'exit 0' >/dev/null 2>&1; then
        ok "Conda environment: ${env_name}"
        return 0
    fi

    missing "Conda environment '${env_name}' (create it or update the corresponding environment variable in config/config.sh)"
    return 1
}

check_env_command() {
    local env_name="$1"
    local tool="$2"
    local source_file="${3:-}"

    if ! command -v conda >/dev/null 2>&1; then
        missing "${tool} (Conda environment ${env_name}; conda unavailable)"
        return
    fi

    if conda run -n "$env_name" bash -c '
        source_file="$1"
        tool="$2"
        if [[ -n "$source_file" ]]; then
            [[ -r "$source_file" ]] || exit 2
            # ILRA path_to_source exports the additional tool paths used by ILRA.
            source "$source_file"
        fi
        command -v "$tool" >/dev/null 2>&1
    ' _ "$source_file" "$tool" >/dev/null 2>&1; then
        ok "${tool} (Conda environment ${env_name})"
    else
        missing "${tool} (Conda environment ${env_name})"
    fi
}

check_conda_env_commands() {
    local env_name="$1"
    local source_file="$2"
    shift 2

    if check_conda_env "$env_name"; then
        local tool
        for tool in "$@"; do
            check_env_command "$env_name" "$tool" "$source_file"
        done
    else
        for tool in "$@"; do
            missing "${tool} (not checked because Conda environment ${env_name} is unavailable)"
        done
    fi
}

check_ilra_configuration() {
    check_executable_file "$ILRA_SH" "ILRA executable"

    if [[ -n "${ILRA_PATH:-}" ]]; then
        check_file "$ILRA_PATH" "ILRA path_to_source"
    else
        warning "ILRA_PATH is empty. ILRA may still work if all its dependencies are already available in ILRA_ENV."
    fi
}

check_references() {
    check_file "$REFERENCE_FASTA" "Reference FASTA"
    check_file "$REFERENCE_GFF" "Reference GFF"
}

check_pilon_stage() {
    printf '\nPilon polishing dependencies:\n'
    check_file "$SAMPLE_TABLE" "Sample table"
    check_file "$REFERENCE_FASTA" "Reference FASTA"
    check_file "$PILON_JAR" "Pilon JAR"
    check_conda_env_commands "$ASSEMBLY_ENV" "" bwa samtools java
}

check_main_stage() {
    printf '\nMain PacBio HiFi workflow dependencies:\n'
    check_file "$SAMPLE_TABLE" "Sample table"
    check_file "$CANU_SPEC" "Canu specification"
    check_references
    check_command fastqc "global PATH (01_prepare_reads.sh)"
    check_conda_env_commands "$ASSEMBLY_ENV" "" canu
    check_conda_env "$ILRA_ENV" || true
    check_ilra_configuration

    # 01_prepare_reads.sh only needs bam2fastq when it finds BAM files
    # that have not yet been converted to the expected FASTQ output.
    local raw_dir bam needs_bam2fastq=0 found_bam=0
    for raw_dir in \
        "${DATA_ROOT}/Clones_ship1/RAW" \
        "${DATA_ROOT}/clones_ship2/RAW" \
        "${DATA_ROOT}/Clones_ship3/RAW"; do
        if [[ ! -d "$raw_dir" ]]; then
            missing "Read directory required by 01_prepare_reads.sh: ${raw_dir}"
            continue
        fi
        shopt -s nullglob
        local bam_files=("${raw_dir}"/*.bam)
        shopt -u nullglob
        for bam in "${bam_files[@]}"; do
            found_bam=1
            if [[ ! -e "${bam}.fastq" && ! -e "${bam}.fastq.fastq" ]]; then
                needs_bam2fastq=1
                break
            fi
        done
    done

    if ((needs_bam2fastq)); then
        check_command bam2fastq "global PATH (BAM inputs need conversion)"
    elif ((found_bam)); then
        ok "bam2fastq not required: FASTQ outputs already exist for detected BAM files"
    else
        ok "bam2fastq not required: no BAM inputs detected in the read directories"
    fi

    # In the current run_main.sh, RUN_ILLUMINA_CORRECTION=1 invokes
    # 04_pilon_main.sh. Check those dependencies if that option is enabled.
    if [[ "${RUN_ILLUMINA_CORRECTION:-0}" == "1" ]]; then
        warning "In the current run_main.sh, RUN_ILLUMINA_CORRECTION=1 invokes 04_pilon_main.sh (Pilon), not 05_illumina_correction.sh."
        check_pilon_stage
    fi
}

check_illumina_stage() {
    printf '\nIllumina-assisted ILRA dependencies:\n'
    check_file "$SAMPLE_TABLE" "Sample table"
    check_references
    check_ilra_configuration
    check_conda_env_commands "$ILRA_ENV" "$ILRA_PATH" bowtie2-build bowtie2 samtools pilon
}

check_coverage_stage() {
    printf '\nCoverage/downsampling dependencies:\n'
    check_file "$SAMPLE_TABLE" "Sample table"
    check_file "$DOWNSAMPLING_TABLE" "Downsampling table"
    check_file "$CANU_SPEC" "Canu specification"
    check_conda_env_commands "$ASSEMBLY_ENV" "" canu seqtk
}

check_3d7_stage() {
    printf '\n3D7 control experiment dependencies:\n'
    check_references
    check_file "$CANU_SPEC" "Canu specification"
    check_file "$THREE_D7_PACBIO_FASTQ" "3D7 PacBio FASTQ"
    if compgen -G "${THREE_D7_ILLUMINA_PREFIX}*" >/dev/null; then
        ok "3D7 Illumina files for prefix: ${THREE_D7_ILLUMINA_PREFIX}"
    else
        missing "3D7 Illumina files for prefix: ${THREE_D7_ILLUMINA_PREFIX}*"
    fi
    check_conda_env_commands "$ASSEMBLY_ENV" "" canu
    check_ilra_configuration
    check_conda_env_commands "$ILRA_ENV" "$ILRA_PATH" bowtie2-build bowtie2 samtools pilon
}

check_telomere_stage() {
    printf '\nTelomere analysis dependencies:\n'
    check_file "$SAMPLE_TABLE" "Sample table"
    check_file "$REFERENCE_FASTA" "Reference FASTA"
    check_command tidk
    check_command samtools
    check_command Rscript
    if command -v Rscript >/dev/null 2>&1; then
        if Rscript --vanilla -e 'if (!requireNamespace("ggplot2", quietly=TRUE)) quit(status=1)' >/dev/null 2>&1; then
            ok "R package: ggplot2"
        else
            missing "R package: ggplot2 (install from CRAN in the R environment used by Rscript)"
        fi
    fi
}

check_motif_stage() {
    printf '\nMotif analysis dependencies:\n'
    check_file "$SAMPLE_TABLE" "Sample table"
    check_file "${REPO_ROOT}/data/motifs/VAR.patter.LARS.fasta" "LARS motif FASTA"
    check_file "${REPO_ROOT}/data/motifs/Stevor.motif.fasta" "Stevor motif FASTA"
    check_file "${REPO_ROOT}/data/motifs/Rifin.motif.fasta" "Rifin motif FASTA"
    check_command makeblastdb
    check_command tblastn
}

check_frameshift_stage() {
    printf '\nFrameshift analysis dependencies:\n'
    check_file "$FRAMESHIFT_TABLE" "Frameshift sample table"
    check_file "$REFERENCE_GFF" "Reference GFF"
    check_command Rscript
}

check_resistance_stage() {
    printf '\nResistance mutation analysis dependencies:\n'
    check_command Rscript
    if command -v Rscript >/dev/null 2>&1; then
        if Rscript --vanilla -e '
            pkgs <- c("Biostrings", "argparse", "stringr", "dplyr", "tidyr")
            missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
            if (length(missing)) {
                cat(paste(missing, collapse=" "), "\n", file=stderr())
                quit(status=1)
            }
        ' >/dev/null 2>&1; then
            ok "R packages: Biostrings, argparse, stringr, dplyr, tidyr"
        else
            missing "One or more R packages required by 08_resistant_mutation_identification.R (Biostrings, argparse, stringr, dplyr, tidyr)"
        fi
    fi
}

case "$STAGE" in
    main)
        check_main_stage
        ;;
    pilon)
        check_pilon_stage
        ;;
    illumina)
        check_illumina_stage
        ;;
    coverage)
        check_coverage_stage
        ;;
    3d7)
        check_3d7_stage
        ;;
    telomeres)
        check_telomere_stage
        ;;
    motifs)
        check_motif_stage
        ;;
    frameshifts)
        check_frameshift_stage
        ;;
    resistance)
        check_resistance_stage
        ;;
    all)
        check_main_stage
        check_pilon_stage
        check_illumina_stage
        check_coverage_stage
        check_3d7_stage
        check_telomere_stage
        check_motif_stage
        check_frameshift_stage
        check_resistance_stage
        ;;
esac

printf '\nDependency check summary for stage: %s\n' "$STAGE"
printf '  Available checks passed: %d\n' "$OK_COUNT"
printf '  Missing requirements: %d\n' "${#MISSING[@]}"
printf '  Warnings: %d\n' "${#WARNINGS[@]}"

if ((${#MISSING[@]} > 0)); then
    printf '\nERROR: one or more required dependencies or files are missing.\n' >&2
    echo "Install/fix the missing items above, then run this check again." >&2
    exit 1
fi

echo "Dependency check passed for stage '${STAGE}'."
