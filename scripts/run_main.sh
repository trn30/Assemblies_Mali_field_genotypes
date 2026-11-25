#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# 1. Check dependencies
bash "${SCRIPT_DIR}/00_check_dependencies.sh"

# 2. Prepare PacBio HiFi reads
bash "${SCRIPT_DIR}/01_prepare_reads.sh"

# 3. De novo assembly using Canu HiFi
bash "${SCRIPT_DIR}/02_canu_main.sh"

# 4. HiFi-only ILRA processing
bash "${SCRIPT_DIR}/03_ilra_main.sh"

# 5. pilon polishing one iteration 
if [[ "${RUN_ILLUMINA_CORRECTION:-0}" == "1" ]]; then
    bash "${SCRIPT_DIR}/04_pilon_main.sh"
fi

echo "Pipeline completed."
