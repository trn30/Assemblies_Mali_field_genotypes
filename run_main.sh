#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

bash "${SCRIPT_DIR}/00_check_dependencies.sh"
bash "${SCRIPT_DIR}/01_prepare_reads.sh"
bash "${SCRIPT_DIR}/02_canu_main.sh"
bash "${SCRIPT_DIR}/03_ilra_main.sh"
bash "${SCRIPT_DIR}/04_pilon_main.sh"
