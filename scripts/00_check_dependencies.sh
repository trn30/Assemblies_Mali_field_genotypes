#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/.." && pwd)/config/config.sh"

tools=(bam2fastq fastqc canu bwa samtools seqtk java)

for tool in "${tools[@]}"; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "OK: $tool -> $(command -v "$tool")"
    else
        echo "MISSING: $tool" >&2
    fi
done

if [[ -x "$ILRA_SH" ]]; then
    echo "OK: ILRA -> $ILRA_SH"
else
    echo "CHECK: ILRA_SH=$ILRA_SH" >&2
fi

require_file "$REFERENCE_FASTA"
require_file "$REFERENCE_GFF"

if [[ ! -f "$PILON_JAR" ]]; then
    echo "CHECK: pilon.jar is expected in the current environment: $PILON_JAR" >&2
else
    echo "OK: Pilon -> $PILON_JAR"
fi
