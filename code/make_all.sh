#!/usr/bin/env bash
# Regenerate every figure and in-text number for the manuscript.
# Usage: bash make_all.sh   (set V2X_DATA_DIR first; see README.md)
set -euo pipefail
cd "$(dirname "$0")"

python3 -c "import config; config.ensure_dirs(); config.check_data()"

for s in \
    dataset_stats.py \
    element1_radio_throughput.py \
    element2_dependence.py \
    element3_redundancy.py \
    element4_berlin.py
do
    if [ -f "$s" ]; then
        echo "=== $s ==="
        python3 "$s"
    else
        echo "=== $s (not present) — skipping ==="
    fi
done

echo "Done. Figures -> ../figures/ ; numbers -> numbers/"
