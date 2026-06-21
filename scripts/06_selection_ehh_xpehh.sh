#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

OUT="${RESULT_DIR}/selection"
mkdir -p "${OUT}"

# This script documents the selscan workflow used for local EHH and XP-EHH.
# It assumes phased VCFs have already been converted to selscan .hap/.map files.

TWB_HAP="${TWB_HAP:-${OUT}/TWB_LRRK2.hap}"
TWB_MAP="${TWB_MAP:-${OUT}/TWB_LRRK2.map}"

"${SELSCAN}" --ehh "${G2385R_ID}" --hap "${TWB_HAP}" --map "${TWB_MAP}" --out "${OUT}/${G2385R_ID}.ehh"
"${SELSCAN}" --ehh "${R1628P_ID}" --hap "${TWB_HAP}" --map "${TWB_MAP}" --out "${OUT}/${R1628P_ID}.ehh"

# XP-EHH example. Repeat with each 1000G EAS population reference panel.
# Required environment variables for each comparison:
#   REF_HAP, REF_MAP, COMPARISON_LABEL
if [[ -n "${REF_HAP:-}" && -n "${REF_MAP:-}" && -n "${COMPARISON_LABEL:-}" ]]; then
  "${SELSCAN}" \
    --xpehh \
    --hap "${TWB_HAP}" \
    --map "${TWB_MAP}" \
    --ref "${REF_HAP}" \
    --out "${OUT}/xpehh_TWB_vs_${COMPARISON_LABEL}"

  "${SELSCAN}" norm \
    --xpehh \
    --files "${OUT}/xpehh_TWB_vs_${COMPARISON_LABEL}.xpehh.out" \
    --out "${OUT}/xpehh_TWB_vs_${COMPARISON_LABEL}.norm"
fi

echo "Selection outputs written to ${OUT}"

