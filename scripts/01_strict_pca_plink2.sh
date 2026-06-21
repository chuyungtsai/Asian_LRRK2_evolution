#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

OUT="${RESULT_DIR}/strict_common_pca"
mkdir -p "${OUT}"

# Strict common-SNP PCA set:
# - common SNPs only
# - high call rate
# - low sample missingness
# - HWE consistency
# - LD pruning
# - long-range LD exclusion
# - chromosome 12 exclusion to avoid local LRRK2 signal driving PCA

"${PLINK2}" \
  --bfile "${MERGED_PLINK_PREFIX}" \
  --maf 0.05 \
  --geno 0.005 \
  --mind 0.02 \
  --hwe 1e-6 \
  --snps-only just-acgt \
  --not-chr 12 \
  --exclude range "${LONG_RANGE_LD_FILE}" \
  --indep-pairwise 200 50 0.1 \
  --out "${OUT}/merged_strict"

"${PLINK2}" \
  --bfile "${MERGED_PLINK_PREFIX}" \
  --extract "${OUT}/merged_strict.prune.in" \
  --pca 20 approx \
  --out "${OUT}/merged_strict_pca"

echo "PCA output written to ${OUT}"

