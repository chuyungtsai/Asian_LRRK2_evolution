#!/usr/bin/env bash
set -euo pipefail

# Edit these paths before running the workflow.
PROJECT_DIR="${PROJECT_DIR:-/path/to/LRRK2_analysis}"
RESULT_DIR="${RESULT_DIR:-${PROJECT_DIR}/results/publication_lrrk2}"

TWB_LRRK2_PHASED_VCF="${TWB_LRRK2_PHASED_VCF:-${PROJECT_DIR}/data/TWB_LRRK2_chr12.phased.vcf.gz}"
KG_LRRK2_PHASED_VCF="${KG_LRRK2_PHASED_VCF:-${PROJECT_DIR}/data/1000G_EAS_LRRK2.phased.vcf.gz}"
KG_SAMPLE_METADATA="${KG_SAMPLE_METADATA:-${PROJECT_DIR}/data/1000G_samples.tsv}"

MERGED_PLINK_PREFIX="${MERGED_PLINK_PREFIX:-${PROJECT_DIR}/data/merged_TWB_1000G}"
LONG_RANGE_LD_FILE="${LONG_RANGE_LD_FILE:-${PROJECT_DIR}/data/long_range_ld_regions.txt}"

PLINK2="${PLINK2:-plink2}"
BEAGLE_JAR="${BEAGLE_JAR:-beagle.jar}"
SELSCAN="${SELSCAN:-selscan}"

CHR="12"
G2385R_ID="rs34778348"
G2385R_POS="40363526"
G2385R_REF="G"
G2385R_ALT="A"

R1628P_ID="rs33949390"
R1628P_POS="40320043"
R1628P_REF="G"
R1628P_ALT="C"

mkdir -p "${RESULT_DIR}"

