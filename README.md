# LRRK2 East Asian Haplotype Analysis Code

This folder contains publication-oriented scripts for the TWB WGS + 1000 Genomes East Asian analysis of **LRRK2 p.G2385R / rs34778348** and **p.R1628P / rs33949390**.

The scripts are organized by analysis module rather than by exploratory run order. They are intended to document and reproduce the main analyses underlying the manuscript figures, tables, and supplements.

## Data Requirements

The scripts assume access to:

- Phased TWB WGS VCF across the LRRK2 region, GRCh38.
- 1000 Genomes GRCh38 East Asian phased VCF or extracted LRRK2-region VCF.
- 1000 Genomes sample metadata with population labels: CHB, CHS, JPT, CDX, KHV.
- PLINK2 binary genotype files for merged TWB + 1000G EAS PCA.
- Optional selscan input files for EHH and XP-EHH.

Raw genotype data are not included in this code deposit.

## Software

- R >= 4.0
- Python >= 3.9
- PLINK2
- Beagle for phasing
- selscan for EHH / XP-EHH

R packages used by plotting and table scripts include:

- data.table
- dplyr
- ggplot2
- ape

## Script Overview

| Module | Script | Main output / manuscript use |
|---|---|---|
| Configuration | `scripts/00_config.sh` | Shared paths and variant constants |
| PCA/QC | `scripts/01_strict_pca_plink2.sh` | Figure 1C-D, Supplementary Table 2 |
| Carrier rates | `scripts/02_carrier_rates_by_population.R` | Table 2, Figure 1A-B, Supplementary Table 3 |
| Dense haplotypes | `scripts/03_dense_haplotype_reconstruction.py` | Figure 2C-D, Table 3, Supplementary Tables 5-6 |
| Cross-population exact haplotypes | `scripts/04_cross_population_exact_haplotype.R` | Figure 3, Supplementary Tables 7-8 |
| Haplotype clustering | `scripts/05_carrier_haplotype_clustering.R` | Figure 4, Supplementary Table 9 |
| EHH / XP-EHH | `scripts/06_selection_ehh_xpehh.sh` | Figure 5, Supplementary Table 11 |
| Manuscript tables | `scripts/07_build_manuscript_tables.R` | Tables 1-3 and supplementary table source files |

## Recommended Run Order

1. Edit `scripts/00_config.sh` to match the local data paths.
2. Run PCA/QC:

```bash
bash scripts/01_strict_pca_plink2.sh
```

3. Compute carrier rates:

```bash
Rscript scripts/02_carrier_rates_by_population.R
```

4. Reconstruct dense SNP haplotypes:

```bash
python3 scripts/03_dense_haplotype_reconstruction.py --variant G2385R
python3 scripts/03_dense_haplotype_reconstruction.py --variant R1628P
```

5. Compare exact haplotypes across TWB and 1000G EAS:

```bash
Rscript scripts/04_cross_population_exact_haplotype.R G2385R
Rscript scripts/04_cross_population_exact_haplotype.R R1628P
```

6. Cluster carrier haplotypes:

```bash
Rscript scripts/05_carrier_haplotype_clustering.R G2385R
Rscript scripts/05_carrier_haplotype_clustering.R R1628P
```

7. Run selection analyses:

```bash
bash scripts/06_selection_ehh_xpehh.sh
```

8. Build manuscript-ready table source files:

```bash
Rscript scripts/07_build_manuscript_tables.R
```

## Manuscript Mapping

- **Figure 1 / Table 2**: variant carrier rates and PCA overlays.
- **Figure 2 / Table 3**: dense SNP-defined haplotype separation and major carrier haplotypes.
- **Figure 3**: exact carrier haplotype geography across TWB and 1000G EAS.
- **Figure 4**: compact carrier haplotype trees and haplotype diversification.
- **Figure 5**: EHH and XP-EHH haplotype-based selection analyses.
- **Supplementary Tables 1-11**: variant coordinates, QC, global frequencies, PC GLM, dense haplotype windows, exact haplotype frequencies, clustering, age estimates, and XP-EHH summaries.

## Notes

- The dense haplotype reconstruction is custom code. It does not rely on a dedicated founder-haplotype package.
- Carrier chromosomes are defined by the phased ALT risk allele at the core variant.
- Exact haplotypes are represented as ordered dense SNP allele strings.
- 1000G matching is performed by chromosome, position, REF, and ALT allele.
- R1628P does not have a previously published historical founder interval analogous to the G2385R founder interval; R1628P windows are therefore empirical, centered on rs33949390.

