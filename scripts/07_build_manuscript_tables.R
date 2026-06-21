#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

result_dir <- Sys.getenv("RESULT_DIR")
if (result_dir == "") stop("Set RESULT_DIR.", call. = FALSE)

out_dir <- file.path(result_dir, "manuscript_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

carrier_file <- file.path(result_dir, "carrier_rate_by_population.tsv")
if (file.exists(carrier_file)) {
  carrier <- fread(carrier_file)
  table2 <- dcast(
    carrier,
    population ~ variant,
    value.var = c("carriers", "n", "carrier_rate", "allele_frequency")
  )
  fwrite(table2, file.path(out_dir, "Table2_carrier_frequencies_source.tsv"), sep = "\t")
}

g_file <- file.path(result_dir, "dense_haplotype", "G2385R_dense_haplotype_window_summary.tsv")
r_file <- file.path(result_dir, "dense_haplotype", "R1628P_dense_haplotype_window_summary.tsv")
if (file.exists(g_file) && file.exists(r_file)) {
  g <- fread(g_file)[window == "historical_239kb_G2385R_to_D12S2517_offset"]
  r <- fread(r_file)[window == "239kb_centered"]
  table3 <- rbind(
    data.table(
      variant = "G2385R",
      interval = "Historical 239 kb analogue",
      dense_snps = g$variant_count,
      carrier_haplotype = paste0(g$exact_major_carrier_hap_count, "/", g$carrier_chromosomes),
      noncarrier_haplotype = paste0(g$exact_same_hap_noncarrier_count, "/", g$noncarrier_chromosomes),
      shared_block_bp = g$shared_size_bp
    ),
    data.table(
      variant = "R1628P",
      interval = "Centered 239 kb window",
      dense_snps = r$variant_count,
      carrier_haplotype = paste0(r$exact_major_carrier_hap_count, "/", r$carrier_chromosomes),
      noncarrier_haplotype = paste0(r$exact_same_hap_noncarrier_count, "/", r$noncarrier_chromosomes),
      shared_block_bp = r$shared_size_bp
    )
  )
  fwrite(table3, file.path(out_dir, "Table3_dense_major_carrier_haplotypes_source.tsv"), sep = "\t")
}

cat("Manuscript table source files written to ", out_dir, "\n", sep = "")

