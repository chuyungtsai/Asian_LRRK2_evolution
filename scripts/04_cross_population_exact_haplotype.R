#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: 04_cross_population_exact_haplotype.R <G2385R|R1628P>", call. = FALSE)

variant <- args[[1]]
result_dir <- Sys.getenv("RESULT_DIR")
kg_vcf <- Sys.getenv("KG_LRRK2_PHASED_VCF")
sample_meta_file <- Sys.getenv("KG_SAMPLE_METADATA")

if (result_dir == "" || kg_vcf == "" || sample_meta_file == "") {
  stop("Set RESULT_DIR, KG_LRRK2_PHASED_VCF, and KG_SAMPLE_METADATA.", call. = FALSE)
}

target_pops <- c("TW", "CHB", "CHS", "JPT", "CDX", "KHV")

read_vcf_selected <- function(path, keep_pos) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  keep_pos <- as.integer(keep_pos)
  samples <- NULL
  rows <- list()
  repeat {
    line <- readLines(con, n = 1)
    if (!length(line)) break
    if (startsWith(line, "##")) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (startsWith(line, "#CHROM")) {
      samples <- fields[10:length(fields)]
      next
    }
    pos <- as.integer(fields[2])
    if (!(pos %in% keep_pos)) next
    gt <- sub(":.*$", "", fields[10:length(fields)])
    alleles <- as.integer(unlist(strsplit(gt, "[|/]")))
    rows[[length(rows) + 1L]] <- list(
      chrom = fields[1], pos = pos, id = fields[3], ref = fields[4], alt = fields[5], alleles = alleles
    )
  }
  variants <- rbindlist(lapply(rows, function(z) data.table(chrom = z$chrom, pos = z$pos, id = z$id, ref = z$ref, alt = z$alt)))
  mat <- do.call(rbind, lapply(rows, `[[`, "alleles"))
  list(samples = samples, variants = variants, hap_by_variant = mat)
}

if (variant == "G2385R") {
  hap_file <- file.path(result_dir, "dense_haplotype", "G2385R_major_haplotype_snps.tsv")
  core_pos <- 40363526L
} else if (variant == "R1628P") {
  hap_file <- file.path(result_dir, "dense_haplotype", "R1628P_major_haplotype_snps.tsv")
  core_pos <- 40320043L
} else {
  stop("Unknown variant: ", variant, call. = FALSE)
}

if (!file.exists(hap_file)) {
  stop("Missing major haplotype SNP file: ", hap_file, call. = FALSE)
}

founder <- fread(hap_file)
founder[, chrom_1kg := sub("^chr", "", chrom)]
founder[, key := paste(chrom_1kg, pos, ref, alt, sep = ":")]

kg <- read_vcf_selected(kg_vcf, founder$pos)
kg_vars <- copy(kg$variants)
kg_vars[, chrom_1kg := sub("^chr", "", chrom)]
kg_vars[, key := paste(chrom_1kg, pos, ref, alt, sep = ":")]
idx <- match(founder$key, kg_vars$key)
matched <- !is.na(idx)
if (!any(matched)) stop("No dense haplotype SNPs matched 1000G by chrom:pos:ref:alt.", call. = FALSE)

kg_mat <- kg$hap_by_variant[idx[matched], , drop = FALSE]
founder_codes <- founder$shared_allele_code[matched]
match_by_hap <- colSums(kg_mat == founder_codes) == length(founder_codes)

meta <- fread(sample_meta_file, fill = TRUE)
meta <- meta[, 1:4]
setnames(meta, c("sample", "pop", "super_pop", "gender"))
hap_sample <- rep(kg$samples, each = 2)
hap_pop <- meta$pop[match(hap_sample, meta$sample)]

core_row <- match(core_pos, kg_vars$pos[idx[matched]])
core_alt_by_hap <- if (!is.na(core_row)) kg_mat[core_row, ] == 1L else rep(NA, ncol(kg_mat))

rows <- rbindlist(lapply(setdiff(target_pops, "TW"), function(pop) {
  keep <- hap_pop == pop
  data.table(
    population = pop,
    matched_haplotype_snps = sum(matched),
    total_chromosomes = sum(keep, na.rm = TRUE),
    alt_chromosomes = sum(core_alt_by_hap[keep], na.rm = TRUE),
    exact_haplotype_count = sum(match_by_hap[keep], na.rm = TRUE),
    exact_haplotype_frequency = sum(match_by_hap[keep], na.rm = TRUE) / sum(keep, na.rm = TRUE),
    exact_haplotype_among_carriers = ifelse(
      sum(core_alt_by_hap[keep], na.rm = TRUE) > 0,
      sum(match_by_hap[keep], na.rm = TRUE) / sum(core_alt_by_hap[keep], na.rm = TRUE),
      NA_real_
    )
  )
}))

out_dir <- file.path(result_dir, "cross_population_exact_haplotype")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(rows, file.path(out_dir, paste0(variant, "_exact_haplotype_frequency_by_population.tsv")), sep = "\t")
fwrite(data.table(founder, matched_1000g = matched), file.path(out_dir, paste0(variant, "_haplotype_snps_matched_to_1000g.tsv")), sep = "\t")
print(rows)

