#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: 05_carrier_haplotype_clustering.R <G2385R|R1628P>", call. = FALSE)

variant <- args[[1]]
result_dir <- Sys.getenv("RESULT_DIR")
if (result_dir == "") stop("Set RESULT_DIR.", call. = FALSE)

hap_file <- file.path(result_dir, "carrier_haplotypes", paste0(variant, "_carrier_haplotype_matrix.tsv"))
if (!file.exists(hap_file)) {
  stop("Missing carrier haplotype matrix: ", hap_file, "\nRows should be carrier chromosomes; columns should be dense SNP allele codes.", call. = FALSE)
}

x <- fread(hap_file)
meta_cols <- intersect(c("haplotype_id", "population", "source"), names(x))
meta <- x[, ..meta_cols]
mat <- as.matrix(x[, setdiff(names(x), meta_cols), with = FALSE])
storage.mode(mat) <- "integer"
rownames(mat) <- meta$haplotype_id

dist_hamming <- dist(mat, method = "manhattan") / ncol(mat)
hc <- hclust(dist_hamming, method = "average")
phy <- as.phylo(hc)

hap_key <- apply(mat, 1, paste0, collapse = "")
unique_keys <- unique(hap_key)
unique_first_idx <- match(unique_keys, hap_key)
unique_mat <- mat[unique_first_idx, , drop = FALSE]
unique_ids <- paste0("H", seq_along(unique_keys))
rownames(unique_mat) <- unique_ids

unique_dist_hamming <- dist(unique_mat, method = "manhattan") / ncol(unique_mat)
unique_hc <- hclust(unique_dist_hamming, method = "average")
unique_phy <- as.phylo(unique_hc)

unique_summary <- data.table(
  unique_haplotype_id = unique_ids,
  hap_key = unique_keys,
  representative_haplotype_id = rownames(mat)[unique_first_idx]
)
unique_summary <- merge(
  unique_summary,
  data.table(hap_key = hap_key, haplotype_id = rownames(mat), population = meta$population),
  by = "hap_key",
  allow.cartesian = TRUE
)

unique_counts <- unique_summary[, .(
  carrier_chromosomes = .N,
  populations = paste(sort(unique(population)), collapse = ",")
), by = .(unique_haplotype_id, representative_haplotype_id)]

out_dir <- file.path(result_dir, "haplotype_clustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(unique_counts, file.path(out_dir, paste0(variant, "_unique_carrier_haplotype_summary.tsv")), sep = "\t")
saveRDS(hc, file.path(out_dir, paste0(variant, "_full_haplotype_hclust.rds")))
saveRDS(unique_hc, file.path(out_dir, paste0(variant, "_collapsed_unique_haplotype_hclust.rds")))

pdf(file.path(out_dir, paste0(variant, "_collapsed_unique_haplotype_tree.pdf")), width = 10, height = 6)
plot(unique_phy, type = "phylogram", direction = "rightwards", show.tip.label = FALSE)
axisPhylo()
dev.off()

cat("Wrote clustering outputs to ", out_dir, "\n", sep = "")

