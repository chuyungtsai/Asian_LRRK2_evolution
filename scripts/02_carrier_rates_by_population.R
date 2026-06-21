#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source_config <- function() {
  cfg <- Sys.getenv(c("RESULT_DIR"))
  if (cfg[["RESULT_DIR"]] == "") stop("Set RESULT_DIR or source scripts/00_config.sh before running.", call. = FALSE)
  cfg
}

args <- commandArgs(trailingOnly = TRUE)
carrier_file <- ifelse(length(args) >= 1, args[[1]], file.path(Sys.getenv("RESULT_DIR"), "lrrk2_carrier_status_long.tsv"))
out_dir <- ifelse(Sys.getenv("RESULT_DIR") != "", Sys.getenv("RESULT_DIR"), dirname(carrier_file))

if (!file.exists(carrier_file)) {
  stop("Missing carrier file: ", carrier_file, "\nExpected columns: variant, rsid, population, sample_id, alt_count", call. = FALSE)
}

x <- fread(carrier_file)
required <- c("variant", "rsid", "population", "sample_id", "alt_count")
if (!all(required %in% names(x))) stop("Carrier file missing required columns: ", paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)

summary <- x[, .(
  n = uniqueN(sample_id),
  carriers = sum(alt_count > 0, na.rm = TRUE),
  alt_count = sum(alt_count, na.rm = TRUE),
  carrier_rate = sum(alt_count > 0, na.rm = TRUE) / uniqueN(sample_id),
  allele_frequency = sum(alt_count, na.rm = TRUE) / (2 * uniqueN(sample_id))
), by = .(variant, rsid, population)]

setorder(summary, variant, population)
fwrite(summary, file.path(out_dir, "carrier_rate_by_population.tsv"), sep = "\t")

eas <- summary[population %in% c("CHB", "CHS", "JPT", "CDX", "KHV"),
  .(
    individuals = sum(n),
    carriers = sum(carriers),
    AC = sum(alt_count),
    AN = 2 * sum(n),
    AF = sum(alt_count) / (2 * sum(n))
  ),
  by = variant
]
fwrite(eas, file.path(out_dir, "carrier_rate_1000g_eas_aggregate.tsv"), sep = "\t")

print(summary)

