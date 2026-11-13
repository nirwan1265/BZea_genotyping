#!/usr/bin/env Rscript
# =============================================================================
# 15_collect_hmm.R — gather the per-pool HMM results (base R)
#   Usage: Rscript 15_collect_hmm.R <D14> <pools.tsv> <outdir>
#   Outputs: pool_summary.tsv (one row per pool), tracts_all.tsv.gz (all pools' tracts)
#   Prints the median donor share per taxon: ~0.125 expected for BC2-derived lines, ~0 for B73.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: Rscript 15_collect_hmm.R <D14> <pools.tsv> <outdir>")
d14 <- args[1]; pools <- read.delim(args[2], stringsAsFactors = FALSE); outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
sf <- list.files(d14, pattern = "\\.summary\\.tsv$", full.names = TRUE)
if (!length(sf)) stop("no *.summary.tsv in ", d14)
s <- do.call(rbind, lapply(sf, read.delim, stringsAsFactors = FALSE))
s <- merge(pools[, c("SAMPLE", "ROLE", "GROUP", "ACCESSION")], s, by = "SAMPLE")
write.table(s, file.path(outdir, "pool_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tr <- do.call(rbind, lapply(s$SAMPLE, function(x) {
  f <- file.path(d14, paste0(x, ".tracts.bed.gz"))
  if (!file.exists(f)) return(NULL)
  t <- read.delim(gzfile(f), stringsAsFactors = FALSE); if (!nrow(t)) return(NULL)
  cbind(SAMPLE = x, t)
}))
gz <- gzfile(file.path(outdir, "tracts_all.tsv.gz"), "w"); write.table(tr, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

cat("pools with results:", nrow(s), "of", sum(pools$ROLE %in% c("control", "donor")), "\n\n")
print(do.call(rbind, lapply(split(s, s$GROUP), function(d) data.frame(
  GROUP = d$GROUP[1], N = nrow(d),
  MEDIAN_DONOR_SHARE = median(d$DONOR_SHARE, na.rm = TRUE),
  MEDIAN_FRAC_DONOR = median(d$FRAC_DONOR, na.rm = TRUE),
  MEDIAN_FRAC_SEG = median(d$FRAC_SEG, na.rm = TRUE),
  MEDIAN_TRACTS = median(d$N_TRACTS)))), row.names = FALSE)
cat("\nWrote:", file.path(outdir, c("pool_summary.tsv", "tracts_all.tsv.gz")), sep = "\n  ")
