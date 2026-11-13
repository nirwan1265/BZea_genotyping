#!/usr/bin/env Rscript
# =============================================================================
# 12_signal_check.R — are the teosinte-panel sites alone enough, per taxon? (base R)
#   Usage: Rscript 12_signal_check.R <D11 dir> <genome size in Mb> <output dir>
#
#   For each site class (panel / non-panel PASS sites):
#     background = ALT fraction of the B73 controls (split-half, from step 11)
#     excess     = pool ALT fraction - background
#     d          = excess / 0.125  (a BC2-derived line is ~12.5% donor, so this is about
#                  the extra ALT fraction inside a fully-donor region)
#   For a typical (median) pool of each taxon, per 1 Mb window:
#     reads/Mb, ALT reads expected in a B73 region (reads x background) and in a
#     fully-donor region (reads x (background + d)), and
#     z_fixed = (ALT donor - ALT B73) / sqrt(ALT B73 + 1)
#   z_fixed >= 10 on panel sites -> panel sites alone are enough for that taxon.
#   (A segregating window gives roughly half of that.)
#
#   Outputs: taxon_signal.tsv, pool_signal.tsv, site_filter_summary.tsv
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: Rscript 12_signal_check.R <D11> <genome_Mb> <OUTDIR>")
d11 <- args[1]; genome_mb <- as.numeric(args[2]); outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
donor_share <- 0.125

rd <- function(pattern) {
  f <- list.files(d11, pattern = pattern, full.names = TRUE)
  if (!length(f)) stop("no ", pattern, " in ", d11)
  do.call(rbind, lapply(f, read.delim, stringsAsFactors = FALSE))
}
pools <- read.delim(file.path(d11, "pools.tsv"), stringsAsFactors = FALSE)
reads <- rd("\\.pool_reads\\.tsv$")
bgt   <- rd("\\.ctrl_background\\.tsv$")
summ  <- rd("\\.summary\\.tsv$")
write.table(summ, file.path(outdir, "site_filter_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cols <- c("PANEL_REF", "PANEL_ALT", "NONPANEL_REF", "NONPANEL_ALT")
p <- aggregate(reads[, cols], by = list(SAMPLE = reads$SAMPLE), FUN = sum)
p <- merge(pools, p, by = "SAMPLE")

bg <- c(PANEL    = sum(bgt$PANEL_ALT)    / max(sum(bgt$PANEL_REF + bgt$PANEL_ALT), 1),
        NONPANEL = sum(bgt$NONPANEL_ALT) / max(sum(bgt$NONPANEL_REF + bgt$NONPANEL_ALT), 1))
for (cl in names(bg)) {
  p[[paste0(cl, "_READS")]]    <- p[[paste0(cl, "_REF")]] + p[[paste0(cl, "_ALT")]]
  p[[paste0(cl, "_ALT_FRAC")]] <- p[[paste0(cl, "_ALT")]] / pmax(p[[paste0(cl, "_READS")]], 1)
  p[[paste0(cl, "_EXCESS")]]   <- p[[paste0(cl, "_ALT_FRAC")]] - bg[[cl]]
}
p <- p[order(p$ROLE, p$GROUP, p$SAMPLE), ]
write.table(p, file.path(outdir, "pool_signal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

don <- p[p$ROLE == "donor", ]
tx <- do.call(rbind, lapply(split(don, don$GROUP), function(d) {
  r <- data.frame(TAXON = d$GROUP[1], N_POOLS = nrow(d))
  for (cl in names(bg)) {
    rpm  <- median(d[[paste0(cl, "_READS")]]) / genome_mb
    exc  <- median(d[[paste0(cl, "_EXCESS")]])
    b73  <- rpm * bg[[cl]]
    dnr  <- rpm * (bg[[cl]] + max(exc, 0) / donor_share)
    r[[paste0(cl, "_READS_PER_MB")]]     <- round(rpm, 1)
    r[[paste0(cl, "_MEDIAN_EXCESS")]]    <- signif(exc, 3)
    r[[paste0(cl, "_ALT_PER_MB_B73")]]   <- round(b73, 2)
    r[[paste0(cl, "_ALT_PER_MB_DONOR")]] <- round(dnr, 1)
    r[[paste0(cl, "_Z_FIXED")]]          <- round((dnr - b73) / sqrt(b73 + 1), 1)
  }
  r$PANEL_ENOUGH <- ifelse(r$PANEL_Z_FIXED >= 10, "yes", "NO - add non-panel sites")
  r
}))
write.table(tx, file.path(outdir, "taxon_signal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tot <- colSums(summ[, c("SITES", "CTRL_ALT", "PANEL_MISMATCH", "DEPTH_OR_THIRD_BASE", "PASS", "PASS_PANEL", "PASS_NONPANEL")])
cat("Sites, all chromosomes:\n"); print(tot)
cat(sprintf("\nB73 background ALT fraction: panel sites %.5f | non-panel sites %.5f\n", bg[["PANEL"]], bg[["NONPANEL"]]))
cat("Pools used:", paste(names(table(p$ROLE)), table(p$ROLE), sep = "=", collapse = " "), "\n\n")
options(width = 250)
print(tx[, c("TAXON", "N_POOLS", "PANEL_READS_PER_MB", "PANEL_ALT_PER_MB_B73", "PANEL_ALT_PER_MB_DONOR", "PANEL_Z_FIXED",
             "NONPANEL_ALT_PER_MB_B73", "NONPANEL_ALT_PER_MB_DONOR", "NONPANEL_Z_FIXED", "PANEL_ENOUGH")], row.names = FALSE)
cat("\nWrote:", file.path(outdir, c("taxon_signal.tsv", "pool_signal.tsv", "site_filter_summary.tsv")), sep = "\n  ")
