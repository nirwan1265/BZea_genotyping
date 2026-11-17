#!/usr/bin/env Rscript
## =============================================================================
## 32c — PCA of the pools' teosinte k-mer profiles (32b): does k-mer content alone separate the donor taxa?
##   Input: kmer_taxon_profile.tsv (one row per pool; f_<taxon> = share of that taxon's taxon-specific k-mer
##   hashes found in the pool). Each pool's f values are divided by their sum (so the depth and the amount of
##   donor DNA cancel out; only *which* teosinte the donor DNA resembles is left), log-ratio transformed
##   (centred log-ratio) and put into prcomp(). Points are coloured by the pedigree taxon. Using the taxon-level
##   sets (not one column per assembly) gives every taxon the same weight (Zh and Zl have one assembly each).
##   Also: a stacked bar per pool of the taxon-specific shares f_<taxon>, and the agreement table
##   (k-mer taxon vs pedigree taxon).
## Usage:  Rscript 32c_kmer_pca.R <kmer_taxon_profile.tsv> <out dir> [metadata_BZea.csv]
##   EXCLUDE (env, default Zm): pedigree groups left out (landraces carry no teosinte-specific k-mers).
##   With the metadata: per accession, the share of pools whose k-mer taxon differs from the pedigree taxon,
##   against the accession's elevation (kmer_mismatch_elevation.png, kmer_accession_agreement.tsv).
##   POOL_SUMMARY (env): step-21 pool_summary.tsv. Then the pool's teosinte k-mer content (f_teo, corrected for
##   depth: residual of log f_teo on log n_hashes within taxon) is compared with the HMM donor share
##   (kmer_vs_hmm.png, kmer_vs_hmm.tsv), pools flagged NOISY left out.
## Outputs: kmer_pca.png, kmer_taxon_shares.png, kmer_taxon_agreement.tsv, kmer_pca_scores.tsv
## Needs: R with data.table, ggplot2, patchwork.
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); if (length(a) < 2) stop("usage: Rscript 32c_kmer_pca.R <kmer_taxon_profile.tsv> <out dir>")
d <- fread(a[1]); out <- a[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
EXCL <- strsplit(Sys.getenv("EXCLUDE", "Zm"), ",")[[1]]
d <- d[role == "donor" & !label %in% EXCL]
cc <- grep("^c_", names(d), value = TRUE); ff <- grep("^f_Z", names(d), value = TRUE)
X <- as.matrix(d[, ..ff]); X <- X / rowSums(X)
keep <- is.finite(rowSums(X)) & rowSums(as.matrix(d[, ..ff])) > 0
eps <- min(X[X > 0], na.rm = TRUE) / 2
L <- log(pmax(X[keep, , drop = FALSE], eps)); L <- L - rowMeans(L)           # centred log-ratio
pc <- prcomp(L, center = TRUE, scale. = FALSE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
sc <- data.table(sample = d$sample[keep], taxon = d$label[keep], accession = d$accession[keep],
                 PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = if (ncol(pc$x) > 2) pc$x[, 3] else NA)
fwrite(sc, file.path(out, "kmer_pca_scores.tsv"), sep = "\t")
cols <- c(Zd = "#E69F00", Zh = "#56B4E9", Zl = "#009E73", Zv = "#CC79A7", Zx = "#0072B2")
labs_ <- c(Zd = "Zd diploperennis", Zh = "Zh huehuetenangensis", Zl = "Zl luxurians", Zv = "Zv parviglumis", Zx = "Zx mexicana")
sc[, PC4 := if (ncol(pc$x) > 3) pc$x[, 4] else NA]
pan <- function(a, b) {
  l <- data.table(asm = sub("^f_", "", ff), x = pc$rotation[, a], y = pc$rotation[, b])
  k <- 0.8 * max(abs(c(sc[[paste0("PC", a)]], sc[[paste0("PC", b)]]))) / max(abs(c(l$x, l$y)))
  ggplot(sc, aes(.data[[paste0("PC", a)]], .data[[paste0("PC", b)]], colour = taxon)) +
    geom_point(size = 1.4, alpha = 0.75) +
    geom_segment(data = l, aes(0, 0, xend = k * x, yend = k * y), inherit.aes = FALSE, colour = "grey40",
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_text(data = l, aes(k * x * 1.12, k * y * 1.12, label = asm), inherit.aes = FALSE, size = 3, colour = "grey20") +
    scale_colour_manual(values = cols, labels = labs_, name = "Pedigree taxon") +
    labs(x = sprintf("PC%d (%.1f%%)", a, ve[a]), y = sprintf("PC%d (%.1f%%)", b, ve[b])) + theme_classic(base_size = 11)
}
p1 <- (pan(1, 2) | pan(3, 4)) + plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("Teosinte k-mer profiles of %d donor pools (sourmash, k = 31; arrows = taxon-specific k-mer sets)", nrow(sc)))
ggsave(file.path(out, "kmer_pca.png"), p1, width = 13, height = 5.5, dpi = 200)

sh <- melt(d[, c("sample", "label", ff), with = FALSE], id.vars = c("sample", "label"), variable.name = "kt", value.name = "f")
sh[, kt := sub("^f_", "", kt)]; sh[, f := f / sum(f), by = sample]
sh[, sample := factor(sample, levels = d[order(label, sample)]$sample)]
p2 <- ggplot(sh, aes(sample, f, fill = kt)) + geom_col(width = 1) + facet_grid(~ label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cols, labels = labs_, name = "k-mer taxon") +
  labs(x = "Donor pools, grouped by pedigree taxon", y = "Share of taxon-specific k-mer hits") +
  theme_classic(base_size = 11) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(out, "kmer_taxon_shares.png"), p2, width = 11, height = 4, dpi = 200)

ag <- d[, .N, by = .(pedigree = label, kmer = assigned)][order(pedigree, -N)]
fwrite(dcast(ag, pedigree ~ kmer, value.var = "N", fill = 0), file.path(out, "kmer_taxon_agreement.tsv"), sep = "\t")
print(dcast(ag, pedigree ~ kmer, value.var = "N", fill = 0))
message(sprintf("agreement: %d / %d donor pools", sum(d$assigned == d$label), nrow(d)))
if (length(a) >= 3 && file.exists(a[3])) {
  m <- fread(a[3], encoding = "UTF-8"); setnames(m, 1, sub("^\ufeff", "", names(m)[1]))
  acc <- d[, .(pools = .N, agree = sum(assigned == label), other = { o <- assigned[assigned != label]; if (length(o)) names(which.max(table(o))) else "" }),
           by = .(taxon = label, accession)]
  acc <- merge(acc, m[, .(accession = `accession-id`, race, state, elevation)], by = "accession", all.x = TRUE)
  acc[, share_other := round(1 - agree / pools, 3)]
  fwrite(acc[order(taxon, elevation)], file.path(out, "kmer_accession_agreement.tsv"), sep = "\t")
  p3 <- ggplot(acc, aes(elevation, share_other, colour = taxon, size = pools)) + geom_point(alpha = 0.8) +
    scale_colour_manual(values = cols, labels = labs_, name = "Pedigree taxon") + scale_size_area(max_size = 5, name = "Pools") +
    labs(x = "Elevation of the accession (m)", y = "Share of pools with another k-mer taxon",
         title = "Per accession: pools whose k-mer taxon differs from the pedigree taxon") + theme_classic(base_size = 11)
  ggsave(file.path(out, "kmer_mismatch_elevation.png"), p3, width = 8, height = 4.5, dpi = 200)
  message("WROTE ", out, "/kmer_mismatch_elevation.png, kmer_accession_agreement.tsv")
}
message("WROTE ", out, "/kmer_pca.png, kmer_taxon_shares.png, kmer_taxon_agreement.tsv, kmer_pca_scores.tsv")
PS <- Sys.getenv("POOL_SUMMARY", "")
if (PS != "" && file.exists(PS)) {
  ps <- fread(PS)[, .(sample = SAMPLE, FLAG, DONOR_SHARE)]
  v <- merge(d[, .(sample, taxon = label, f_teo, n_hashes)], ps, by = "sample")[FLAG == "ok" & f_teo > 0]
  v[, f_teo_depth_corrected := resid(lm(log(f_teo) ~ log(n_hashes))), by = taxon]
  rr <- v[, .(pools = .N, spearman_rho = round(cor(f_teo_depth_corrected, DONOR_SHARE, method = "spearman"), 3)), by = taxon][order(taxon)]
  fwrite(rr, file.path(out, "kmer_vs_hmm.tsv"), sep = "\t"); print(rr)
  p4 <- ggplot(v, aes(100 * DONOR_SHARE, f_teo_depth_corrected, colour = taxon)) + geom_point(size = 1, alpha = 0.6) +
    facet_wrap(~ taxon, nrow = 1) + scale_colour_manual(values = cols, guide = "none") +
    geom_text(data = rr, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f", spearman_rho)), hjust = -0.1, vjust = 1.3,
              colour = "black", size = 3.2, inherit.aes = FALSE) +
    labs(x = "Donor alleles per pool, k-mer HMM (%)", y = "Teosinte k-mers in the reads\n(log, corrected for depth)",
         title = "Two measures of the donor share from the same reads: teosinte-specific k-mers vs the HMM") +
    theme_classic(base_size = 11)
  ggsave(file.path(out, "kmer_vs_hmm.png"), p4, width = 12, height = 3.6, dpi = 200)
  message("WROTE ", out, "/kmer_vs_hmm.png, kmer_vs_hmm.tsv")
}
