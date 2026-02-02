## ============================================================
##  BZea ANGSD QC Figures (Coverage + Missingness + AFS)
##  Inputs:
##    1) allchr.sample_qc.tsv  (sample x chr coverage QC)
##    2) chr*.afs_bins.tsv     (AFS bins per chr)
##  Outputs: PNGs in outdir
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

## ---- paths (EDIT if needed) ----
cov_file <- "/Users/nirwantandukar/Documents/Research/data/BZea/angsd_genotyping/coverage/allchr.sample_qc.tsv"
afs_dir  <- "/Users/nirwantandukar/Documents/Research/data/BZea/angsd_genotyping/afs"
outdir   <- "/Users/nirwantandukar/Documents/Research/data/BZea/angsd_genotyping/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## -------------------------------
## 1) Coverage / Missingness QC
## -------------------------------
cov <- fread(cov_file)
# enforce types
cov[, mean_dp := as.numeric(mean_dp)]
cov[, frac_missing := as.numeric(frac_missing)]
cov[, frac_dpge2 := as.numeric(frac_dpge2)]
cov[, n_sites := as.numeric(n_sites)]

# Keep chr1..chr10 order
cov[, chr := factor(chr, levels = paste0("chr", 1:10))]

# Genome-wide per-sample summary (weighted by n_sites)
cov_sample <- cov[, .(
  mean_dp_w       = weighted.mean(mean_dp, w = n_sites, na.rm = TRUE),
  frac_missing_w  = weighted.mean(frac_missing, w = n_sites, na.rm = TRUE),
  frac_dpge2_w    = weighted.mean(frac_dpge2, w = n_sites, na.rm = TRUE),
  total_sites     = sum(n_sites, na.rm = TRUE)
), by = sample]

# GG THEME
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 24,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 24,      # X‐axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 24,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),
    
    # ---- grids: light grey major grids for x + y, no minor grids ----
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 16),
    
    plot.margin    = margin(15, 15, 15, 15)
    
  )

# (A) Mean depth across samples (genome-wide)
mu  <- mean(cov_sample$mean_dp_w, na.rm = TRUE)
med <- median(cov_sample$mean_dp_w, na.rm = TRUE)

p_depth_hist <- ggplot(cov_sample, aes(x = mean_dp_w)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = mu,  linetype = "dotted", color = "red",  linewidth = 0.8) +
  geom_vline(xintercept = med, linetype = "dotted", color = "blue", linewidth = 0.8) +
  annotate("text", x = mu,  y = Inf, label = sprintf("mean = %.3f", mu),
           vjust = 1.2, hjust = -0.05, color = "red") +
  annotate("text", x = med, y = Inf, label = sprintf("median = %.3f", med),
           vjust = 2.6, hjust = -0.05, color = "blue") +
  labs(x = "Mean depth at teosinte SNP panel (weighted across chr1–chr10)",
       y = "Number of samples",
       title = "Coverage distribution across samples") +
  theme_bw() +
  plot_theme


quartz()
p_depth_hist


ggsave(file.path(outdir, "QC_meanDepth_hist.png"), p_depth_hist, width = 7.5, height = 4.5, dpi = 300)

# (B) Missingness across samples (genome-wide)
mu  <- mean(cov_sample$frac_missing_w, na.rm = TRUE)
med <- median(cov_sample$frac_missing_w, na.rm = TRUE)

p_miss_hist <- ggplot(cov_sample, aes(x = frac_missing_w)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = mu,  linetype = "dotted", color = "red",  linewidth = 0.8) +
  geom_vline(xintercept = med, linetype = "dotted", color = "blue", linewidth = 0.8) +
  annotate("text", x = mu,  y = Inf, label = sprintf("mean = %.1f%%", 100 * mu),
           vjust = 1.2, hjust = -0.05, color = "red") +
  annotate("text", x = med, y = Inf, label = sprintf("median = %.1f%%", 100 * med),
           vjust = 2.6, hjust = -0.05, color = "blue") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Fraction missing (DP=0) at teosinte SNP panel (weighted)",
       y = "Number of samples",
       title = "Missingness distribution across samples") +
  theme_bw() +
  plot_theme


quartz()
p_miss_hist
ggsave(file.path(outdir, "QC_missingness_hist.png"), p_miss_hist, width = 7.5, height = 4.5, dpi = 300)

# (C) Fraction of sites with DP>=2 across samples (genome-wide)
p_dpge2_hist <- ggplot(cov_sample, aes(x = frac_dpge2_w)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Fraction of teosinte SNPs with DP ≥ 2 (weighted)",
       y = "Number of samples",
       title = "Information content (DP≥2) across samples") +
  theme_bw() + plot_theme

quartz()
p_dpge2_hist
ggsave(file.path(outdir, "QC_dpge2_hist.png"), p_dpge2_hist, width = 7.5, height = 4.5, dpi = 300)

# (D) Per-chromosome distributions (boxplots)
p_chr_depth <- ggplot(cov, aes(x = chr, y = mean_dp)) +
  geom_boxplot(outlier.size = 0.4) +
  labs(x = "Chromosome", y = "Mean depth at teosinte SNP panel",
       title = "Depth by chromosome (per sample)") +
  theme_bw() + plot_theme
 
quartz()
p_chr_depth
ggsave(file.path(outdir, "QC_meanDepth_byChr_box.png"), p_chr_depth, width = 8, height = 4.5, dpi = 300)

p_chr_missing <- ggplot(cov, aes(x = chr, y = frac_missing)) +
  geom_boxplot(outlier.size = 0.4) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Chromosome", y = "Fraction missing (DP=0)",
       title = "Missingness by chromosome (per sample)") +
  theme_bw() + plot_theme

quartz()
p_chr_missing

ggsave(file.path(outdir, "QC_missingness_byChr_box.png"), p_chr_missing, width = 8, height = 4.5, dpi = 300)


## -------------------------------
## 2) AFS (Folded MAF spectrum)
## -------------------------------
afs_files <- list.files(afs_dir, pattern = "^chr[0-9]+\\.afs_bins\\.tsv$", full.names = TRUE)

afs_list <- lapply(afs_files, function(f) {
  dt <- fread(f)
  chr <- sub("\\.afs_bins\\.tsv$", "", basename(f))
  dt[, chr := factor(chr, levels = paste0("chr", 1:10))]
  dt
})
afs <- rbindlist(afs_list, use.names = TRUE, fill = TRUE)

afs[, bin_left := as.numeric(bin_left)]
afs[, count := as.numeric(count)]

# (E) AFS per chromosome (lines), log scale on y (counts are huge)
p_afs_lines <- ggplot(afs, aes(x = bin_left, y = count, group = chr, color = chr)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_y_log10(labels = label_number()) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  labs(x = "Folded MAF bin (left edge; bin width = 0.01)",
       y = "Count of SNPs (log10 scale)",
       title = "Folded allele frequency spectrum (per chromosome)") +
  theme_bw() +
  theme(legend.position = "right") + plot_theme

quartz()
p_afs_lines
ggsave(file.path(outdir, "AFS_byChr_lines_logY.png"), p_afs_lines, width = 9.0, height = 5.0, dpi = 300)

# (F) Genome-wide AFS (sum across chromosomes)
afs_genome <- afs[, .(count = sum(count, na.rm = TRUE)), by = bin_left]

p_afs_genome <- ggplot(afs_genome, aes(x = bin_left, y = count)) +
  geom_col() +
  scale_y_log10(labels = label_number()) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  labs(x = "Folded MAF bin (left edge; bin width = 0.01)",
       y = "Count of SNPs (log10 scale)",
       title = "Folded allele frequency spectrum (genome-wide; chr1–chr10 summed)") +
  theme_bw() + plot_theme

quartz()
p_afs_genome
ggsave(file.path(outdir, "AFS_genomewide_bar_logY.png"), p_afs_genome, width = 8.5, height = 4.8, dpi = 300)

## -------------------------------
## 3) Optional: simple QC table summary
## -------------------------------
qc_summary <- cov_sample[, .(
  n_samples = .N,
  mean_depth_mean = mean(mean_dp_w, na.rm = TRUE),
  mean_depth_median = median(mean_dp_w, na.rm = TRUE),
  miss_mean = mean(frac_missing_w, na.rm = TRUE),
  miss_median = median(frac_missing_w, na.rm = TRUE),
  dpge2_mean = mean(frac_dpge2_w, na.rm = TRUE),
  dpge2_median = median(frac_dpge2_w, na.rm = TRUE)
)]


fwrite(qc_summary, file.path(outdir, "QC_summary_stats.tsv"), sep = "\t")

message("Wrote figures to: ", outdir)












## Poisson expectation for missingness under low-pass sequencing
## Inputs you gave:
##   mean depth (lambda) = 0.822
##   observed missingness = 0.58  (DP=0 fraction)
##
## This script:
##   1) computes Poisson P(DP=0) = exp(-lambda)
##   2) compares expected vs observed
##   3) explains the deviation as mapping/baseQ filters + targeted SNP panel + ref bias
##
## Optional: if you have per-sample mean_dp_w and frac_missing_w in cov_sample,
## it also computes expected per sample and plots expected vs observed.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

lambda <- 0.822
obs_miss <- 0.58

exp_miss <- exp(-lambda)                 # P(DP=0)
exp_cov  <- 1 - exp_miss                 # P(DP>=1)

cat(sprintf("lambda (mean DP)          = %.3f\n", lambda))
cat(sprintf("Expected missing (DP=0)   = exp(-lambda) = %.3f (%.1f%%)\n", exp_miss, 100*exp_miss))
cat(sprintf("Expected non-missing      = 1-exp(-lambda) = %.3f (%.1f%%)\n", exp_cov, 100*exp_cov))
cat(sprintf("Observed missing (DP=0)   = %.3f (%.1f%%)\n", obs_miss, 100*obs_miss))
cat(sprintf("Observed - Expected       = %.3f (%.1f percentage points)\n",
            obs_miss - exp_miss, 100*(obs_miss - exp_miss)))

## Interpretation helper:
if (obs_miss < exp_miss) {
  cat("\nObserved missingness is LOWER than Poisson expectation.\n")
  cat("This can happen if DP is computed only at a targeted SNP panel and/or after filters,\n")
  cat("or if lambda reflects only callable sites (not genome-wide), making DP=0 less frequent.\n")
} else {
  cat("\nObserved missingness is HIGHER than Poisson expectation.\n")
  cat("This is common when you impose MAPQ/baseQ filters (e.g., minMapQ 30, minQ 20)\n")
  cat("and when aligning divergent haplotypes to a single reference (B73), which increases DP=0.\n")
}

## ---- OPTIONAL: if you have cov_sample with mean_dp_w and frac_missing_w ----
## Expectation per sample: exp(-mean_dp_w)
## Then compare observed vs expected.
if (exists("cov_sample")) {
  if (all(c("mean_dp_w", "frac_missing_w") %in% names(cov_sample))) {
    
    cov_sample <- as.data.table(cov_sample)
    cov_sample[, miss_exp := exp(-mean_dp_w)]
    cov_sample[, diff := frac_missing_w - miss_exp]
    
    p <- ggplot(cov_sample, aes(x = miss_exp, y = frac_missing_w)) +
      geom_point(alpha = 0.35, size = 1.2) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        x = "Expected missingness (Poisson: exp(-mean DP))",
        y = "Observed missingness (DP=0 fraction)",
        title = "Observed vs expected missingness under Poisson low-pass model",
        subtitle = "Dashed line: y = x"
      ) +
      theme_bw() +
      plot_theme
    quartz()
    print(p)
    
    cat("\nPer-sample summary (Observed - Expected):\n")
    print(cov_sample[, .(
      mean_diff = mean(diff, na.rm = TRUE),
      median_diff = median(diff, na.rm = TRUE),
      q05 = quantile(diff, 0.05, na.rm = TRUE),
      q95 = quantile(diff, 0.95, na.rm = TRUE)
    )])
  } else {
    message("cov_sample exists but does not have mean_dp_w and frac_missing_w columns.")
  }
}


# Given a mean sequencing depth of 0.82×, the Poisson expectation for zero coverage is ~44%. We observed a higher fraction of missing genotypes (~58%), consistent with the application of stringent mapping and base-quality filters and alignment of divergent teosinte haplotypes to the B73 reference genome. Across samples, observed missingness scaled linearly with Poisson expectations, indicating that missing data are driven by depth-dependent sampling rather than technical artifacts.