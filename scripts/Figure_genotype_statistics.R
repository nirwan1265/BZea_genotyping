library(tidyverse)

# -------------------------
# Helper: parse a bcftools stats section using its own header line
# (works even if column layout differs by bcftools version)
# -------------------------
read_bcftools_section <- function(stats_file, section = "PSC") {
  lines <- readLines(stats_file, warn = FALSE)
  
  # --- header: must be tab-delimited ("# PSC\t...")
  hdr_idx <- which(grepl(paste0("^#\\s*", section, "\\t"), lines))
  if (length(hdr_idx) == 0) {
    stop("No tab-delimited header found for section ", section,
         ". Quick check:\n",
         "  grep -n '^# ", section, "' ", stats_file, " | head\n",
         "  grep -n '^", section, "\\t' ", stats_file, " | head")
  }
  hdr <- lines[hdr_idx[1]]
  
  # strip leading "# " then split
  hdr_fields <- strsplit(sub("^#\\s*", "", hdr), "\t")[[1]]
  
  # drop the first field (the section label itself)
  hdr_fields <- hdr_fields[-1]
  
  # remove "[n]" prefixes, trim, and guarantee non-empty unique names
  coln <- gsub("^\\[[0-9]+\\]", "", hdr_fields)
  coln <- trimws(coln)
  bad  <- is.na(coln) | coln == ""
  if (any(bad)) coln[bad] <- paste0("V", which(bad))
  coln <- make.names(coln, unique = TRUE)
  
  # --- data lines
  dat_lines <- lines[grepl(paste0("^", section, "\\t"), lines)]
  if (length(dat_lines) == 0) stop("No data lines found for section ", section)
  
  df <- readr::read_tsv(I(dat_lines), col_names = FALSE, show_col_types = FALSE)
  
  # first column is the section label; drop it
  df <- df[, -1, drop = FALSE]
  
  # if columns don't match header, pad names
  if (ncol(df) > length(coln)) {
    coln <- c(coln, paste0("extra_", seq_len(ncol(df) - length(coln))))
  }
  names(df) <- coln[seq_len(ncol(df))]
  df
}

# -------------------------
# Load sample-level stats (PSC) from unimputed + imputed
# -------------------------
psc_unimp <- read_bcftools_section("data/qc/unimputed.bcftools.stats.txt", "PSC") %>%
  mutate(dataset = "unimputed")

psc_imp <- read_bcftools_section("data/qc/imputed.bcftools.stats.txt", "PSC") %>%
  mutate(dataset = "imputed")

psc <- bind_rows(psc_unimp, psc_imp)

# Inspect available columns (bcftools version dependent)
print(names(psc))

# Try to compute common QC metrics if the columns exist
# Typical PSC columns often include: sample, nHets, nNonRefHom, nRefHom, average_depth, nMissing
get_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

col_sample  <- get_col(psc, c("sample", "Sample", "INDV"))
col_hets    <- get_col(psc, c("nHets", "nHET", "nHets.1"))
col_altHom  <- get_col(psc, c("nNonRefHom", "nAltHom", "nNonRefHoms"))
col_refHom  <- get_col(psc, c("nRefHom", "nRefHoms"))
col_depth   <- get_col(psc, c("average_depth", "avg_depth", "average.depth", "mean_depth", "dp"))
col_missing <- get_col(psc, c("nMissing", "nMiss", "missing"))

if (is.null(col_sample)) stop("Couldn't find a sample column in PSC. Columns are: ", paste(names(psc), collapse=", "))

qc_samp <- psc %>%
  rename(sample = all_of(col_sample)) %>%
  mutate(
    nHets       = if (!is.null(col_hets)) as.numeric(.data[[col_hets]]) else NA_real_,
    nNonRefHom  = if (!is.null(col_altHom)) as.numeric(.data[[col_altHom]]) else NA_real_,
    nRefHom     = if (!is.null(col_refHom)) as.numeric(.data[[col_refHom]]) else NA_real_,
    meanDP      = if (!is.null(col_depth)) as.numeric(.data[[col_depth]]) else NA_real_,
    nMissing    = if (!is.null(col_missing)) as.numeric(.data[[col_missing]]) else NA_real_
  ) %>%
  mutate(
    # alt allele burden proxy (if counts exist)
    alt_alleles = ifelse(!is.na(nNonRefHom) & !is.na(nHets),
                         2*nNonRefHom + nHets, NA_real_),
    # heterozygosity rate among called genotypes (if counts exist)
    called = ifelse(!is.na(nRefHom) & !is.na(nNonRefHom) & !is.na(nHets),
                    nRefHom + nNonRefHom + nHets, NA_real_),
    het_rate = ifelse(!is.na(called) & called > 0 & !is.na(nHets),
                      nHets / called, NA_real_)
  )



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


# =========================
# INDIVIDUAL
# =========================

# Add derived QC metrics once (so plots work for both datasets)
qc_samp <- qc_samp %>%
  mutate(
    total_sites  = nRefHom + nNonRefHom + nHets + nMissing,
    missing_rate = ifelse(!is.na(total_sites) & total_sites > 0, nMissing / total_sites, NA_real_),
    
    total_alleles   = 2 * (nRefHom + nNonRefHom + nHets),
    alt_allele_freq = ifelse(!is.na(total_alleles) & total_alleles > 0 & !is.na(alt_alleles),
                             alt_alleles / total_alleles, NA_real_)
  )

# ---- plot function: makes SEPARATE plots per dataset and saves them ----
make_qc_plots_one <- function(df, ds, out_dir = "data/qc/figs_by_dataset",
                              prefix = NULL, open_quartz = FALSE) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(prefix)) prefix <- ds
  
  d <- df %>% filter(dataset == ds)
  
  # 1) Depth distribution
  p_depth <- d %>%
    filter(!is.na(meanDP)) %>%
    ggplot(aes(x = meanDP)) +
    geom_histogram(bins = 60, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Mean DP (across variants)",
      y = "Samples",
      title = paste0("Per-sample depth: ", ds)
    ) + plot_theme
  
  # 2) Residual heterozygosity
  p_het <- d %>%
    filter(!is.na(het_rate)) %>%
    ggplot(aes(x = het_rate)) +
    geom_histogram(bins = 60, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Residual heterozygosity (nHets / called)",
      y = "Samples",
      title = paste0("Residual heterozygosity: ", ds)
    ) + plot_theme
  
  # 3) Hets vs alt burden
  p_het_alt <- d %>%
    filter(!is.na(het_rate), !is.na(alt_alleles)) %>%
    ggplot(aes(x = alt_alleles, y = het_rate)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Alt allele burden (2*AltHom + Het)",
      y = "Residual heterozygosity",
      title = paste0("Hets vs alt burden: ", ds)
    ) + plot_theme
  
  # 4) Missing rate distribution
  p_missing <- d %>%
    filter(!is.na(missing_rate)) %>%
    ggplot(aes(x = missing_rate)) +
    geom_histogram(bins = 60, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Missing genotype rate",
      y = "Samples",
      title = paste0("Missingness: ", ds)
    ) +
    scale_x_continuous(labels = scales::percent) + plot_theme
  
  # 5) Per-sample alt allele frequency (your “alts” summary)
  p_afs <- d %>%
    filter(!is.na(alt_allele_freq)) %>%
    ggplot(aes(x = alt_allele_freq)) +
    geom_density(alpha = 0.7) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Per-sample alternate allele frequency",
      y = "Density",
      title = paste0("Per-sample alt allele frequency: ", ds)
    ) +
    scale_x_continuous(labels = scales::percent) + plot_theme
  
  # 6) Ts/Tv (only if PSC has those columns; otherwise skip quietly)
  if (all(c("nTransitions", "nTransversions") %in% names(d))) {
    p_tstv <- d %>%
      mutate(
        ts_tv_ratio = nTransitions / nTransversions,
        ts_tv_ratio = ifelse(is.infinite(ts_tv_ratio), NA_real_, ts_tv_ratio)
      ) %>%
      filter(!is.na(ts_tv_ratio)) %>%
      ggplot(aes(x = ts_tv_ratio)) +
      geom_histogram(bins = 60, alpha = 0.8) +
      theme_minimal(base_size = 14) +
      labs(
        x = "Ts/Tv ratio",
        y = "Samples",
        title = paste0("Ts/Tv: ", ds) 
      ) + plot_theme
  } else {
    p_tstv <- NULL
  }
  
  # ---- Save PDFs ----
  ggsave(file.path(out_dir, paste0(prefix, "_depth.pdf")),    p_depth,   width = 12, height = 6)
  ggsave(file.path(out_dir, paste0(prefix, "_het.pdf")),      p_het,     width = 12, height = 6)
  ggsave(file.path(out_dir, paste0(prefix, "_het_vs_alt.pdf")), p_het_alt, width = 12, height = 6)
  ggsave(file.path(out_dir, paste0(prefix, "_missing.pdf")),  p_missing, width = 12, height = 6)
  ggsave(file.path(out_dir, paste0(prefix, "_alt_af.pdf")),   p_afs,     width = 12, height = 6)
  
  if (!is.null(p_tstv)) {
    ggsave(file.path(out_dir, paste0(prefix, "_tstv.pdf")), p_tstv, width = 12, height = 6)
  }
  
  # optional: pop up plots
  if (open_quartz) {
    quartz(); print(p_depth)
    quartz(); print(p_het)
    quartz(); print(p_het_alt)
    quartz(); print(p_missing)
    quartz(); print(p_afs)
    if (!is.null(p_tstv)) { quartz(); print(p_tstv) }
  }
  
  # return plots (handy in RStudio)
  invisible(list(
    depth = p_depth,
    het = p_het,
    het_vs_alt = p_het_alt,
    missing = p_missing,
    alt_af = p_afs,
    tstv = p_tstv
  ))
}

# ---- run for each dataset (SEPARATE outputs) ----
make_qc_plots_one(qc_samp, "unimputed", out_dir = "data/qc/figs_by_dataset", open_quartz = TRUE)
make_qc_plots_one(qc_samp, "imputed",   out_dir = "data/qc/figs_by_dataset", open_quartz = TRUE)

# ---- optional: separate summary table per dataset ----
summary_stats <- qc_samp %>%
  group_by(dataset) %>%
  summarise(
    n_samples     = n_distinct(sample),
    mean_depth    = mean(meanDP, na.rm = TRUE),
    sd_depth      = sd(meanDP, na.rm = TRUE),
    mean_missing  = mean(missing_rate, na.rm = TRUE),
    mean_het      = mean(het_rate, na.rm = TRUE),
    mean_alt_af   = mean(alt_allele_freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

print(summary_stats)































# -------------------------
# -------------------------
#### IMPUTED AND UNIMPUTED SIDE BY SIDE
# -------------------------
# -------------------------

# -------------------------
# Plot 1: Mean depth per sample (unimputed vs imputed)
# -------------------------
p_depth <- qc_samp %>%
  filter(!is.na(meanDP)) %>%
  ggplot(aes(x = dataset, y = meanDP)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5) +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Mean DP (across variants)", title = "Per-sample depth (VCF DP)")

quartz()
p_depth

# -------------------------
# Plot 2: Residual heterozygosity distribution
# -------------------------
p_het <- qc_samp %>%
  filter(!is.na(het_rate)) %>%
  ggplot(aes(x = het_rate, fill = dataset)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  theme_minimal(base_size = 12) +
  labs(x = "Residual heterozygosity (nHets / called)", y = "Samples",
       title = "Residual heterozygosity (BC2S3 sanity check)")

quartz()
p_het

# -------------------------
# Plot 3: Het vs alt allele burden (your “hets vs alts” plot)
# -------------------------
p_het_alt <- qc_samp %>%
  filter(!is.na(het_rate), !is.na(alt_alleles)) %>%
  ggplot(aes(x = alt_alleles, y = het_rate, color = dataset)) +
  geom_point(alpha = 0.7) +
  theme_minimal(base_size = 12) +
  labs(x = "Alt allele burden (2*AltHom + Het)", y = "Residual heterozygosity",
       title = "Heterozygosity vs alternate allele burden")



quartz()
p_het_alt



# -------------------------
# Plot 4: Missing Data Rates Comparison
# -------------------------

# Calculate missing rate
qc_samp <- qc_samp %>%
  mutate(
    total_sites = nRefHom + nNonRefHom + nHets + nMissing,
    missing_rate = nMissing / total_sites
  )

p_missing <- ggplot(qc_samp, aes(x = dataset, y = missing_rate)) +
  geom_violin(fill = "lightblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Missing Genotype Rate", 
       title = "Reduction in Missing Data After Imputation",
       subtitle = "1,500 maize B73×teosinte lines (0.8x coverage)") +
  scale_y_continuous(labels = scales::percent)

quartz()
print(p_missing)

# -------------------------
# Plot 5: Genotype Concordance/Consistency Plot
# -------------------------

# Since you have both datasets for same samples
qc_combined <- qc_samp %>%
  select(sample, dataset, het_rate, missing_rate, meanDP) %>%
  pivot_wider(names_from = dataset, values_from = c(het_rate, missing_rate, meanDP))

p_concordance <- ggplot(qc_combined, 
                        aes(x = het_rate_unimputed, y = het_rate_imputed)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(x = "Heterozygosity (Unimputed)", 
       y = "Heterozygosity (Imputed)",
       title = "Genotype Concordance After Imputation",
       subtitle = "Each point = one maize line") +
  coord_equal()

quartz()
print(p_concordance)


# -------------------------
# Plot 6: Allele Frequency Spectrum Comparison
# -------------------------
# Calculate per-sample alternate allele frequency
qc_samp <- qc_samp %>%
  mutate(
    total_alleles = 2 * (nRefHom + nNonRefHom + nHets),
    alt_allele_freq = alt_alleles / total_alleles
  )

p_afs <- ggplot(qc_samp, aes(x = alt_allele_freq, fill = dataset)) +
  geom_density(alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(x = "Alternate Allele Frequency (per sample)",
       y = "Density",
       title = "Allele Frequency Spectrum",
       subtitle = "Comparison of unimputed vs imputed datasets") +
  scale_x_continuous(labels = scales::percent)

quartz()
print(p_afs)



# -------------------------
# Plot 7: Transition/Transversion Ratio (Ts/Tv)
# -------------------------
# Calculate per-sample Ts/Tv ratio
qc_samp <- qc_samp %>%
  mutate(
    ts_tv_ratio = nTransitions / nTransversions,
    ts_tv_ratio = ifelse(is.infinite(ts_tv_ratio), NA, ts_tv_ratio)
  )

p_tstv <- ggplot(qc_samp %>% filter(!is.na(ts_tv_ratio)), 
                 aes(x = dataset, y = ts_tv_ratio)) +
  geom_violin(fill = "lightgreen", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Ts/Tv Ratio",
       title = "Transition/Transversion Ratio",
       subtitle = "Expected ~2.0-2.1 for maize") +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "red")

quartz()
print(p_tstv)


library(kableExtra)

summary_stats <- qc_samp %>%
  group_by(dataset) %>%
  summarise(
    n_samples = n(),
    mean_depth = mean(meanDP, na.rm = TRUE),
    sd_depth = sd(meanDP, na.rm = TRUE),
    mean_missing = mean(missing_rate, na.rm = TRUE),
    mean_het = mean(het_rate, na.rm = TRUE),
    mean_ts_tv = mean(ts_tv_ratio, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

kable(summary_stats, caption = "Summary Statistics Before and After Imputation") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))


# Calculate imputation success rate
imputation_gain <- qc_samp %>%
  group_by(dataset) %>%
  summarise(
    mean_called_sites = mean(called, na.rm = TRUE),
    mean_missing_sites = mean(nMissing, na.rm = TRUE)
  ) %>%
  mutate(
    gain_percent = (mean_called_sites[2] - mean_called_sites[1]) / mean_called_sites[1] * 100
  )

print(imputation_gain)



# -------------------------
# Plot 4: Imputation quality DR2/R2 vs AF
# -------------------------
imp_info <- read_tsv("data/qc/imputed.AF.DR2.tsv", show_col_types = FALSE, col_names = c("CHROM","POS","AF","R2")) %>%
  mutate(AF = as.numeric(AF), R2 = as.numeric(R2)) %>%
  filter(!is.na(AF), !is.na(R2))

if (nrow(imp_info) == 0) {
  # fallback if your tag was R2 not DR2 (rename the file in bash or change here)
  message("qc/imputed.AF.DR2.tsv empty/not found? If your tag was R2, load qc/imputed.AF.R2.tsv instead.")
}

p_r2 <- imp_info %>%
  ggplot(aes(x = AF, y = R2)) +
  geom_point(alpha = 0.15) +
  theme_minimal(base_size = 12) +
  labs(x = "Allele frequency (AF)", y = "Imputation quality (DR2/R2)",
       title = "Imputation quality vs allele frequency")

# -------------------------
# Save figures (paper-friendly)
# -------------------------
ggsave("qc/Fig_QC_depth_by_dataset.pdf", p_depth, width = 6, height = 4)
ggsave("qc/Fig_QC_residual_heterozygosity.pdf", p_het, width = 6, height = 4)
ggsave("qc/Fig_QC_het_vs_alt_burden.pdf", p_het_alt, width = 6, height = 4)
ggsave("qc/Fig_QC_imputation_R2_vs_AF.pdf", p_r2, width = 6, height = 4)
