#!/usr/bin/env Rscript
# =============================================================================
# 28b — Spatially corrected BLUEs (SpATS) for predicted leaf nitrogen, 2025 field B5, block 2.
#
# Input: data/phenotypes/B5_block2_predictedN.csv, one row per scan (5 per plot):
#   Plot, Row, Range, Predicted_N, M_Distance, Genotype, Accession, Taxa
#   Predicted_N = nitrogen predicted from the leaf spectrum; M_Distance = Mahalanobis distance of the
#   spectrum to the calibration set (large = spectrum unlike the calibration, prediction unreliable).
# Steps:
#   1. Scans with M_Distance > MDIST (default 3, the usual cut for spectral outliers) are dropped.
#   2. Plot value = mean of the remaining scans.
#   3. SpATS, as in step 28: N = genotype (fixed) + f(Range, Row) + RowF + RangeF + e (PSANOVA surface);
#      one outlier pass (|robust z| > OUTLIER_SD, on the genotype-random fit) and refit; genotype random for H2 (Oakey).
#      B73 (125 plots) and the purple check are in the fit and give the spatial trend; they are not written out.
# Output (data/phenotypes/):
#   BZea_spats_BLUEs_N_2025.csv   new_genotype, PredN, Sequencing_ID (one row per sequenced pool; step 25 reads it)
#   spats_plots/2025_PredN.pdf    SpATS diagnostic plots
#   BZea_spats_summary_N.tsv      plots, lines, scans dropped, outliers removed, H2
# Run (from the repo root, R with SpATS, dplyr):
#   Rscript final_script/28b_spats_nitrogen.R
# Options (environment): IN, OUT, MDIST=3, OUTLIER_SD=3.5
# =============================================================================
suppressPackageStartupMessages({ library(SpATS); library(dplyr) })
env <- function(k, d) { v <- Sys.getenv(k, ""); if (v == "") d else v }
IN         <- env("IN", "data/phenotypes/B5_block2_predictedN.csv")
OUT        <- env("OUT", "data/phenotypes")
RAW        <- env("RAW", "data/phenotypes/raw")
MDIST      <- as.numeric(env("MDIST", "3"))
OUTLIER_SD <- as.numeric(env("OUTLIER_SD", "3.5"))
CHECKS     <- c("B73", "PurpleCheck")
fix_geno <- function(x) { x <- trimws(x); x[tolower(gsub(" ", "", x)) == "purplecheck"] <- "PurpleCheck"; x }

s <- read.csv(IN, stringsAsFactors = FALSE)
s$geno <- fix_geno(s$Genotype)
n_scan <- nrow(s); drop <- s$M_Distance > MDIST
message(sprintf("%d scans on %d plots; %d scans (%.1f%%) dropped with M_Distance > %g",
                n_scan, length(unique(s$Plot)), sum(drop), 100 * mean(drop), MDIST))
d <- s[!drop, ] %>% group_by(Plot, Row, Range, geno) %>%
  summarise(PredN = mean(Predicted_N), scans = n(), .groups = "drop") %>% as.data.frame()
stopifnot(!anyDuplicated(d$Plot))
d$geno <- factor(d$geno); d$RowF <- factor(d$Row); d$RangeF <- factor(d$Range)
message(sprintf("%d plots, %d lines (B73 %d plots)", nrow(d), nlevels(d$geno), sum(d$geno == "B73")))

nseg <- c(max(8, min(20, round(length(unique(d$Range)) / 3))), max(8, min(20, round(length(unique(d$Row)) / 3))))
fit_spats <- function(d, geno_random)
  SpATS(response = "PredN", genotype = "geno", genotype.as.random = geno_random,
        spatial = ~ PSANOVA(Range, Row, nseg = nseg), random = ~ RowF + RangeF, data = d,
        control = list(tolerance = 1e-03, monitoring = 0))
# Outlier pass on the genotype-random fit: most lines have one plot, so with genotype fixed their residuals
# are exactly 0 and only the checks would have residuals to judge.
fr <- fit_spats(d, TRUE)
r <- residuals(fr); z <- r / (1.4826 * mad(r, na.rm = TRUE)); out <- which(abs(z) > OUTLIER_SD)
message(length(out), " plots with |robust z| > ", OUTLIER_SD, if (length(out)) paste0(": ",
        paste(sprintf("%s (%.2f)", d$geno[out], d$PredN[out]), collapse = ", ")) else "")
if (length(out)) { d <- d[-out, ]; d$geno <- droplevels(d$geno); d$RowF <- droplevels(d$RowF); d$RangeF <- droplevels(d$RangeF)
                   fr <- fit_spats(d, TRUE) }
fit <- fit_spats(d, FALSE)
p <- predict(fit, which = "geno")
h2 <- tryCatch(unname(getHeritability(fr)), error = function(e) NA)
dir.create(file.path(OUT, "spats_plots"), showWarnings = FALSE, recursive = TRUE)
pdf(file.path(OUT, "spats_plots", "2025_PredN.pdf"), width = 11, height = 8); plot(fit, main = "2025 PredN"); dev.off()

nc <- read.delim(file.path(RAW, "name_change_final.txt"), colClasses = "character")
nc$new_genotype <- fix_geno(nc$new_genotype)
b <- data.frame(new_genotype = as.character(p$geno), PredN = p$predicted.values) %>% filter(!new_genotype %in% CHECKS)
miss <- setdiff(b$new_genotype, nc$new_genotype)
if (length(miss)) message(length(miss), " lines without a Sequencing_ID: ", paste(head(miss, 5), collapse = ", "))
w <- inner_join(b, nc[, c("new_genotype", "Sequencing_ID")], by = "new_genotype")
write.csv(w, file.path(OUT, "BZea_spats_BLUEs_N_2025.csv"), row.names = FALSE)
S <- data.frame(year = "2025", trait = "PredN", scans = n_scan, scans_dropped_Mdist = sum(drop), plots = nrow(d),
                lines = nrow(p), outliers = length(out), H2 = round(h2, 3), pools_written = nrow(w))
write.table(S, file.path(OUT, "BZea_spats_summary_N.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(S, row.names = FALSE)
message("WROTE ", file.path(OUT, "BZea_spats_BLUEs_N_2025.csv"), " (", nrow(w), " pools)")
