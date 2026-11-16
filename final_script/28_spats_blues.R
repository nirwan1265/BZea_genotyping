#!/usr/bin/env Rscript
# =============================================================================
# 28 — Spatially corrected BLUEs (SpATS) for the BZea BC2S3 field trials, 2023 and 2025.
#
# Stage 1, per year and trait (SpATS, Rodriguez-Alvarez et al. 2018):
#     y = Rep + Border + genotype (fixed) + f(col, row) + RowF + ColF + e
#   f = 2-D P-spline surface (PSANOVA). RowF/ColF = random row and column effects.
#   Border = 1 for plots next to a "B" (border) cell in the field map (2023 only; 2025 has none).
#   One outlier pass (continuous traits only): plots with |standardised residual| > OUTLIER_SD are set to NA and the model refit.
#   Genotype as random gives the broad-sense heritability H2 (Oakey) for the summary.
#   Binary traits (ST, StPu, StPi) are not fitted: the BLUE is the line mean (share of plots = 1).
# Stage 2, across years (two-stage, weighted):
#   year effect from a weighted lm(BLUE ~ geno + year), weights 1/SE^2; combined BLUE =
#   weighted mean of the year-adjusted BLUEs, on the REF_YEAR (2023) scale. A line seen in one year keeps its (adjusted) BLUE.
#
# Inputs (RAW dir): CLY23_D4_fieldbook.csv, CLY25_B5_fieldbook.csv, CLY23_field_map.csv,
#                   CLY25_field_map.csv, name_change_final.txt (Sequencing_ID -> new_genotype)
# Outputs (OUT dir):
#   BZea_spats_BLUEs_2023.csv, _2025.csv, _combined.csv   new_genotype, traits..., Sequencing_ID
#                                                         (one row per sequenced pool; step 25 reads these)
#   BZea_spats_BLUEs_long.tsv     year, trait, new_genotype, BLUE, SE (checks included)
#   BZea_spats_summary.tsv        year, trait, plots, lines, outliers removed, H2
#   spats_plots/<year>_<trait>.pdf  SpATS diagnostic plots (raw, fitted, spatial trend, residuals)
#   spats_cache/                  one RDS per year x trait; a rerun skips cached traits (delete to refit)
#
# Run (from the repo root, R with SpATS, dplyr, tidyr):
#   Rscript final_script/28_spats_blues.R
# Options (environment):
#   PLANT_2023=2023-04-04   planting date (DTS = DOS - planting). Matches the fieldbook DTS.
#   PLANT_2025=             planting date 2025. Empty = days since 1 Jan 2025: BLUEs are then shifted by a
#                           constant, which the year effect absorbs (GWAS unchanged).
#   EXCLUDE_NOTES=regex     plots whose notes match are set to NA for all traits
#   OUTLIER_SD=3.5          0 = no outlier pass
#   TRAITS_ONLY="DTS PH"    run only these traits
#   REF_YEAR=2023           scale of the combined BLUEs (2023 has a known planting date)
#   H2=1                    0 = skip the random-genotype fit (half the run time)
# =============================================================================
suppressPackageStartupMessages({
  library(SpATS); library(dplyr); library(tidyr)
})

env <- function(k, d) { v <- Sys.getenv(k, ""); if (v == "") d else v }
RAW           <- env("RAW", "data/phenotypes/raw")
OUT           <- env("OUT", "data/phenotypes")
PLANT_2023    <- as.Date(env("PLANT_2023", "2023-04-04"))
PLANT_2025    <- env("PLANT_2025", "")
EXCLUDE_NOTES <- env("EXCLUDE_NOTES", "^NG$|\\bNG\\b|hybrid|wrong|B73\\?|afflicted|sickly")
OUTLIER_SD    <- as.numeric(env("OUTLIER_SD", "3.5"))
TRAITS_ONLY   <- strsplit(env("TRAITS_ONLY", ""), " +")[[1]]
DO_H2         <- env("H2", "1") == "1"
REF_YEAR      <- env("REF_YEAR", "2023")
BINARY        <- c("ST", "StPu", "StPi")
DISCRETE      <- c("EN", "Prolif", "Kinki", "NBR", "LAE")
CHECKS        <- c("B73", "PurpleCheck")
CACHE         <- file.path(OUT, "spats_cache")    # one RDS per year x trait; delete to refit
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
PLOTS         <- file.path(OUT, "spats_plots")     # SpATS diagnostic plot per year x trait
dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)

num <- function(x) suppressWarnings(as.numeric(x))
fix_geno <- function(x) { x <- trimws(x); x[tolower(gsub(" ", "", x)) == "purplecheck"] <- "PurpleCheck"; x }

# ---- field map: grid -> (plot, row, col, Border) --------------------------------------------------
# Row labels are in the first column, column numbers in the first line. "B" = border, "X" = alley.
read_map <- function(f) {
  g <- as.matrix(read.csv(f, header = FALSE, colClasses = "character", encoding = "UTF-8"))   # arrows/BOM left as text
  g[is.na(g)] <- ""; g <- trimws(gsub("\ufeff", "", g, useBytes = TRUE))
  is_int <- function(v) grepl("^[0-9]+$", v)
  cols <- which(is_int(g[1, ])); cols <- cols[cols > 1]
  rows <- which(is_int(g[, 1])); rows <- rows[rows > 1]
  cells <- expand.grid(i = rows, j = cols)
  cells$row  <- as.integer(g[cbind(cells$i, 1)])
  cells$col  <- as.integer(g[cbind(1, cells$j)])
  cells$cell <- g[cbind(cells$i, cells$j)]
  bord <- cells[cells$cell == "B", c("row", "col")]
  key  <- function(r, c) paste(r, c)
  bk   <- key(bord$row, bord$col)
  plots <- cells[is_int(cells$cell), ]
  plots$plot <- as.integer(plots$cell)
  plots$Border <- as.integer(key(plots$row + 1, plots$col) %in% bk | key(plots$row - 1, plots$col) %in% bk |
                             key(plots$row, plots$col + 1) %in% bk | key(plots$row, plots$col - 1) %in% bk)
  stopifnot(!anyDuplicated(plots$plot))
  plots[, c("plot", "row", "col", "Border")]
}

# ---- fieldbooks -> one row per plot ----------------------------------------------------------------
days_after <- function(d, plant, fmt) as.numeric(as.Date(d, format = fmt) - plant)

read_2023 <- function() {
  fb <- read.csv(file.path(RAW, "CLY23_D4_fieldbook.csv"), colClasses = "character", check.names = FALSE)
  d <- data.frame(plot = as.integer(fb$CLY23_D4), Rep = paste0("Rep", fb$Rep),
                  geno = fix_geno(fb$genotype_id), notes = fb$Notes,
                  DTS = days_after(fb$DOS, PLANT_2023, "%m/%d/%Y"),
                  DTA = days_after(fb$DOA, PLANT_2023, "%m/%d/%Y"))
  d$DTS[is.na(d$DTS)] <- num(fb$DTS)[is.na(d$DTS)]
  d$DTA[is.na(d$DTA)] <- num(fb$DTA)[is.na(d$DTA)]
  d$ASI <- d$DTS - d$DTA
  for (t in c("PH", "EH", "EN", "BW", "BL", "SL", "ST", "StPu", "StPi", "Kinki", "Prolif", "NBR")) d[[t]] <- num(fb[[t]])
  d$SPAD  <- num(fb$SPAD1)    # mid-July reading (closest to the 2025 timing); SPAD2 kept separately
  d$SPAD2 <- num(fb$SPAD2)
  d
}

read_2025 <- function() {
  fb <- read.csv(file.path(RAW, "CLY25_B5_fieldbook.csv"), colClasses = "character", check.names = FALSE)
  fb <- fb[grepl("^[0-9]+$", trimws(fb[[1]])), ]            # drops the "Completion:" line
  plant <- if (PLANT_2025 == "") as.Date("2025-01-01") else as.Date(PLANT_2025)
  if (PLANT_2025 == "") message("NOTE: PLANT_2025 not set: 2025 DTS/DTA are days since 1 Jan 2025 (constant shift only).")
  d <- data.frame(plot = as.integer(fb[[1]]), Rep = sub("^BZea-", "", fb[["Who/What"]]),
                  geno = fix_geno(fb$Description), notes = paste(fb$notes1, fb$notes2),
                  DTS = days_after(fb$DOS, plant, "%m/%d/%y"),
                  DTA = days_after(fb$DOA, plant, "%m/%d/%y"))
  d$ASI <- d$DTS - d$DTA
  for (t in c("PH", "EH", "EN", "Prolif", "StPi", "StPu", "NBR", "LAE", "SPAD")) d[[t]] <- num(fb[[t]])
  d
}

# ---- SpATS ---------------------------------------------------------------------------------------
nseg_for <- function(d) c(max(8, min(20, round(length(unique(d$col)) / 3))),
                          max(8, min(20, round(length(unique(d$row)) / 3))))

fit_spats <- function(d, trait, geno_random, use_rep = TRUE, use_border = TRUE) {
  d$Rep <- droplevels(d$Rep); d$Border <- droplevels(d$Border)
  terms <- c(if (use_rep && nlevels(d$Rep) > 1) "Rep", if (use_border && nlevels(d$Border) > 1) "Border")
  fx <- if (length(terms)) paste("~", paste(terms, collapse = " + ")) else NULL
  if (!is.null(fx)) fx <- as.formula(fx)
  SpATS(response = trait, genotype = "geno", genotype.as.random = geno_random,
        spatial = ~ PSANOVA(col, row, nseg = nseg_for(d)),
        fixed = fx, random = ~ RowF + ColF, data = d,
        control = list(tolerance = 1e-03, monitoring = 0))
}

run_year <- function(year, d, map) {
  d <- merge(d, map, by = "plot", all.x = TRUE)
  if (any(is.na(d$row))) stop(year, ": ", sum(is.na(d$row)), " fieldbook plots not in the field map")
  bad <- grepl(EXCLUDE_NOTES, d$notes, ignore.case = TRUE) & !is.na(d$notes)
  message(year, ": ", nrow(d), " plots, ", length(unique(d$geno)), " lines, ", sum(bad),
          " plots excluded by notes (", paste(head(unique(trimws(d$notes[bad])), 8), collapse = "; "), ")")
  traits <- setdiff(names(d), c("plot", "Rep", "geno", "notes", "row", "col", "Border"))
  if (length(TRAITS_ONLY) && TRAITS_ONLY[1] != "") traits <- intersect(traits, TRAITS_ONLY)
  d[bad, traits] <- NA
  d$Rep <- factor(d$Rep); d$Border <- factor(d$Border); d$geno <- factor(d$geno)
  d$RowF <- factor(d$row); d$ColF <- factor(d$col)

  blues <- list(); summ <- list()
  for (t in traits) {
    dt <- d[!is.na(d[[t]]), ]
    if (nrow(dt) < 100 || length(unique(dt[[t]])) < 2) { message("  ", t, ": skipped (too few values)"); next }
    if (t %in% BINARY) {
      b <- dt %>% group_by(geno) %>% summarise(BLUE = mean(.data[[t]]), SE = sd(.data[[t]]) / sqrt(n()), .groups = "drop")
      b$SE[is.na(b$SE) | b$SE == 0] <- NA
      blues[[t]] <- data.frame(year = as.character(year), trait = t, new_genotype = as.character(b$geno), BLUE = b$BLUE, SE = b$SE)
      summ[[t]] <- data.frame(year = as.character(year), trait = t, plots = nrow(dt), lines = nrow(b), outliers = 0, H2 = NA,
                              method = "line mean (binary)")
      message(sprintf("  %-6s line mean (binary)", t)); next
    }
    cache <- file.path(CACHE, sprintf("%s_%s.rds", year, t))
    if (file.exists(cache)) { r <- readRDS(cache); blues[[t]] <- r$b; summ[[t]] <- r$s
      message(sprintf("  %-6s (cached)", t)); next }
    t0 <- Sys.time()
    use_rep <- TRUE
    # Rep can be confounded with genotype when a trait was scored mostly in one rep: then drop Rep
    try_fit <- function(dt) {
      f <- tryCatch(fit_spats(dt, t, FALSE, TRUE), error = function(e) NULL)
      if (!is.null(f)) return(list(fit = f, rep = TRUE, border = TRUE))
      f <- tryCatch(fit_spats(dt, t, FALSE, FALSE), error = function(e) NULL)
      if (!is.null(f)) return(list(fit = f, rep = FALSE, border = TRUE))
      # scored in one rep with ~1 plot per line: Border is then only estimable from checks; drop it too
      f <- tryCatch(fit_spats(dt, t, FALSE, FALSE, FALSE), error = function(e) { message("  ", t, ": SpATS failed: ", conditionMessage(e)); NULL })
      list(fit = f, rep = FALSE, border = FALSE)
    }
    tf <- try_fit(dt); fit <- tf$fit; use_rep <- tf$rep; use_border <- tf$border
    if (is.null(fit)) next
    n_out <- 0
    if (OUTLIER_SD > 0 && !t %in% DISCRETE) {   # count/score traits: a robust SD is ~0, no outlier pass
      r <- residuals(fit); z <- r / (1.4826 * mad(r, na.rm = TRUE))
      out <- which(abs(z) > OUTLIER_SD)
      if (length(out)) { n_out <- length(out); dt <- dt[-out, ]; dt$geno <- droplevels(dt$geno)
                        tf <- try_fit(dt); if (is.null(tf$fit)) next; fit <- tf$fit; use_rep <- tf$rep; use_border <- tf$border }
    }
    p <- predict(fit, which = "geno")
    blues[[t]] <- data.frame(year = as.character(year), trait = t, new_genotype = as.character(p$geno),
                             BLUE = p$predicted.values, SE = p$standard.errors)
    h2 <- NA
    if (DO_H2) h2 <- tryCatch(unname(getHeritability(fit_spats(dt, t, TRUE, use_rep, use_border))), error = function(e) NA)
    pdf(file.path(PLOTS, sprintf("%s_%s.pdf", year, t)), width = 11, height = 8)
    try(plot(fit, main = paste(year, t)), silent = TRUE); dev.off()
    summ[[t]] <- data.frame(year = as.character(year), trait = t, plots = nrow(dt), lines = nrow(p), outliers = n_out,
                            H2 = round(h2, 3), method = paste0("SpATS", if (!use_rep) ", no Rep" else "", if (!use_border) ", no Border" else ""))
    saveRDS(list(b = blues[[t]], s = summ[[t]]), cache)
    message(sprintf("  %-6s plots %4d  lines %4d  outliers %3d  H2 %s  (%.0f s)", t, nrow(dt), nrow(p), n_out,
                    ifelse(is.na(h2), "NA", sprintf("%.2f", h2)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  list(blues = bind_rows(blues), summary = bind_rows(summ))
}

# ---- stage 2: combine years ----------------------------------------------------------------------
combine_years <- function(L) {
  out <- list()
  for (t in unique(L$trait)) {
    b <- L[L$trait == t & !is.na(L$BLUE), ]
    yrs <- unique(b$year)
    if (length(yrs) == 1) { out[[t]] <- b %>% transmute(trait, new_genotype, BLUE, SE, years = 1L); next }
    b$w <- ifelse(is.na(b$SE) | b$SE <= 0, 1 / median(b$SE^2, na.rm = TRUE), 1 / b$SE^2)
    if (all(is.na(b$SE))) b$w <- 1
    m <- lm(BLUE ~ new_genotype + factor(year), data = b, weights = w)
    ye <- c(0, coef(m)[paste0("factor(year)", yrs[-1])]); names(ye) <- yrs
    ye[is.na(ye)] <- 0; ye <- ye - ye[REF_YEAR]            # combined BLUEs on the REF_YEAR scale
    b$adj <- b$BLUE - ye[as.character(b$year)]
    out[[t]] <- b %>% group_by(new_genotype) %>%
      summarise(trait = t, BLUE = sum(w * adj) / sum(w), SE = sqrt(1 / sum(w)), years = n(), .groups = "drop") %>%
      select(trait, new_genotype, BLUE, SE, years)
    message(sprintf("  %-6s year effects: %s", t, paste(sprintf("%s %+.2f", names(ye), ye), collapse = ", ")))
  }
  bind_rows(out)
}

# ---- write one wide file per scope, one row per Sequencing_ID ---------------------------------------
nc <- read.delim(file.path(RAW, "name_change_final.txt"), colClasses = "character")
nc$new_genotype <- fix_geno(nc$new_genotype)
write_wide <- function(L, file) {
  w <- L %>% filter(!new_genotype %in% CHECKS) %>% select(new_genotype, trait, BLUE) %>%
    pivot_wider(names_from = trait, values_from = BLUE)
  miss <- setdiff(w$new_genotype, nc$new_genotype)
  if (length(miss)) message("  ", basename(file), ": ", length(miss), " phenotyped lines without a Sequencing_ID: ",
                            paste(head(miss, 5), collapse = ", "))
  w <- inner_join(w, nc[!nc$new_genotype %in% CHECKS, c("new_genotype", "Sequencing_ID")], by = "new_genotype")
  write.csv(w, file, row.names = FALSE, na = "")
  message("WROTE ", file, " (", nrow(w), " pools, traits: ", paste(setdiff(names(w), c("new_genotype", "Sequencing_ID")), collapse = " "), ")")
}

# ---- main ---------------------------------------------------------------------------------------
r23 <- run_year(2023, read_2023(), read_map(file.path(RAW, "CLY23_field_map.csv")))
r25 <- run_year(2025, read_2025(), read_map(file.path(RAW, "CLY25_field_map.csv")))
L <- bind_rows(r23$blues, r25$blues)
message("Combining years:")
C <- combine_years(L)

write.table(bind_rows(L, C %>% mutate(year = "combined")) %>% select(year, trait, new_genotype, BLUE, SE, years),
            file.path(OUT, "BZea_spats_BLUEs_long.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
# correlation between years per trait (lines in both years)
cy <- L %>% select(year, trait, new_genotype, BLUE) %>% pivot_wider(names_from = year, values_from = BLUE) %>%
  filter(!new_genotype %in% CHECKS) %>% group_by(trait) %>%
  summarise(r_2023_2025 = tryCatch(round(cor(`2023`, `2025`, use = "complete.obs"), 2), error = function(e) NA_real_),
            .groups = "drop")
S <- left_join(bind_rows(r23$summary, r25$summary), cy, by = "trait")
write.table(S, file.path(OUT, "BZea_spats_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
print(S, row.names = FALSE)

write_wide(r23$blues, file.path(OUT, "BZea_spats_BLUEs_2023.csv"))
write_wide(r25$blues, file.path(OUT, "BZea_spats_BLUEs_2025.csv"))
write_wide(C,         file.path(OUT, "BZea_spats_BLUEs_combined.csv"))
