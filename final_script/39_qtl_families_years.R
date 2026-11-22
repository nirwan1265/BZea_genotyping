#!/usr/bin/env Rscript
## =============================================================================
## 39 — Is each QTL carried by many families, and does it show in both years?
##      For every QTL that holds in the joint model (step 29: joint -log10 p > -log10(0.05 / number of QTL of the trait)):
##   Families: carrier plots (teosinte copies > 0 at the QTL window, step 25), carrier accessions and F1 families
##     (F1 family = accession + F1 plant, e.g. Zv.0200_P3), and the share of carrier plots in the largest family.
##   Effect: y ~ TAXON + dos + polygenic background (LOCO kinship of the QTL's chromosome; h2 by maximum likelihood
##     on a grid, effects by GLS, as in step 38), with all plots and with all plots of the largest carrier family
##     left out. This is the single-QTL effect (as in the scan), not the joint-model effect.
##   Families that agree: in each F1 family with carrier and non-carrier plots, carrier mean minus non-carrier mean;
##     how many of these families differ in the direction of the QTL effect.
##   Years: the same single-QTL model on each year's BLUEs alone (data/phenotypes/BZea_spats_BLUEs_<YEAR>.csv,
##     step 28): effect per copy and -log10 p in 2023 and in 2025 (traits not measured in a year are left empty).
##     The years are the same lines in two environments, not new lines, and the QTL were chosen on both years
##     together, so this checks consistency, not independent replication.
## Usage (from the repo root; step-25 inputs as in step 26):
##   Rscript final_script/39_qtl_families_years.R $HOME/bzea_assoc_fix1 data/introgression_fix1 \
##     "DTS DTA ASI NBR SPAD2 PH EH EN Prolif SL BL BW" data/introgression_fix1/checks/qtl_families_years.tsv data/phenotypes
## Needs: R with data.table.
## =============================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if (length(a) < 5) stop("usage: Rscript 39_qtl_families_years.R <step-25 dir> <results dir> <\"traits\"> <out tsv> <phenotype dir>")
in_dir <- a[1]; base <- a[2]; traits <- strsplit(a[3], " ")[[1]]; out <- a[4]; pdir <- a[5]
YP <- lapply(c(`2023` = "2023", `2025` = "2025"), function(y) fread(file.path(pdir, sprintf("BZea_spats_BLUEs_%s.csv", y)), data.table = FALSE))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE); names(ph)[1] <- "Line"
hdr <- names(fread(file.path(in_dir, "geno_dosage.csv"), nrows = 0))
Kc <- list()
getK <- function(chr) {
  if (is.null(Kc[[chr]])) {
    k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", as.integer(sub("chr", "", chr)))), data.table = FALSE)
    m <- as.matrix(k[, -1]); rownames(m) <- colnames(m) <- k[[1]]; Kc[[chr]] <<- m
  }
  Kc[[chr]]
}
lmmK <- function(d, K) {                     # y ~ TAXON + dos + polygenic background; returns beta, se, p of dos
  X <- model.matrix(~ TAXON + dos, d); y <- d$y; n <- length(y)
  e <- eigen(K[d$Line, d$Line], symmetric = TRUE); Ut <- t(e$vectors); ty <- Ut %*% y; tX <- Ut %*% X
  best <- NULL
  for (h in seq(0, 0.99, 0.01)) {
    w <- 1 / (h * pmax(e$values, 0) + 1 - h)
    XtW <- t(tX * w); A <- XtW %*% tX; b <- solve(A, XtW %*% ty); r <- ty - tX %*% b
    s2 <- sum(w * r^2) / n; ll <- -0.5 * (n * log(2 * pi * s2) - sum(log(w)) + n)
    if (is.null(best) || ll > best$ll) best <- list(ll = ll, b = b, s2 = s2, A = A)
  }
  j <- which(colnames(X) == "dos"); se <- sqrt(solve(best$A)[j, j] * best$s2 * n / (n - ncol(X)))
  c(beta = best$b[j], se = se, p = 2 * pnorm(-abs(best$b[j] / se)))
}
yr <- function(d, trait, year, K) {          # d: plots with dos; the year's BLUE of the trait as y
  P <- YP[[year]]; col <- if (trait %in% names(P)) trait else NA
  if (is.na(col)) return(c(NA, NA))
  d$y <- as.numeric(P[[col]][match(d$Line, P$new_genotype)]); d <- d[!is.na(d$y), ]
  if (nrow(d) < 100 || sum(d$dos > 0) < 10) return(c(NA, NA))
  f <- lmmK(d, K); c(-log10(f["p"]), f["beta"])
}

rows <- list()
for (tr in traits) {
  Q <- fread(file.path(base, sprintf("joint_qtl/noperm/qtl2_%s_qtl.tsv", tr)), data.table = FALSE)
  cut <- -log10(0.05 / nrow(Q))
  for (i in seq_len(nrow(Q))) {
    q <- Q[i, ]
    if (!(q$mlog10p_joint > cut) || !(q$qtl %in% hdr)) next
    g <- fread(file.path(in_dir, "geno_dosage.csv"), select = c(hdr[1], q$qtl), data.table = FALSE)
    d <- ph[!is.na(ph[[tr]]), c("Line", "TAXON", "F1FAMILY", tr)]; names(d)[4] <- "y"
    d$dos <- g[[2]][match(d$Line, g[[1]])]; d <- d[!is.na(d$dos), ]
    d$acc <- sub("_.*", "", d$Line)
    car <- d[d$dos > 0, ]
    fam_n <- sort(table(car$F1FAMILY), decreasing = TRUE)
    top <- names(fam_n)[1]
    K <- getK(q$chr)
    f_all <- lmmK(d, K); f_drop <- lmmK(d[d$F1FAMILY != top, ], K)
    sib <- sapply(split(d, d$F1FAMILY), function(x) if (any(x$dos > 0) && any(x$dos == 0))
      mean(x$y[x$dos > 0]) - mean(x$y[x$dos == 0]) else NA)
    sib <- sib[!is.na(sib)]
    y23 <- yr(d, tr, "2023", K); y25 <- yr(d, tr, "2025", K)
    rows[[length(rows) + 1]] <- data.frame(
      trait = tr, Q = paste0("Q", i), qtl = q$qtl, effect_joint = round(q$effect_per_copy_joint, 3),
      carrier_plots = nrow(car), carrier_accessions = length(unique(car$acc)), carrier_families = length(fam_n),
      largest_family = top, largest_family_share = round(as.numeric(fam_n[1]) / nrow(car), 3),
      effect_all = round(f_all["beta"], 3), mlog10p_all = round(-log10(f_all["p"]), 2),
      effect_without_largest = round(f_drop["beta"], 3), mlog10p_without_largest = round(-log10(f_drop["p"]), 2),
      families_compared = length(sib), families_same_direction = sum(sign(sib) == sign(f_all["beta"])),
      mlog10p_2023 = round(y23[1], 2), effect_2023 = round(y23[2], 3),
      mlog10p_2025 = round(y25[1], 2), effect_2025 = round(y25[2], 3), row.names = NULL)
    message(sprintf("%s Q%d %s: %d plots / %d acc / %d fam (largest %s %.0f%%); effect %.3f -> %.3f without it",
                    tr, i, q$qtl, nrow(car), length(unique(car$acc)), length(fam_n), top, 100 * fam_n[1] / nrow(car),
                    f_all["beta"], f_drop["beta"]))
  }
}
R <- do.call(rbind, rows)
fwrite(R, out, sep = "\t")
cat("wrote", out, ":", nrow(R), "QTL\n")
