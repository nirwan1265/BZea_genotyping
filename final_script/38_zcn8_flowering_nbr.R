#!/usr/bin/env Rscript
## =============================================================================
## 38 — One QTL window, several traits: does the teosinte allele at a window shift flowering time as well as the
##      QTL trait, and is the QTL effect still there once flowering time is in the model?
##      Written for NBR Q6 (chr8 ~126.2-126.7 Mb, ZCN8/PEBP8 Zm00001eb353250 at 126.68 Mb), but any window and traits work.
##
## Model (as in steps 26 and 37): y ~ TAXON + fixed terms + polygenic background (LOCO kinship of the window's
## chromosome). The background variance share h2 is found by maximum likelihood on a grid (0, 0.01, ..., 0.99);
## effects are generalised least squares at that h2 with Wald SEs.
##   1. effect per teosinte copy at the window, for each trait            y ~ TAXON + dos
##   2. the QTL trait adjusted for flowering time                           NBR ~ TAXON + dos + DTS  (and + DTA)
##   3. the phenotypic link between the QTL trait and flowering time       NBR ~ TAXON + DTS
##   4. effect per copy inside each taxon (and Zx subtaxon)                 y ~ TAXON + TAXON:dos
##   5. per accession: mean of carrier plots minus mean of non-carrier plots of the same accession
##
## Usage (from the repo root; <step-25 dir> as in step 26, e.g. $HOME/bzea_assoc_fix1):
##   Rscript final_script/38_zcn8_flowering_nbr.R $HOME/bzea_assoc_fix1 chr8:126650000 NBR "DTS DTA" "NBR DTS DTA PH EH" \
##     data/introgression_fix1/zcn8
## Outputs in <out dir>: zcn8_effects.tsv (1, 2, 3), zcn8_by_taxon.tsv (4), zcn8_by_accession.tsv (5)
## Needs: R with data.table.
## =============================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if (length(a) < 6) stop("usage: Rscript 38_zcn8_flowering_nbr.R <step-25 dir> <window> <QTL trait> <\"mediators\"> <\"traits\"> <out dir>")
in_dir <- a[1]; win <- a[2]; qt <- a[3]; med <- strsplit(a[4], " ")[[1]]; traits <- strsplit(a[5], " ")[[1]]; out <- a[6]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE); names(ph)[1] <- "Line"
g <- fread(file.path(in_dir, "geno_dosage.csv"), select = c(1, which(names(fread(file.path(in_dir, "geno_dosage.csv"), nrows = 0)) == win)),
           data.table = FALSE)
ph$dos <- g[[2]][match(ph$Line, g[[1]])]
ph$acc <- sub("_.*", "", ph$Line)
ph <- ph[!is.na(ph$dos), ]
chr <- sub(":.*", "", win)
k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", as.integer(sub("chr", "", chr)))), data.table = FALSE)
K <- as.matrix(k[, -1]); rownames(K) <- colnames(K) <- k[[1]]
message(sprintf("window %s: %d plots, %d carrier plots (%d SEG, %d DONOR)", win, nrow(ph), sum(ph$dos > 0),
                sum(ph$dos == 1), sum(ph$dos == 2)))

## mixed model by eigen-rotation: V = h2 K + (1 - h2) I
lmmK <- function(d, fml) {
  X <- model.matrix(fml, d); y <- d$y
  e <- eigen(K[d$Line, d$Line], symmetric = TRUE); Ut <- t(e$vectors); ty <- Ut %*% y; tX <- Ut %*% X
  n <- length(y)
  fit <- function(h) {
    w <- 1 / (h * pmax(e$values, 0) + 1 - h)
    XtW <- t(tX * w); b <- solve(XtW %*% tX, XtW %*% ty); r <- ty - tX %*% b
    s2 <- sum(w * r^2) / n
    list(ll = -0.5 * (n * log(2 * pi * s2) - sum(log(w)) + n), b = b, s2 = s2, XtWX = XtW %*% tX, h = h)
  }
  fs <- lapply(seq(0, 0.99, 0.01), fit); f <- fs[[which.max(sapply(fs, `[[`, "ll"))]]
  se <- sqrt(diag(solve(f$XtWX)) * f$s2 * n / (n - ncol(X)))
  data.frame(term = colnames(X), beta = as.vector(f$b), se = se, p = 2 * pnorm(-abs(as.vector(f$b) / se)), h2 = f$h)
}
row <- function(model, y, term, r) {
  r <- r[r$term == term, ]
  data.frame(model = model, trait = y, term = term, plots = NA, beta = signif(r$beta, 4), se = signif(r$se, 3),
             p = signif(r$p, 3), h2 = r$h2)
}
dat <- function(y, extra = NULL) { d <- ph[!is.na(ph[[y]]) & complete.cases(ph[, extra, drop = FALSE]), ]; d$y <- d[[y]]; d }

res <- list()
for (y in traits) {                                         # 1
  d <- dat(y); r <- row("1 y ~ TAXON + dos", y, "dos", lmmK(d, ~ TAXON + dos)); r$plots <- nrow(d); res[[length(res) + 1]] <- r
}
for (m in med) {                                            # 2, 3
  d <- dat(qt, m)
  r0 <- row(sprintf("2 same plots, no %s", m), qt, "dos", lmmK(d, ~ TAXON + dos))
  r1 <- lmmK(d, as.formula(sprintf("~ TAXON + dos + %s", m)))
  r <- rbind(r0, row(sprintf("2 y ~ TAXON + dos + %s", m), qt, "dos", r1), row(sprintf("2 y ~ TAXON + dos + %s", m), qt, m, r1),
             row(sprintf("3 y ~ TAXON + %s", m), qt, m, lmmK(d, as.formula(sprintf("~ TAXON + %s", m)))))
  r$plots <- nrow(d); res[[length(res) + 1]] <- r
}
E <- do.call(rbind, res)
fwrite(E, file.path(out, "zcn8_effects.tsv"), sep = "\t")

bt <- list()                                                # 4
for (y in traits) for (grp in c("TAXON", "SUBTAXON")) {
  d <- dat(y); d$grp <- factor(d[[grp]])
  if (grp == "SUBTAXON") { d$grp <- factor(ifelse(grepl("^Zx", d$SUBTAXON), d$SUBTAXON, d$TAXON)) }
  r <- lmmK(d, ~ grp + grp:dos); r <- r[grepl(":dos$", r$term), ]
  cn <- tapply(d$dos > 0, d$grp, sum); ca <- tapply(d$acc[d$dos > 0], d$grp[d$dos > 0], function(z) length(unique(z)))
  lv <- sub("^grp(.*):dos$", "\\1", r$term)
  bt[[length(bt) + 1]] <- data.frame(trait = y, by = grp, group = lv, carrier_plots = as.vector(cn[lv]),
                                     carrier_accessions = as.vector(ca[lv]), beta = signif(r$beta, 3), se = signif(r$se, 3), p = signif(r$p, 3))
}
B <- do.call(rbind, bt); B <- B[!(B$by == "SUBTAXON" & !grepl("^Zx", B$group)), ]
fwrite(B, file.path(out, "zcn8_by_taxon.tsv"), sep = "\t")

A <- rbindlist(lapply(split(ph, ph$acc), function(d) {     # 5
  if (!any(d$dos > 0) || !any(d$dos == 0)) return(NULL)
  o <- data.table(acc = d$acc[1], TAXON = d$TAXON[1], SUBTAXON = d$SUBTAXON[1], carriers = sum(d$dos > 0), noncarriers = sum(d$dos == 0))
  for (y in traits) o[[paste0("d_", y)]] <- round(mean(d[[y]][d$dos > 0], na.rm = TRUE) - mean(d[[y]][d$dos == 0], na.rm = TRUE), 3)
  o
}))
fwrite(A, file.path(out, "zcn8_by_accession.tsv"), sep = "\t")

options(width = 200)
cat("\n== effects (per teosinte copy; mediator rows per day)\n"); print(E, row.names = FALSE)
cat("\n== per taxon\n"); print(B, row.names = FALSE)
cat("\n== per accession (carrier minus non-carrier plots)\n"); print(A)
