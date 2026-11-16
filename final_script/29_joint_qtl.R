#!/usr/bin/env Rscript
## =============================================================================
## 29 — Joint multi-QTL model with per-taxon QTL effects (k-mer ancestry, LOCO kinship).
##
## Model:  y = taxon + sum_q sum_t b_qt * dosage_q * [taxon = t] + u + e,   u ~ N(0, s2g * K_LOCO), e ~ N(0, s2e * I)
##   dosage 0/1/2 = B73/SEG/DONOR per 100-kb window (step 25). b_qt = effect per teosinte copy of QTL q in taxon t.
##
## 1. Threshold: phenotypes shuffled within taxon, genome scan rerun, largest -log10 p kept; 95th percentile
##    over NPERM permutations (or THRESHOLD=, or the step-26 summary with SUMMARY=).
## 2. Forward selection (multi-locus mixed model, Segura et al. 2012): scan with the selected QTL as cofactors
##    (one pooled effect each), add the top window if it passes the threshold, repeat. Windows within MIN_DIST
##    (5 Mb) of, or with r2 > R2_MAX (0.5) to, a selected QTL are not eligible: neighbouring windows differ in
##    only a few plots, and a pair of them fits those plots with huge opposite effects.
##    The scan uses, for a window on chromosome c, the LOCO K without c; the variance ratio is fitted by ML
##    once per chromosome and step and held fixed across that chromosome's windows (P3D / EMMAX, Kang et al. 2010).
## 3. Backward: the weakest QTL is dropped while its conditional test (exact ML, LOCO K) is below the threshold.
## 4. Final model, per QTL q (LOCO K of q's chromosome, exact ML):
##    - presence: all QTL nested in taxon vs q removed (LRT, df = number of taxon effects of q)
##    - heterogeneity: q nested in taxon vs one pooled effect for q (LRT, df = taxon effects - 1)
##    - pooled effect and SE; per-taxon effects and SEs; PVE (variance of q's fitted contribution / variance of y)
##    - support interval: q's chromosome rescanned with the other QTL as cofactors; windows within DROP (1.5)
##      -log10 p of the top window; B73 v5 genes in it when GFF is given.
##    A taxon gets its own effect at q only with >= MIN_TAXON_CARRIERS carrier plots (dosage > 0); carriers of
##    smaller taxa share one "other" effect.
##
## Usage:  Rscript 29_joint_qtl.R <assoc input dir (step 25)> <TRAIT> <out dir>
## Env:    NPERM (1000) THRESHOLD SUMMARY GFF MAX_QTL (20) MIN_TAXON_CARRIERS (5) DROP (1.5) MIN_CARRIERS (20) SEED (1)
##         R2_MAX (0.5), MIN_DIST (5e6): a new QTL must have r2 <= R2_MAX with, and be >= MIN_DIST bp from, every selected QTL
## Needs:  R >= 4.1 with data.table. No GridLMM.
## Outputs (<out dir>):
##   joint_qtl_<TRAIT>.tsv           one row per QTL: position, interval, pooled effect, p presence, p heterogeneity, PVE
##   joint_qtl_<TRAIT>_by_taxon.tsv  QTL x taxon: carriers (SEG/DONOR), effect per copy, SE, p
##   joint_qtl_<TRAIT>_genes.tsv     B73 v5 genes in each QTL interval (with GFF)
##   joint_qtl_<TRAIT>_steps.tsv     forward/backward log;  joint_qtl_<TRAIT>_perm.tsv  permutation maxima
## =============================================================================
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
if (length(args) < 3) stop("usage: Rscript 29_joint_qtl.R <assoc input dir> <TRAIT> <out dir>")
in_dir <- args[1]; trait <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
envn <- function(k, d) { v <- Sys.getenv(k, ""); if (v == "") d else as.numeric(v) }
NPERM <- envn("NPERM", 1000); MAX_QTL <- envn("MAX_QTL", 20); MIN_TC <- envn("MIN_TAXON_CARRIERS", 5)
DROP <- envn("DROP", 1.5); R2_MAX <- envn("R2_MAX", 0.5); MIN_DIST <- envn("MIN_DIST", 5e6); MIN_CAR <- envn("MIN_CARRIERS", 20); set.seed(envn("SEED", 1))
THRESHOLD <- Sys.getenv("THRESHOLD", ""); SUMMARY <- Sys.getenv("SUMMARY", ""); GFF <- Sys.getenv("GFF", "")
t_start <- Sys.time()
say <- function(...) message(sprintf(...))

## ---------- inputs ----------
ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE)
if (!trait %in% names(ph)) stop("trait ", trait, " not in pheno.csv")
g <- fread(file.path(in_dir, "geno_dosage.csv"), data.table = FALSE)
G <- as.matrix(g[, -1]); storage.mode(G) <- "double"; rownames(G) <- g[[1]]; rm(g)
map <- fread(file.path(in_dir, "window_map.tsv"), data.table = FALSE)
stopifnot(identical(colnames(G), map$Marker))
dat <- data.frame(Line = ph[[1]], y = as.numeric(ph[[trait]]), TAXON = as.character(ph$TAXON), stringsAsFactors = FALSE)
dat <- dat[!is.na(dat$y) & dat$Line %in% rownames(G), ]
G <- G[dat$Line, , drop = FALSE]
keep <- colSums(G > 0) >= MIN_CAR & apply(G, 2, function(x) length(unique(x)) > 1)
G <- G[, keep, drop = FALSE]; map <- map[keep, ]
n <- nrow(dat); taxa <- sort(unique(dat$TAXON))
X0 <- model.matrix(~ TAXON, data = dat)
chrs <- unique(map$chr); win_chr <- map$chr
say("joint QTL | trait %s | plots %d | windows %d | taxa %s", trait, n, ncol(G),
    paste(sprintf("%s %d", taxa, as.integer(table(dat$TAXON)[taxa])), collapse = ", "))

## ---------- LOCO K: eigen-decomposition and rotated genotypes per chromosome ----------
E <- list(); GR <- list()
for (ch in chrs) {
  k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", as.integer(sub("chr", "", ch)))), data.table = FALSE)
  K <- as.matrix(k[, -1]); rownames(K) <- colnames(K) <- k[[1]]
  K <- K[dat$Line, dat$Line]
  e <- eigen(K, symmetric = TRUE); e$values[e$values < 0] <- 0
  E[[ch]] <- list(U = e$vectors, d = e$values)
  GR[[ch]] <- crossprod(e$vectors, G[, win_chr == ch, drop = FALSE])     # U' G_c
}
say("eigen-decompositions done (%.0f s)", as.numeric(difftime(Sys.time(), t_start, units = "secs")))

## ---------- ML for y = X b + u + e with K = U diag(d) U' ----------
## profiled log-likelihood in delta = s2e/s2g:  -n/2 log(RSS/n) - 1/2 sum log(d + delta)   (constants dropped)
gls <- function(yr, Xr, w) {
  sw <- sqrt(w); Xs <- Xr * sw; ys <- yr * sw
  qx <- qr(Xs); b <- qr.coef(qx, ys); b0 <- b; b0[is.na(b0)] <- 0
  r <- ys - Xs %*% b0
  list(b = b, rss = sum(r^2), qx = qx, Xs = Xs, ys = ys)
}
loglik <- function(ld, yr, Xr, d) {
  delta <- exp(ld); w <- 1 / (d + delta)
  f <- gls(yr, Xr, w)
  -n / 2 * log(f$rss / n) - 0.5 * sum(log(d + delta))
}
fit_ml <- function(y, X, ch) {
  U <- E[[ch]]$U; d <- E[[ch]]$d
  yr <- crossprod(U, y); Xr <- crossprod(U, X)
  o <- optimize(loglik, c(-10, 10), yr = yr, Xr = Xr, d = d, maximum = TRUE)
  delta <- exp(o$maximum); w <- 1 / (d + delta); f <- gls(yr, Xr, w)
  s2 <- f$rss / n                                  # ML: s2g (times K) ; s2e = delta * s2g
  ok <- !is.na(f$b); se <- rep(NA_real_, length(f$b))
  se[ok] <- sqrt(pmax(diag(solve(crossprod(f$Xs[, ok, drop = FALSE]))) * s2, 0))
  list(ll = o$objective, b = setNames(drop(f$b), colnames(X)), se = setNames(se, colnames(X)), delta = delta,
       h2 = 1 / (1 + delta), s2g = s2, yr = yr, Xr = Xr, w = w)
}
MASS_ginv <- function(A) { s <- svd(A); tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol; s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos]) }

## ---------- genome scan with cofactors (variance ratio fixed per chromosome) ----------
## returns -log10 p (1-df LRT, profiled) for every window
scan_genome <- function(y, X) {
  out <- numeric(ncol(G))
  for (ch in chrs) {
    f <- fit_ml(y, X, ch); sw <- sqrt(f$w)
    Xs <- f$Xr * sw; ys <- drop(f$yr) * sw; Gs <- GR[[ch]] * sw
    qx <- qr(Xs); ry <- qr.resid(qx, ys); RG <- qr.resid(qx, Gs)
    gg <- colSums(RG^2); gy <- drop(crossprod(RG, ry)); rss0 <- sum(ry^2)
    gg0 <- colSums(Gs^2)
    rss1 <- rss0 - gy^2 / gg
    lrt <- n * log(rss0 / pmax(rss1, 1e-12))
    lrt[gg < 1e-8 * gg0 | !is.finite(lrt)] <- 0                     # window collinear with a cofactor
    out[win_chr == ch] <- -pchisq(lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  }
  out
}

## ---------- threshold ----------
thr <- NA
if (THRESHOLD != "") { thr <- as.numeric(THRESHOLD); say("threshold given: %.2f", thr)
} else if (SUMMARY != "" && file.exists(SUMMARY)) {
  s <- fread(SUMMARY, data.table = FALSE); thr <- as.numeric(s$threshold[s$scope == "ALL"])
  say("threshold from %s: %.2f", SUMMARY, thr)
}
perm_file <- file.path(out_dir, sprintf("joint_qtl_%s_perm.tsv", trait))
if (is.na(thr)) {
  pm <- if (file.exists(perm_file)) fread(perm_file)$max_mlog10p else numeric(0)
  todo <- NPERM - length(pm)
  if (todo > 0) {
    say("%d permutations (within taxon) ...", todo); t0 <- Sys.time()
    idx <- split(seq_len(n), dat$TAXON)
    for (i in seq_len(todo)) {
      yp <- dat$y; for (ix in idx) yp[ix] <- dat$y[ix][sample.int(length(ix))]
      pm <- c(pm, max(scan_genome(yp, X0)))
      if (i %% 50 == 0 || i == todo) {
        fwrite(data.table(perm = seq_along(pm), max_mlog10p = pm), perm_file, sep = "\t")
        say("  %d / %d  (%.1f s per permutation)", i, todo, as.numeric(difftime(Sys.time(), t0, units = "secs")) / i)
      }
    }
  }
  thr <- unname(quantile(pm, 0.95)); say("5%% genome-wide threshold from %d permutations: %.2f", length(pm), thr)
}

## ---------- forward / backward selection ----------
steps <- list(); Q <- character(0)
dose <- function(m) G[, m]
Xq <- function(qs) if (length(qs)) cbind(X0, G[, qs, drop = FALSE]) else X0
chr_of <- function(m) map$chr[match(m, map$Marker)]
for (k in seq_len(MAX_QTL)) {
  p <- scan_genome(dat$y, Xq(Q))
  if (length(Q)) {                  # a new QTL must not be a near-copy of a selected one (collinear neighbours)
    r2 <- suppressWarnings(cor(G, G[, Q, drop = FALSE])^2); r2[is.na(r2)] <- 0
    near <- sapply(Q, function(q) map$chr == chr_of(q) & abs(map$pos - map$pos[match(q, map$Marker)]) < MIN_DIST)
    p[apply(r2, 1, max) > R2_MAX | rowSums(as.matrix(near)) > 0] <- 0
  }
  top <- which.max(p)
  steps[[length(steps) + 1]] <- data.frame(step = length(steps) + 1, action = if (p[top] > thr) "add" else "stop",
                                           window = map$Marker[top], mlog10p = round(p[top], 2), threshold = round(thr, 2))
  say("forward %d: top %s  -log10 p %.2f", k, map$Marker[top], p[top])
  if (p[top] <= thr) break
  Q <- c(Q, map$Marker[top])
}
cond_p <- function(qs, q) {         # exact ML LRT for q given the other QTL, LOCO K of q's chromosome
  ch <- chr_of(q); f1 <- fit_ml(dat$y, Xq(qs), ch); f0 <- fit_ml(dat$y, Xq(setdiff(qs, q)), ch)
  -pchisq(2 * (f1$ll - f0$ll), 1, lower.tail = FALSE, log.p = TRUE) / log(10)
}
repeat {
  if (!length(Q)) break
  cp <- sapply(Q, function(q) cond_p(Q, q)); w <- which.min(cp)
  if (cp[w] > thr) break
  steps[[length(steps) + 1]] <- data.frame(step = length(steps) + 1, action = "drop", window = Q[w],
                                           mlog10p = round(cp[w], 2), threshold = round(thr, 2))
  say("backward: drop %s (conditional -log10 p %.2f)", Q[w], cp[w]); Q <- Q[-w]
}
fwrite(rbindlist(steps), file.path(out_dir, sprintf("joint_qtl_%s_steps.tsv", trait)), sep = "\t")
if (!length(Q)) { say("no QTL above the threshold (%.2f)", thr); quit(save = "no") }
Q <- Q[order(match(chr_of(Q), chrs), map$start[match(Q, map$Marker)])]
say("QTL in the joint model: %s", paste(Q, collapse = ", "))

## ---------- per-taxon design ----------
nested_cols <- function(q) {        # dosage of q split by taxon; small taxa pooled as "other"
  x <- dose(q); car <- tapply(x > 0, dat$TAXON, sum)
  own <- names(car)[car >= MIN_TC]; rest <- setdiff(names(car)[car > 0], own)
  M <- sapply(own, function(t) x * (dat$TAXON == t))
  if (!is.matrix(M)) M <- matrix(M, ncol = length(own))
  colnames(M) <- if (length(own)) paste0(q, "|", own) else character(0)
  if (length(rest)) { o <- x * (dat$TAXON %in% rest); M <- cbind(M, o); colnames(M)[ncol(M)] <- paste0(q, "|other(", paste(rest, collapse = "+"), ")") }
  M
}
NC <- setNames(lapply(Q, nested_cols), Q)
X_full <- function(except = NULL, pooled = NULL) {
  X <- X0
  for (q in Q) {
    if (identical(q, except)) next
    if (identical(q, pooled)) { X <- cbind(X, G[, q, drop = FALSE]); next }
    X <- cbind(X, NC[[q]])
  }
  X
}

## ---------- final model per QTL ----------
yv <- var(dat$y)
rows <- list(); trows <- list(); grows <- list()
genes <- NULL
if (GFF != "") {
  gf <- fread(cmd = sprintf("grep -v '^#' %s | awk -F'\\t' '$3==\"gene\"'", shQuote(GFF)), header = FALSE, sep = "\t", quote = "")
  if (nrow(gf)) {
    a <- gf$V9
    genes <- data.frame(chr = ifelse(grepl("^chr", gf$V1), gf$V1, paste0("chr", gf$V1)), start = gf$V4, end = gf$V5,
                        id = sub(".*(gene_id|ID)=(gene:)?([^;]+).*", "\\3", a),
                        name = ifelse(grepl("Name=", a), sub(".*Name=([^;]+).*", "\\1", a), ""),
                        desc = ifelse(grepl("description=", a), utils::URLdecode(sub(".*description=([^;]+).*", "\\1", a)), ""))
  }
}
for (q in Q) {
  ch <- chr_of(q)
  fF <- fit_ml(dat$y, X_full(), ch); fR <- fit_ml(dat$y, X_full(except = q), ch); fP <- fit_ml(dat$y, X_full(pooled = q), ch)
  kq <- ncol(NC[[q]])
  p_pres <- -pchisq(2 * (fF$ll - fR$ll), kq, lower.tail = FALSE, log.p = TRUE) / log(10)
  p_het <- if (kq > 1) -pchisq(2 * (fF$ll - fP$ll), kq - 1, lower.tail = FALSE, log.p = TRUE) / log(10) else NA
  bq <- fF$b[colnames(NC[[q]])]; bq[is.na(bq)] <- 0
  contrib <- drop(NC[[q]] %*% bq)
  pve <- var(contrib) / yv
  ## support interval: rescan q's chromosome with the other QTL as pooled cofactors
  others <- setdiff(Q, q); Xo <- Xq(others)
  f <- fit_ml(dat$y, Xo, ch); sw <- sqrt(f$w); qx <- qr(f$Xr * sw); ry <- qr.resid(qx, drop(f$yr) * sw)
  Gs <- GR[[ch]] * sw; RG <- qr.resid(qx, Gs); gg <- colSums(RG^2); gy <- drop(crossprod(RG, ry)); rss0 <- sum(ry^2)
  lrt <- n * log(rss0 / pmax(rss0 - gy^2 / gg, 1e-12)); lrt[!is.finite(lrt)] <- 0
  pc <- -pchisq(lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  mc <- map[win_chr == ch, ]; top <- which.max(pc); lo <- top; hi <- top
  while (lo > 1 && pc[lo - 1] >= pc[top] - DROP) lo <- lo - 1
  while (hi < length(pc) && pc[hi + 1] >= pc[top] - DROP) hi <- hi + 1
  bP <- fP$b[q]; seP <- fP$se[q]
  rows[[q]] <- data.frame(trait = trait, qtl = q, chr = ch, pos = map$pos[match(q, map$Marker)],
                          interval_start = mc$start[lo], interval_end = mc$end[hi], conditional_top = mc$Marker[top],
                          conditional_top_mlog10p = round(pc[top], 2), carriers = sum(dose(q) > 0),
                          effect_per_copy = round(bP, 3), se = round(seP, 3), mlog10p_presence = round(p_pres, 2),
                          df_presence = kq, mlog10p_heterogeneity = round(p_het, 2), pve = round(pve, 4),
                          h2_polygenic = round(fF$h2, 3), threshold = round(thr, 2))
  for (cn in colnames(NC[[q]])) {
    tx <- sub("^[^|]*\\|", "", cn); x <- NC[[q]][, cn]; i <- match(cn, names(fF$b))
    z <- fF$b[i] / fF$se[i]
    trows[[cn]] <- data.frame(trait = trait, qtl = q, taxon = tx, carriers = sum(x > 0), seg = sum(x == 1), donor = sum(x == 2),
                              effect_per_copy = round(fF$b[i], 3), se = round(fF$se[i], 3),
                              mlog10p = if (is.na(z)) NA else round(-pchisq(z^2, 1, lower.tail = FALSE, log.p = TRUE) / log(10), 2))
  }
  if (!is.null(genes)) {
    gq <- genes[genes$chr == ch & genes$end >= mc$start[lo] & genes$start <= mc$end[hi], ]
    if (nrow(gq)) grows[[q]] <- data.frame(trait = trait, qtl = q, chr = ch, interval_start = mc$start[lo],
                                           interval_end = mc$end[hi], gene = gq$id, gene_start = gq$start, gene_end = gq$end,
                                           name = gq$name, description = gq$desc)
  }
}
R <- rbindlist(rows); TB <- rbindlist(trows)
ftot <- fit_ml(dat$y, X_full(), chrs[1])
bt <- ftot$b[-seq_len(ncol(X0))]; bt[is.na(bt)] <- 0
tot <- drop(X_full()[, -seq_len(ncol(X0)), drop = FALSE] %*% bt)
fwrite(R, file.path(out_dir, sprintf("joint_qtl_%s.tsv", trait)), sep = "\t")
fwrite(TB, file.path(out_dir, sprintf("joint_qtl_%s_by_taxon.tsv", trait)), sep = "\t")
if (length(grows)) fwrite(rbindlist(grows), file.path(out_dir, sprintf("joint_qtl_%s_genes.tsv", trait)), sep = "\t")
print(as.data.frame(R[, .(qtl, interval_start, interval_end, carriers, effect_per_copy, se, mlog10p_presence, mlog10p_heterogeneity, pve)]), row.names = FALSE)
print(as.data.frame(TB), row.names = FALSE)
say("all QTL together explain %.1f%% of the variance of %s (taxon means not counted). %.1f min.",
    100 * var(tot) / yv, trait, as.numeric(difftime(Sys.time(), t_start, units = "mins")))
say("WROTE %s/joint_qtl_%s{,_by_taxon,_genes,_steps,_perm}.tsv", out_dir, trait)
