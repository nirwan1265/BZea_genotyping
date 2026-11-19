#!/usr/bin/env Rscript
## =============================================================================
## 33 — Fine-mapping inside each QTL support interval with SuSiE (susieR; Wang et al. 2020).
##
## The support interval (1.5 drop in -log10 p) says where a QTL could be, but not how likely each window is.
## SuSiE ("sum of single effects") fits all windows of the region together and returns
##   - PIP per window: posterior probability that the window carries a causal effect;
##   - 95% credible sets (CS): smallest groups of windows that hold one causal effect with probability >= 0.95.
##     One CS per independent effect, so two QTL in one region give two CS.
## Model, the same as the scan (step 26) and the joint model (step 29):
##   y = taxon + other QTL of the trait (their top windows) + polygenic (LOCO kinship) + region windows + e
##   1. h2 of the polygenic part: REML on the model without the region (eigen-decomposition of K, 1-D search).
##   2. y, the region windows and the covariates are multiplied by V^(-1/2), V = h2 K + (1 - h2) I, which makes
##      the residuals independent; the covariates are then projected out of y and the windows.
##   3. susie(X, y, L = SUSIE_L): residual variance and prior effect variance are estimated by susieR.
## Windows inside one tract in every pool have identical dosage, so SuSiE cannot split them: they share the PIP
## and enter a CS together. The CS width is therefore the resolution the tract breakpoints allow.
##
## Usage:  Rscript 33_susie_finemap.R <assoc input dir (step 25)> <TRAIT> <qtl table (step 29 qtl2_<TRAIT>_qtl.tsv)>
##                                    <out dir> [genes table: qtl_genes/<TRAIT>_qtl_genes.tsv]
## Env:    PAD_MB (1: the region is the support interval +- PAD_MB), SUSIE_L (10), COVERAGE (0.95),
##         MIN_ABS_CORR (0.5: a CS whose windows are less correlated than this is dropped, susieR default),
##         FIX_COL (TAXON)
## Needs:  R with susieR and data.table.
## Outputs (<out dir>):
##   susie_<TRAIT>_cs.tsv     one row per credible set: QTL, windows, span, purity, lead window, genes in the CS
##   susie_<TRAIT>_pip.tsv    PIP per window and QTL region (with the CS it belongs to)
##   susie_<TRAIT>_qtl.tsv    one row per QTL: support interval width, number of CS, total CS width
## =============================================================================
suppressPackageStartupMessages({ library(susieR); library(data.table) })
args <- commandArgs(TRUE)
if (length(args) < 4) stop("usage: Rscript 33_susie_finemap.R <assoc dir> <TRAIT> <qtl table> <out dir> [genes table]")
in_dir <- args[1]; trait <- args[2]; qtl_f <- args[3]; out_dir <- args[4]; genes_f <- if (length(args) >= 5) args[5] else ""
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
envn <- function(k, d) { v <- Sys.getenv(k, ""); if (v == "") d else as.numeric(v) }
PAD <- envn("PAD_MB", 1) * 1e6; L <- envn("SUSIE_L", 10); COV <- envn("COVERAGE", 0.95); MINR <- envn("MIN_ABS_CORR", 0.5)
FIX_COL <- Sys.getenv("FIX_COL", "TAXON")

ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE)
g <- fread(file.path(in_dir, "geno_dosage.csv"), data.table = FALSE)
D <- as.matrix(g[, -1]); storage.mode(D) <- "double"; rownames(D) <- g[[1]]; rm(g)
map <- fread(file.path(in_dir, "window_map.tsv"), data.table = FALSE)
ph <- ph[!is.na(as.numeric(ph[[trait]])) & ph[[1]] %in% rownames(D), ]
ids <- ph[[1]]; D <- D[ids, , drop = FALSE]; y0 <- as.numeric(ph[[trait]])
Q <- fread(qtl_f, data.table = FALSE)
cut <- -log10(0.05 / nrow(Q))                                  # Bonferroni over the trait's QTL, as in the README
Q$holds <- Q$mlog10p_joint > cut
GN <- if (nzchar(genes_f) && file.exists(genes_f)) fread(genes_f, data.table = FALSE) else NULL
message(sprintf("SuSiE | trait %s | plots %d | QTL %d (%d hold) | region = interval +- %g Mb | L = %d",
                trait, length(ids), nrow(Q), sum(Q$holds), PAD / 1e6, L))

Kcache <- list()
getK <- function(chr) {
  if (is.null(Kcache[[chr]])) {
    k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", as.integer(sub("chr", "", chr)))), data.table = FALSE)
    K <- as.matrix(k[, -1]); rownames(K) <- k[[1]]; colnames(K) <- k[[1]]
    e <- eigen(K[ids, ids], symmetric = TRUE); e$values <- pmax(e$values, 0)
    Kcache[[chr]] <<- e
  }
  Kcache[[chr]]
}
## REML h2 for y ~ C + polygenic(K) + e, on the eigenbasis of K
reml_h2 <- function(yt, Ct, d) {
  n <- length(yt); p <- ncol(Ct)
  f <- function(h2) {
    w <- h2 * d + (1 - h2); Cw <- Ct / sqrt(w); yw <- yt / sqrt(w)
    XtX <- crossprod(Cw); b <- solve(XtX, crossprod(Cw, yw)); r <- yw - Cw %*% b
    s2 <- sum(r^2) / (n - p)
    0.5 * (sum(log(w)) + (n - p) * log(s2) + as.numeric(determinant(XtX)$modulus))
  }
  optimize(f, c(0.001, 0.999))$minimum
}

cs_rows <- list(); pip_rows <- list(); q_rows <- list()
for (i in seq_len(nrow(Q))) {
  q <- Q[i, ]; chr <- q$chr; lead <- q$qtl
  lo <- q$interval_start_Mb * 1e6 - PAD; hi <- q$interval_end_Mb * 1e6 + PAD
  w <- map[map$chr == chr & map$end > lo & map$start < hi, ]
  X <- D[, w$Marker, drop = FALSE]; keep <- apply(X, 2, function(x) length(unique(x)) > 1)
  X <- X[, keep, drop = FALSE]; w <- w[keep, ]
  ## covariates: intercept, taxon, the other QTL's top windows (outside this region)
  oth <- setdiff(Q$qtl, lead); oth <- oth[!oth %in% w$Marker & oth %in% colnames(D)]
  C <- cbind(1, model.matrix(~ factor(ph[[FIX_COL]]))[, -1, drop = FALSE], D[, oth, drop = FALSE])
  C <- C[, qr(C)$pivot[seq_len(qr(C)$rank)], drop = FALSE]
  e <- getK(chr); U <- e$vectors
  yt <- crossprod(U, y0); Ct <- crossprod(U, C); Xt <- crossprod(U, X)
  h2 <- reml_h2(yt, Ct, e$values); s <- 1 / sqrt(h2 * e$values + (1 - h2))
  yw <- yt * s; Cw <- Ct * s; Xw <- Xt * s
  P <- qr(Cw); yr <- qr.resid(P, yw); Xr <- qr.resid(P, Xw)
  colnames(Xr) <- w$Marker
  fit <- susie(Xr, as.numeric(yr), L = L, coverage = COV, min_abs_corr = MINR, verbose = FALSE)
  pip <- fit$pip; cs <- fit$sets$cs; pur <- fit$sets$purity
  inCS <- rep(NA_character_, ncol(Xr))
  ncs <- length(cs); width_tot <- 0
  if (ncs > 0) for (j in seq_len(ncs)) {
    m <- cs[[j]]; inCS[m] <- names(cs)[j]
    ww <- w[m, ]; span <- c(min(ww$start), max(ww$end))
    width_tot <- width_tot + sum(ww$end - ww$start)
    top <- ww$Marker[which.max(pip[m])]
    genes <- if (!is.null(GN)) {
      gq <- GN[GN$chr == chr, ]
      hit <- vapply(seq_len(nrow(gq)), function(k) any(gq$start[k] < ww$end & gq$end[k] > ww$start), TRUE)
      gq <- unique(gq[hit, c("gene", "gene_name")])
      list(n = nrow(gq), ids = paste(gq$gene, collapse = ";"),
           names = paste(ifelse(is.na(gq$gene_name) | gq$gene_name == "", gq$gene, gq$gene_name), collapse = "; "))
    } else list(n = NA, ids = "", names = "")
    cs_rows[[length(cs_rows) + 1]] <- data.frame(
      trait = trait, qtl = lead, holds_joint_model = q$holds, cs = names(cs)[j], n_windows = length(m),
      cs_Mb = round(sum(ww$end - ww$start) / 1e6, 2), span_start_Mb = span[1] / 1e6, span_end_Mb = span[2] / 1e6,
      purity_min_abs_corr = round(pur[j, "min.abs.corr"], 3), cs_pip_sum = round(sum(pip[m]), 3),
      lead_window = top, lead_pip = round(max(pip[m]), 3), scan_top_window_in_cs = lead %in% w$Marker[m],
      n_genes = genes$n, genes = genes$ids, gene_names = genes$names)
  }
  pip_rows[[i]] <- data.frame(trait = trait, qtl = lead, Marker = w$Marker, chr = chr, start = w$start, end = w$end,
                              pip = signif(pip, 4), cs = inCS)
  q_rows[[i]] <- data.frame(trait = trait, qtl = lead, holds_joint_model = q$holds, mlog10p_joint = q$mlog10p_joint,
                            interval_Mb = round(q$interval_end_Mb - q$interval_start_Mb, 2),
                            region_windows = ncol(Xr), h2_polygenic = round(h2, 3), n_cs = ncs,
                            cs_Mb_total = round(width_tot / 1e6, 2))
  message(sprintf("  %-18s interval %5.1f Mb | %3d windows | h2 %.2f | %d CS, %.1f Mb", lead,
                  q$interval_end_Mb - q$interval_start_Mb, ncol(Xr), h2, ncs, width_tot / 1e6))
}
fwrite(rbindlist(q_rows), file.path(out_dir, sprintf("susie_%s_qtl.tsv", trait)), sep = "\t")
fwrite(rbindlist(pip_rows), file.path(out_dir, sprintf("susie_%s_pip.tsv", trait)), sep = "\t")
if (length(cs_rows)) fwrite(rbindlist(cs_rows), file.path(out_dir, sprintf("susie_%s_cs.tsv", trait)), sep = "\t")
message("WROTE ", out_dir, "/susie_", trait, "_{qtl,pip,cs}.tsv")
