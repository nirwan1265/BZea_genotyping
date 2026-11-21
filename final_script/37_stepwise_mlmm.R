#!/usr/bin/env Rscript
## =============================================================================
## 37 — Stepwise multi-locus mixed model (MLMM-style; Segura et al. 2012) on the k-mer dosage: a check of the
##      QTL of steps 26-29 with a search that adds AND removes QTL while keeping the kinship term.
##
## Model: y ~ taxon + QTL dosages (fixed) + polygenic background (random, LOCO kinship of the tested chromosome),
##        fitted with GridLMM (ML), as in step 26. Genotypes, plots and kinship: step 25.
## 1. Forward:  scan the genome with the current QTL as covariates; if the best window passes the genome-wide
##              threshold (the step-26 permutation threshold, -log10 p), add it; repeat until none passes.
## 2. Backward: test each QTL given all the others; if the weakest is below the same threshold, drop it, refit,
##              repeat until every QTL passes.
## 3. Refine:   move each QTL to the best window of its chromosome given all the other QTL (one pass).
## The same genome-wide threshold is used for adding and keeping a QTL, so the final model has no QTL that the
## others explain. This is stricter than the Bonferroni cut of the joint test in step 29.
##
## Usage (from final_script/, step-25 inputs and the step-26 summary):
##   Rscript 37_stepwise_mlmm.R <step-25 dir> <TRAIT> <step-26 summary_<TRAIT>.tsv> <out dir> [cores]
## e.g. Rscript 37_stepwise_mlmm.R $HOME/bzea_assoc_fix1 DTS ../data/introgression_fix1/assoc_gridlmm_spats/summary_DTS.tsv \
##        ../data/introgression_fix1/stepwise 2
## Outputs: stepwise_<TRAIT>_steps.tsv (every add/drop/move), stepwise_<TRAIT>_qtl.tsv (final model: window,
##          -log10 p given the others, effect per teosinte copy, SE)
## Needs: R with GridLMM and data.table.
## =============================================================================
suppressPackageStartupMessages({ library(GridLMM); library(data.table) })
args <- commandArgs(TRUE)
if (length(args) < 4) stop("usage: Rscript 37_stepwise_mlmm.R <step-25 dir> <TRAIT> <summary tsv> <out dir> [cores]")
in_dir <- args[1]; trait <- args[2]; thr <- fread(args[3])$threshold[1]; out <- args[4]
cores <- if (length(args) >= 5) as.integer(args[5]) else 1
dir.create(out, showWarnings = FALSE, recursive = TRUE)

ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE)
dat <- data.frame(Line = ph[[1]], y = as.numeric(ph[[trait]]), TAXON = factor(ph$TAXON), stringsAsFactors = FALSE)
g <- fread(file.path(in_dir, "geno_dosage.csv"), data.table = FALSE)
D <- as.matrix(g[, -1]); storage.mode(D) <- "double"; rownames(D) <- g[[1]]
map <- fread(file.path(in_dir, "window_map.tsv"), data.table = FALSE)
dat <- dat[!is.na(dat$y) & dat$Line %in% rownames(D), ]
D <- D[dat$Line, , drop = FALSE]
D <- D[, apply(D, 2, function(x) length(unique(x)) > 1), drop = FALSE]
map <- map[map$Marker %in% colnames(D), ]
chrs <- unique(map$chr)
K <- lapply(chrs, function(chr) {
  k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", as.integer(sub("chr", "", chr)))), data.table = FALSE)
  m <- as.matrix(k[, -1]); rownames(m) <- colnames(m) <- k[[1]]; m[dat$Line, dat$Line]
}); names(K) <- chrs
chr_of <- setNames(map$chr, map$Marker)
message(sprintf("%s: %d plots, %d windows, threshold -log10 p %.2f", trait, nrow(dat), ncol(D), thr))

## scan the windows `mk` of one chromosome with the QTL `qtl` as covariates (the QTL windows themselves are not
## retested: on a QTL's own chromosome every other window is tested given that QTL)
scan_chr <- function(chr, qtl, mk = map$Marker[map$chr == chr]) {
  d <- dat; cov <- qtl; mk <- setdiff(mk, qtl)
  if (!length(mk)) return(NULL)
  for (q in cov) d[[make.names(q)]] <- D[, q]
  fml <- as.formula(paste("y ~ TAXON", paste(c("", make.names(cov)), collapse = " + "), "+ (1 | Line)"))
  r <- GridLMM_GWAS(formula = fml, test_formula = ~1, reduced_formula = ~0, data = d, X = D[, mk, drop = FALSE],
                    X_ID = "Line", relmat = list(Line = K[[chr]]), method = "ML", centerX = FALSE, scaleX = FALSE,
                    fillNAX = FALSE, mc.cores = cores, verbose = FALSE)
  nb <- ncol(r$setup$X_cov)
  data.frame(Marker = mk, chr = chr, mlog10p = -log10(r$results$p_value_ML), beta = r$results[[paste0("beta.", nb + 1)]])
}
scan_all <- function(qtl) do.call(rbind, lapply(chrs, scan_chr, qtl = qtl))
## each QTL given all the others: its own window, scanned with the others as covariates
test_each <- function(qtl) do.call(rbind, lapply(qtl, function(q) scan_chr(chr_of[[q]], setdiff(qtl, q), q)))

steps <- list(); qtl <- character(0); s <- 0
repeat {                                   # 1. forward
  sc <- scan_all(qtl); b <- sc[which.max(sc$mlog10p), ]
  s <- s + 1; message(sprintf("forward %d: best %s, -log10 p %.2f", s, b$Marker, b$mlog10p))
  if (b$mlog10p < thr) break
  qtl <- c(qtl, b$Marker)
  steps[[length(steps) + 1]] <- data.frame(step = s, action = "add", Marker = b$Marker, mlog10p = b$mlog10p, beta = b$beta)
  if (length(qtl) >= 20) break
}
repeat {                                   # 2. backward
  if (!length(qtl)) break
  te <- test_each(qtl); w <- te[which.min(te$mlog10p), ]
  if (w$mlog10p >= thr) break
  message(sprintf("backward: drop %s (-log10 p %.2f given the others)", w$Marker, w$mlog10p))
  qtl <- setdiff(qtl, w$Marker)
  steps[[length(steps) + 1]] <- data.frame(step = NA, action = "drop", Marker = w$Marker, mlog10p = w$mlog10p, beta = w$beta)
}
for (i in seq_along(qtl)) {                # 3. refine
  q <- qtl[i]; sc <- scan_chr(chr_of[[q]], setdiff(qtl, q)); b <- sc[which.max(sc$mlog10p), ]
  if (b$Marker != q) {
    message(sprintf("refine: %s -> %s", q, b$Marker))
    steps[[length(steps) + 1]] <- data.frame(step = NA, action = paste("move from", q), Marker = b$Marker, mlog10p = b$mlog10p, beta = b$beta)
    qtl[i] <- b$Marker; chr_of[[b$Marker]] <- b$chr
  }
}
fin <- if (length(qtl)) test_each(qtl) else data.frame()
if (nrow(fin)) {
  fin$pos_Mb <- map$pos[match(fin$Marker, map$Marker)] / 1e6
  fin <- fin[order(as.integer(sub("chr", "", fin$chr)), fin$pos_Mb), c("chr", "pos_Mb", "Marker", "mlog10p", "beta")]
  names(fin)[4:5] <- c("mlog10p_given_others", "effect_per_copy")
}
fwrite(do.call(rbind, steps), file.path(out, sprintf("stepwise_%s_steps.tsv", trait)), sep = "\t")
fwrite(fin, file.path(out, sprintf("stepwise_%s_qtl.tsv", trait)), sep = "\t")
print(fin)
message("WROTE ", file.path(out, sprintf("stepwise_%s_qtl.tsv", trait)))
