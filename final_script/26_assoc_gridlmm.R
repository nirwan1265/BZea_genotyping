#!/usr/bin/env Rscript
## ============================================================
## 26 — Association mapping on the k-mer introgression dosage (GridLMM, LOCO kinship)
##
## Based on run_bzea_gridlmm_scan (scripts/Analyze_QTL_Founder2.R), with:
##   - genotype = step-25 window dosage (0 B73, 1 SEG, 2 DONOR) for all 10 chromosomes
##   - LOCO kinship from the same dosage (step 25)
##   - scopes: ALL (every plot, taxon as fixed effect) + one scan per group (GROUP_COL, default
##     TAXON; F1FAMILY also possible) with at least MIN_N plots; within a group no fixed effect
##   - additive scan:   y ~ [taxon] + (1|Line) + b*dosage            (LRT, 1 df); b = days per teosinte copy
##   - dominance scan:  y ~ [taxon] + (1|Line) + b*dosage + c*SEG    vs without c (LRT, 1 df)
##       plot means: B73 plot = 0, DONOR plot = 2a, SEG plot = 1/4*0 + 1/2*(a+d) + 1/4*2a = a + d/2,
##       so c = d/2: c > 0 teosinte allele (partly) dominant, c < 0 (partly) recessive.
##       Only where the scope has >= 5 SEG and >= 5 DONOR plots (otherwise SEG and dosage coincide).
##   - windows per scope: >= MIN_CARRIERS plots of the scope carry donor (SEG or DONOR)
##   - permutation threshold per scope: phenotypes shuffled within taxon (ALL) or within the group,
##     genome-wide max -log10 p per permutation; threshold = 95th percentile (mode = "perm", in chunks)
##
## Usage:
##   Rscript 26_assoc_gridlmm.R scan    <in_dir> <trait> <out_dir> [cores]
##   Rscript 26_assoc_gridlmm.R perm    <in_dir> <trait> <out_dir> <n_perm> <seed> [cores]
##   Rscript 26_assoc_gridlmm.R summary <in_dir> <trait> <out_dir>
## Environment (optional): GROUP_COL (TAXON | F1FAMILY | none; default TAXON), MIN_N (default 20),
##   MIN_CARRIERS (per-group scans, default 10)
## in_dir = step-25 output (pheno.csv, geno_dosage.csv, geno_seg.csv, window_map.tsv, K_matrix_chr*.txt)
## ============================================================
suppressPackageStartupMessages({
  library(GridLMM)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; in_dir <- args[2]; trait <- args[3]; out_dir <- args[4]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
GROUP_COL <- Sys.getenv("GROUP_COL", "TAXON")
MIN_N <- as.integer(Sys.getenv("MIN_N", "20"))
MIN_CARRIERS <- as.integer(Sys.getenv("MIN_CARRIERS", "10"))
cores <- as.integer(if (mode == "perm") (if (length(args) >= 7) args[7] else 1) else (if (length(args) >= 5) args[5] else 1))

## ---------- inputs ----------
ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE)
dat <- data.frame(Line = ph[[1]], y = as.numeric(ph[[trait]]), TAXON = factor(ph$TAXON), stringsAsFactors = FALSE)
if (GROUP_COL != "none") dat$GROUP <- as.character(ph[[GROUP_COL]])
read_geno <- function(f) {
  g <- fread(f, data.table = FALSE)
  X <- as.matrix(g[, -1, drop = FALSE]); storage.mode(X) <- "double"; rownames(X) <- g[[1]]; X
}
D <- read_geno(file.path(in_dir, "geno_dosage.csv"))
S <- read_geno(file.path(in_dir, "geno_seg.csv"))
map <- fread(file.path(in_dir, "window_map.tsv"), data.table = FALSE)
dat <- dat[!is.na(dat$y) & dat$Line %in% rownames(D), ]
D <- D[dat$Line, , drop = FALSE]; S <- S[dat$Line, , drop = FALSE]
chrs <- unique(map$chr)
Kall <- lapply(chrs, function(chr) {
  k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", as.integer(sub("chr", "", chr)))), data.table = FALSE)
  K <- as.matrix(k[, -1]); rownames(K) <- k[[1]]; colnames(K) <- k[[1]]; K
}); names(Kall) <- chrs

scopes <- "ALL"
if (GROUP_COL != "none") {
  n_g <- table(dat$GROUP); scopes <- c(scopes, sort(names(n_g)[n_g >= MIN_N]))
}
scope_rows <- function(sc) if (sc == "ALL") seq_len(nrow(dat)) else which(dat$GROUP == sc)
safe <- function(sc) gsub("[^A-Za-z0-9._-]", "_", sc)
message(sprintf("%s | trait %s | plots %d | windows %d | scopes: %s", mode, trait, nrow(dat), nrow(map),
                paste(sapply(scopes, function(sc) sprintf("%s(%d)", sc, length(scope_rows(sc)))), collapse = " ")))

## windows tested in a scope: enough carriers, and the dosage varies
scope_windows <- function(rows, sc) {
  Dk <- D[rows, , drop = FALSE]
  car <- colSums(Dk > 0); v <- apply(Dk, 2, function(x) length(unique(x)) > 1)
  keep <- v & car >= (if (sc == "ALL") 1 else MIN_CARRIERS)
  dom_ok <- colSums(S[rows, , drop = FALSE]) >= 5 & colSums(Dk == 2) >= 5
  list(keep = keep, dom_ok = dom_ok, car = car)
}

## one scope, one chromosome: additive (+ dominance) scan, reusing the kinship decomposition
scan_chr <- function(d, rows, chr, keep, dom_ok, dominance, V_setup = NULL) {
  mk <- map$Marker[map$chr == chr & keep[map$Marker]]
  if (length(mk) < 1) return(NULL)
  fml <- if (length(unique(d$TAXON)) > 1) y ~ TAXON + (1 | Line) else y ~ 1 + (1 | Line)
  K <- Kall[[chr]][d$Line, d$Line]
  add <- GridLMM_GWAS(formula = fml, test_formula = ~1, reduced_formula = ~0, data = d,
                      X = D[d$Line, mk, drop = FALSE], X_ID = "Line", relmat = list(Line = K),
                      method = "ML", centerX = FALSE, scaleX = FALSE, fillNAX = FALSE,
                      V_setup = V_setup, mc.cores = cores, verbose = FALSE)
  nb <- ncol(add$setup$X_cov)                          # fixed covariates; the marker effect is beta.(nb+1)
  res <- data.frame(Marker = mk, p_add = add$results$p_value_ML, beta_add = add$results[[paste0("beta.", nb + 1)]],
                    h2 = add$results$Line.ML, p_dom = NA_real_, c_dom = NA_real_)
  md <- mk[dom_ok[mk]]
  if (dominance && length(md) > 0) {
    st <- add$setup
    Xd <- D[d$Line, md, drop = FALSE]; Xs <- S[d$Line, md, drop = FALSE]
    dom <- GridLMM:::run_GridLMM_GWAS(st$Y, st$X_cov, list(Xd, Xs), list(Xd), rep(0, nb + 2), NULL,
                                      st$V_setup, st$h2_start, 0.01, 0.99, st$proximal_markers, st$proximal_Xs, "ML", cores, FALSE)
    i <- match(md, res$Marker)
    res$p_dom[i] <- dom$p_value_ML
    res$c_dom[i] <- dom[[paste0("beta.", nb + 2)]]     # c = d/2 (effect of SEG beyond the additive dosage)
  }
  list(res = res, V_setup = add$setup$V_setup)
}

if (mode == "scan") {
  for (sc in scopes) {
    t0 <- Sys.time(); rows <- scope_rows(sc); d <- dat[rows, ]
    w <- scope_windows(rows, sc); out <- list()
    for (chr in chrs) out[[chr]] <- scan_chr(d, rows, chr, w$keep, w$dom_ok, dominance = TRUE)$res
    res <- merge(map[, c("Marker", "chr", "start", "end", "pos")], rbindlist(out), by = "Marker", sort = FALSE)
    res$n_plots <- nrow(d); res$n_carriers <- w$car[res$Marker]
    res$n_seg <- colSums(S[rows, res$Marker, drop = FALSE]); res$n_donor <- colSums(D[rows, res$Marker, drop = FALSE] == 2)
    res$mlog10p_add <- -log10(res$p_add); res$mlog10p_dom <- -log10(res$p_dom)
    f <- file.path(out_dir, if (sc == "ALL") sprintf("scan_%s.tsv", trait) else sprintf("scan_%s_%s.tsv", trait, safe(sc)))
    fwrite(res, f, sep = "\t")
    message(sprintf("%s: %d plots, %d windows, %.1f min -> %s", sc, nrow(d), nrow(res),
                    as.numeric(difftime(Sys.time(), t0, units = "mins")), f))
  }
}

if (mode == "perm") {
  n_perm <- as.integer(args[5]); seed <- as.integer(args[6]); set.seed(seed)
  W <- lapply(scopes, function(sc) scope_windows(scope_rows(sc), sc)); names(W) <- scopes
  Vs <- list(); out <- list()
  for (b in seq_len(n_perm)) {
    for (sc in scopes) {
      rows <- scope_rows(sc); d <- dat[rows, ]
      blk <- if (sc == "ALL") d$TAXON else rep(1, nrow(d))           # shuffle within taxon (ALL) or within the group
      for (g in unique(blk)) { i <- which(blk == g); d$y[i] <- d$y[i][sample.int(length(i))] }
      m <- 0
      for (chr in chrs) {
        key <- paste(sc, chr)
        r <- scan_chr(d, rows, chr, W[[sc]]$keep, W[[sc]]$dom_ok, dominance = FALSE, V_setup = Vs[[key]])
        if (is.null(r)) next
        Vs[[key]] <- r$V_setup                          # the kinship decomposition doesn't depend on y: reuse it
        m <- max(m, -log10(min(r$res$p_add)))
      }
      out[[length(out) + 1]] <- data.frame(seed = seed, perm = b, scope = sc, max_mlog10p_add = m)
    }
  }
  fwrite(rbindlist(out), file.path(out_dir, sprintf("perm_%s_seed%d.tsv", trait, seed)), sep = "\t")
}

if (mode == "summary") {
  pf <- list.files(out_dir, sprintf("^perm_%s_seed.*\\.tsv$", trait), full.names = TRUE)
  P <- if (length(pf)) rbindlist(lapply(pf, fread)) else NULL
  summ <- list()
  for (sc in scopes) {
    f <- file.path(out_dir, if (sc == "ALL") sprintf("scan_%s.tsv", trait) else sprintf("scan_%s_%s.tsv", trait, safe(sc)))
    if (!file.exists(f)) next
    res <- fread(f, data.table = FALSE)
    pp <- if (!is.null(P)) P$max_mlog10p_add[P$scope == sc] else numeric(0)
    thr <- if (length(pp) >= 20) unname(quantile(pp, 0.95)) else NA
    top <- res[which.min(res$p_add), ]
    summ[[sc]] <- data.frame(trait = trait, scope = sc, plots = top$n_plots, windows = nrow(res), permutations = length(pp),
                             threshold = round(thr, 2), top_window = top$Marker, top_mlog10p = round(top$mlog10p_add, 2),
                             top_beta = round(top$beta_add, 2), n_above = if (is.finite(thr)) sum(res$mlog10p_add > thr) else NA)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      suppressPackageStartupMessages(library(ggplot2))
      res$chr_n <- as.integer(sub("chr", "", res$chr)); res <- res[order(res$chr_n, res$pos), ]
      cl <- tapply(map$end, as.integer(sub("chr", "", map$chr)), max)      # same x axis for every scope
      off <- c(0, cumsum(cl))[seq_along(cl)]; names(off) <- names(cl)
      res$x <- res$pos + off[as.character(res$chr_n)]
      mids <- off + cl / 2
      for (what in c("add", "dom")) {
        if (all(is.na(res[[paste0("mlog10p_", what)]]))) next
        p <- ggplot(res[!is.na(res[[paste0("mlog10p_", what)]]), ],
                    aes(x = x / 1e6, y = .data[[paste0("mlog10p_", what)]], colour = factor(chr_n %% 2))) +
          geom_point(size = 0.6, show.legend = FALSE) + scale_colour_manual(values = c("#2a78d6", "#0b3a73")) +
          scale_x_continuous(breaks = mids / 1e6, labels = names(mids), limits = c(0, sum(cl) / 1e6)) +
          labs(x = "Chromosome", y = expression(-log[10](p)),
               title = sprintf("%s, %s effect, %s (%d plots)", trait, if (what == "add") "additive" else "dominance",
                               if (sc == "ALL") "all taxa" else sc, top$n_plots)) +
          theme_classic(base_size = 11)
        if (what == "add" && is.finite(thr)) p <- p + geom_hline(yintercept = thr, colour = "#0b0b0b") +
            annotate("text", x = 0, y = thr, label = "5% permutation", hjust = 0, vjust = -0.4, size = 3)
        ggsave(file.path(out_dir, sprintf("manhattan_%s_%s_%s.png", trait, safe(sc), what)), p, width = 10, height = 3.6, dpi = 200)
      }
    }
  }
  s <- rbindlist(summ); print(s)
  fwrite(s, file.path(out_dir, sprintf("summary_%s.tsv", trait)), sep = "\t")
}
