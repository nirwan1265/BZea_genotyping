#!/usr/bin/env Rscript
## ============================================================
## 26 — Association mapping on the k-mer introgression dosage (GridLMM, LOCO kinship)
##
## Based on run_bzea_gridlmm_scan (scripts/Analyze_QTL_Founder2.R), with:
##   - genotype = step-25 window dosage (0 B73, 1 SEG, 2 DONOR) for all 10 chromosomes
##   - LOCO kinship from the same dosage (step 25)
##   - taxon (Zd/Zl/Zv/Zx) as fixed effect; "per-taxon" scans optional
##   - additive scan:   y ~ taxon + (1|Line) + b*dosage            (LRT, 1 df)
##   - dominance scan:  y ~ taxon + (1|Line) + b*dosage + c*SEG    vs without c (LRT, 1 df)
##       plot means: B73 plot = 0, DONOR plot = 2a, SEG plot = 1/4*0 + 1/2*(a+d) + 1/4*2a = a + d/2,
##       so c = d/2: c > 0 teosinte allele (partly) dominant, c < 0 (partly) recessive
##   - permutation threshold: phenotypes shuffled within taxon, genome-wide max -log10 p per
##     permutation; threshold = 95th percentile (run in chunks: mode = "perm")
##
## Usage:
##   Rscript 26_assoc_gridlmm.R scan <in_dir> <trait> <out_dir> [cores]
##   Rscript 26_assoc_gridlmm.R perm <in_dir> <trait> <out_dir> <n_perm> <seed> [cores]
##   Rscript 26_assoc_gridlmm.R summary <in_dir> <trait> <out_dir>
## in_dir = step-25 output (pheno.csv, geno_dosage.csv, geno_seg.csv, window_map.tsv, K_matrix_chr*.txt)
## ============================================================
suppressPackageStartupMessages({
  library(GridLMM)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; in_dir <- args[2]; trait <- args[3]; out_dir <- args[4]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load_inputs <- function(in_dir, trait) {
  ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE)
  dat <- data.frame(Line = ph[[1]], y = as.numeric(ph[[trait]]), TAXON = factor(ph$TAXON), stringsAsFactors = FALSE)
  dat <- dat[!is.na(dat$y), ]
  read_geno <- function(f) {
    g <- fread(f, data.table = FALSE)
    X <- as.matrix(g[, -1, drop = FALSE]); storage.mode(X) <- "double"; rownames(X) <- g[[1]]; X
  }
  D <- read_geno(file.path(in_dir, "geno_dosage.csv"))
  S <- read_geno(file.path(in_dir, "geno_seg.csv"))
  map <- fread(file.path(in_dir, "window_map.tsv"), data.table = FALSE)
  dat <- dat[dat$Line %in% rownames(D), ]
  list(dat = dat, D = D[dat$Line, , drop = FALSE], S = S[dat$Line, , drop = FALSE], map = map)
}

read_K <- function(in_dir, chr_num, ids) {
  k <- fread(file.path(in_dir, sprintf("K_matrix_chr%d.txt", chr_num)), data.table = FALSE)
  K <- as.matrix(k[, -1]); rownames(K) <- k[[1]]; colnames(K) <- k[[1]]
  K[ids, ids]
}

## one chromosome: additive scan (+ optional dominance scan), reusing the V decomposition
scan_chr <- function(dat, D, S, map, chr, K, cores, dominance = TRUE, V_setup = NULL) {
  mk <- map$Marker[map$chr == chr]
  fml <- y ~ TAXON + (1 | Line)
  add <- GridLMM_GWAS(formula = fml, test_formula = ~1, reduced_formula = ~0, data = dat,
                      X = D[, mk, drop = FALSE], X_ID = "Line", relmat = list(Line = K),
                      method = "ML", centerX = FALSE, scaleX = FALSE, fillNAX = FALSE,
                      V_setup = V_setup, mc.cores = cores, verbose = FALSE)
  nb <- ncol(add$setup$X_cov)                         # fixed covariates; the marker effect is beta.(nb+1)
  res <- data.frame(Marker = mk, p_add = add$results$p_value_ML, beta_add = add$results[[paste0("beta.", nb + 1)]],
                    h2 = add$results$Line.ML)
  if (dominance) {
    st <- add$setup
    Xd <- D[match(dat$Line, rownames(D)), mk, drop = FALSE]; Xs <- S[match(dat$Line, rownames(S)), mk, drop = FALSE]
    dom <- GridLMM:::run_GridLMM_GWAS(st$Y, st$X_cov, list(Xd, Xs), list(Xd), rep(0, ncol(st$X_cov) + 2), NULL,
                                      st$V_setup, st$h2_start, 0.01, 0.99, st$proximal_markers, st$proximal_Xs, "ML", cores, FALSE)
    res$p_dom <- dom$p_value_ML
    res$c_dom <- dom[[paste0("beta.", nb + 2)]]      # c = d/2 (effect of SEG beyond the additive dosage)
  }
  list(res = res, V_setup = add$setup$V_setup)
}

inp <- load_inputs(in_dir, trait)
dat <- inp$dat; map <- inp$map
chrs <- unique(map$chr)
cores <- as.integer(if (mode == "perm") (if (length(args) >= 7) args[7] else 1) else (if (length(args) >= 5) args[5] else 1))
message(sprintf("%s | trait %s | plots %d | windows %d | taxa %s", mode, trait, nrow(dat), nrow(map),
                paste(names(table(dat$TAXON)), table(dat$TAXON), sep = ":", collapse = " ")))

if (mode == "scan") {
  t0 <- Sys.time(); out <- list()
  for (chr in chrs) {
    K <- read_K(in_dir, as.integer(sub("chr", "", chr)), dat$Line)
    out[[chr]] <- scan_chr(dat, inp$D, inp$S, map, chr, K, cores, dominance = TRUE)$res
    message(chr, " done ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
  res <- merge(map, rbindlist(out), by = "Marker", sort = FALSE)
  res$mlog10p_add <- -log10(res$p_add); res$mlog10p_dom <- -log10(res$p_dom)
  fwrite(res, file.path(out_dir, sprintf("scan_%s.tsv", trait)), sep = "\t")
  message("WROTE ", file.path(out_dir, sprintf("scan_%s.tsv", trait)))
}

if (mode == "perm") {
  n_perm <- as.integer(args[5]); seed <- as.integer(args[6]); set.seed(seed)
  Ks <- lapply(chrs, function(chr) read_K(in_dir, as.integer(sub("chr", "", chr)), dat$Line)); names(Ks) <- chrs
  Vs <- list(); maxp <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    d2 <- dat
    for (tx in levels(d2$TAXON)) { i <- which(d2$TAXON == tx); d2$y[i] <- d2$y[sample(i)] }   # shuffle within taxon
    m <- 0
    for (chr in chrs) {
      r <- scan_chr(d2, inp$D, inp$S, map, chr, Ks[[chr]], cores, dominance = FALSE, V_setup = Vs[[chr]])
      Vs[[chr]] <- r$V_setup                       # the kinship decomposition doesn't depend on y: reuse it
      m <- max(m, -log10(min(r$res$p_add)))
    }
    maxp[b] <- m
  }
  fwrite(data.frame(seed = seed, perm = seq_len(n_perm), max_mlog10p_add = maxp),
         file.path(out_dir, sprintf("perm_%s_seed%d.tsv", trait, seed)), sep = "\t")
}

if (mode == "summary") {
  res <- fread(file.path(out_dir, sprintf("scan_%s.tsv", trait)), data.table = FALSE)
  pf <- list.files(out_dir, sprintf("^perm_%s_seed.*\\.tsv$", trait), full.names = TRUE)
  thr <- if (length(pf)) quantile(rbindlist(lapply(pf, fread))$max_mlog10p_add, 0.95) else NA
  n_eff_bonf <- -log10(0.05 / nrow(res))
  message(sprintf("permutations: %d | 5%% genome-wide threshold -log10p = %s | Bonferroni (all windows) = %.2f",
                  if (length(pf)) nrow(rbindlist(lapply(pf, fread))) else 0, format(thr, digits = 3), n_eff_bonf))
  res <- res[order(res$p_add), ]
  print(head(res[, c("Marker", "chr", "pos", "n_carriers", "beta_add", "mlog10p_add", "c_dom", "mlog10p_dom")], 15))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    res$chr_n <- as.integer(sub("chr", "", res$chr)); res <- res[order(res$chr_n, res$pos), ]
    off <- cumsum(c(0, tapply(res$pos, res$chr_n, max)))[seq_along(unique(res$chr_n))]
    res$x <- res$pos + off[match(res$chr_n, sort(unique(res$chr_n)))]
    mids <- tapply(res$x, res$chr_n, function(v) mean(range(v)))
    for (what in c("add", "dom")) {
      p <- ggplot(res, aes(x = x / 1e6, y = .data[[paste0("mlog10p_", what)]], colour = factor(chr_n %% 2))) +
        geom_point(size = 0.6, show.legend = FALSE) + scale_colour_manual(values = c("#2a78d6", "#0b3a73")) +
        scale_x_continuous(breaks = mids / 1e6, labels = names(mids)) +
        labs(x = "Chromosome", y = expression(-log[10](p)),
             title = sprintf("%s, %s effect (k-mer introgression dosage, GridLMM, LOCO K)", trait,
                             if (what == "add") "additive" else "dominance")) +
        theme_classic(base_size = 11)
      if (what == "add" && is.finite(thr)) p <- p + geom_hline(yintercept = thr, colour = "#0b0b0b") +
          annotate("text", x = 0, y = thr, label = "5% permutation", hjust = 0, vjust = -0.4, size = 3)
      ggsave(file.path(out_dir, sprintf("manhattan_%s_%s.png", trait, what)), p, width = 10, height = 3.6, dpi = 200)
    }
  }
}
