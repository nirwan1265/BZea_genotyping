#!/usr/bin/env Rscript
## =============================================================================
## 29 — QTL scan, joint QTL model and per-taxon QTL effects with R/qtl2 (Broman et al. 2019).
##
## All statistics are R/qtl2 functions; this script only puts the k-mer calls into R/qtl2's format.
##   Genotypes: teosinte allele probability per 100-kb window = dosage / 2 (B73 0, SEG 0.5, DONOR 1),
##              as R/qtl2 allele probabilities with two alleles, B (B73) and T (teosinte).
##   Kinship:   calc_kinship(probs, "loco")            leave-one-chromosome-out
##   Scan:      scan1(probs, y, kinship, addcovar = taxon)          additive, 1 df
##   Threshold: scan1perm(..., perm_strata = taxon)                 5% genome-wide, permutations within taxon
##   Peaks:     find_peaks(threshold, peakdrop, drop = 1.5)         QTL and 1.5-LOD support intervals
##   Taxa differ? scan1(..., intcovar = taxon)                      LOD(per-taxon effects) - LOD(one effect)
##   Joint:     fit1 at each QTL with taxon + the other QTL as covariates: conditional LOD, effect per copy, SE
##   Per taxon: fit1 at each QTL within each taxon (other QTL as covariates): effect per copy, SE
##   Effect per teosinte copy: fit1(contrasts = cbind(mu = c(1, 1), a = c(-1, 1))), coefficient a.
##   More QTL in a region? Conditional scans (forward selection): scan1 with taxon + all QTL so far as covariates;
##              the best window above the threshold is added as a QTL, and the scan is repeated until none is left.
##              Only windows that are significant in the ordinary scan and not in LD with a QTL already in the
##              model (r^2 <= MAX_R2_WITH_QTL) can be added, so this splits significant regions, it does not
##              create new ones. This is how two peaks inside one long significant region are told apart.
##
## Usage:  Rscript 29_joint_qtl.R <assoc input dir (step 25)> <TRAIT> <out dir>
## Without permutations (NPERM=0), the QTL come from:
##   PEAKS=<step-27 peaks_<TRAIT>_ALL.tsv>  the k-mer GridLMM scan's peaks: the same top windows and support intervals,
##                                          so both analyses describe the same QTL (29.2_joint_qtl_noperm.sh); or
##   THRESHOLD_LOD=<LOD>                    peaks of the R/qtl2 additive scan above this LOD
## Env:    TAXON_COL (TAXON | SUBTAXON), NPERM (1000), CORES (1), PEAKDROP (Inf: one QTL per chromosome per scan; e.g. 2 allows more, but tracts
##         are long and neighbouring peaks on a chromosome are rarely separable), DROP (1.5),
##         MIN_TAXON_CARRIERS (5), GFF (B73 v5 GFF3: genes in the intervals),
##         COND_MLOG10P (genome-wide -log10 p threshold for the conditional scans; default: the permutation threshold;
##         29.2 passes the step-26 one), MAX_ADDED_QTL (10; 0 = no conditional scans), MAX_R2_WITH_QTL (0.5),
##         GRIDLMM_SCAN (step-26 scan_<TRAIT>.tsv: GridLMM -log10 p and effect for added QTL)
## Needs:  R with qtl2 and data.table (both on CRAN).
## Outputs (<out dir>):
##   qtl2_<TRAIT>_scan.tsv          LOD per window: additive, per-taxon (intcovar), difference
##   qtl2_<TRAIT>_perm.tsv          permutation maxima (additive and per-taxon models)
##   qtl2_<TRAIT>_qtl.tsv           one row per QTL: interval, LOD, conditional LOD (joint), effect ± SE, LOD taxa-differ, PVE
##   qtl2_<TRAIT>_by_taxon.tsv      QTL x taxon: carriers (SEG/DONOR), effect per copy ± SE, LOD
##   qtl2_<TRAIT>_genes.tsv         B73 v5 genes in each QTL interval (with GFF)
##   qtl2_<TRAIT>_added.tsv         QTL added by the conditional scans: round, window, conditional LOD and -log10 p
##   qtl2_<TRAIT>_conditional.tsv   the last conditional scan (LOD per window given all QTL)
##   qtl2_<TRAIT>_manhattan.png     additive and per-taxon LOD, thresholds
## =============================================================================
suppressPackageStartupMessages({ library(qtl2); library(data.table) })
args <- commandArgs(TRUE)
if (length(args) < 3) stop("usage: Rscript 29_joint_qtl.R <assoc input dir> <TRAIT> <out dir>")
in_dir <- args[1]; trait <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
envn <- function(k, d) { v <- Sys.getenv(k, ""); if (v == "") d else as.numeric(v) }
NPERM <- envn("NPERM", 1000); CORES <- envn("CORES", 1); PEAKDROP <- envn("PEAKDROP", Inf); DROP <- envn("DROP", 1.5)
MIN_TC <- envn("MIN_TAXON_CARRIERS", 5); GFF <- Sys.getenv("GFF", "")
PEAKS <- Sys.getenv("PEAKS", ""); THRESHOLD_LOD <- Sys.getenv("THRESHOLD_LOD", "")
say <- function(...) message(sprintf(...)); t0 <- Sys.time()
lod2p <- function(lod, df) -pchisq(2 * log(10) * lod, df, lower.tail = FALSE, log.p = TRUE) / log(10)   # -log10 p

## ---------- k-mer calls -> R/qtl2 allele probabilities ----------
ph <- fread(file.path(in_dir, "pheno.csv"), data.table = FALSE)
if (!trait %in% names(ph)) stop("trait ", trait, " not in pheno.csv")
g <- fread(file.path(in_dir, "geno_dosage.csv"), data.table = FALSE)
G <- as.matrix(g[, -1]); storage.mode(G) <- "double"; rownames(G) <- g[[1]]; rm(g)
map <- fread(file.path(in_dir, "window_map.tsv"), data.table = FALSE)
ph <- ph[!is.na(as.numeric(ph[[trait]])) & ph[[1]] %in% rownames(G), ]
ids <- ph[[1]]; G <- G[ids, , drop = FALSE]
y <- matrix(as.numeric(ph[[trait]]), ncol = 1, dimnames = list(ids, trait))
TAXON_COL <- Sys.getenv("TAXON_COL", "TAXON")   # TAXON, or SUBTAXON (Zx split by race): covariate and per-group effects
taxon <- setNames(ph[[TAXON_COL]], ids); taxa <- sort(unique(taxon))
covar <- model.matrix(~ factor(taxon))[, -1, drop = FALSE]; rownames(covar) <- ids     # taxon (no intercept)
intc <- covar                                                                            # per-taxon QTL effects
chrs <- unique(map$chr); cn <- sub("^chr", "", chrs)
probs <- setNames(lapply(chrs, function(ch) {
  m <- map$Marker[map$chr == ch]; p <- G[, m, drop = FALSE] / 2
  a <- array(NA_real_, c(length(ids), 2, length(m)), dimnames = list(ids, c("B", "T"), m))
  a[, "B", ] <- 1 - p; a[, "T", ] <- p; a
}), cn)
attr(probs, "crosstype") <- "risib"; attr(probs, "is_x_chr") <- setNames(rep(FALSE, length(cn)), cn)
attr(probs, "alleles") <- c("B", "T"); attr(probs, "alleleprobs") <- TRUE
class(probs) <- c("calc_genoprob", "list")
pmap <- setNames(lapply(chrs, function(ch) { x <- map[map$chr == ch, ]; setNames(x$pos / 1e6, x$Marker) }), cn)
say("qtl2 | trait %s | plots %d | windows %d | taxa %s", trait, length(ids), ncol(G),
    paste(sprintf("%s %d", taxa, as.integer(table(taxon)[taxa])), collapse = ", "))

## ---------- kinship, scans, permutations ----------
K <- calc_kinship(probs, "loco", cores = CORES)
add <- scan1(probs, y, K, addcovar = covar, cores = CORES)
int <- scan1(probs, y, K, addcovar = covar, intcovar = intc, cores = CORES)
say("scans done (%.1f min)", as.numeric(difftime(Sys.time(), t0, units = "mins")))
perm_file <- file.path(out_dir, sprintf("qtl2_%s_perm.tsv", trait))
thr_add <- NA; thr_int <- NA
if (NPERM > 0) {
  if (file.exists(perm_file) && nrow(fread(perm_file)) >= NPERM) {
    pm <- fread(perm_file, data.table = FALSE)
  } else {
    strata <- setNames(as.integer(factor(taxon)), ids)
    pa <- scan1perm(probs, y, K, addcovar = covar, n_perm = NPERM, perm_strata = strata, cores = CORES)
    pi <- scan1perm(probs, y, K, addcovar = covar, intcovar = intc, n_perm = NPERM, perm_strata = strata, cores = CORES)
    pm <- data.frame(perm = seq_len(NPERM), max_lod_additive = pa[, 1], max_lod_per_taxon = pi[, 1])
    fwrite(pm, perm_file, sep = "\t")
  }
  thr_add <- unname(quantile(pm$max_lod_additive, 0.95)); thr_int <- unname(quantile(pm$max_lod_per_taxon, 0.95))
  say("5%% thresholds (%d permutations within taxon): additive LOD %.2f (-log10 p %.2f), per-taxon LOD %.2f",
      nrow(pm), thr_add, lod2p(thr_add, 1), thr_int)
} else if (THRESHOLD_LOD != "") {
  thr_add <- as.numeric(THRESHOLD_LOD); thr_int <- NA
  say("no permutations; additive threshold given: LOD %.2f (-log10 p %.2f)", thr_add, lod2p(thr_add, 1))
} else if (PEAKS == "") stop("NPERM=0 needs PEAKS=<step-27 peaks file> or THRESHOLD_LOD=")

sc <- data.frame(Marker = rownames(add), lod_additive = round(add[, 1], 3), lod_per_taxon = round(int[, 1], 3))
sc$lod_taxa_differ <- round(sc$lod_per_taxon - sc$lod_additive, 3)
sc$mlog10p_additive <- round(lod2p(sc$lod_additive, 1), 3)
sc <- merge(map[, c("Marker", "chr", "start", "end", "pos", "n_carriers")], sc, by = "Marker", sort = FALSE)
fwrite(sc, file.path(out_dir, sprintf("qtl2_%s_scan.tsv", trait)), sep = "\t")

## ---------- QTL ----------
if (PEAKS != "") {
  ## the same QTL as the k-mer scan (step 26/27): its top windows and support intervals
  if (!file.exists(PEAKS)) { say("no peaks file %s: nothing to fit", PEAKS); quit(save = "no") }
  kp <- fread(PEAKS, data.table = FALSE)
  if (!nrow(kp)) { say("peaks file %s is empty", PEAKS); quit(save = "no") }
  bcol <- intersect(c("BETA_PER_COPY", "BETA_DAYS_PER_COPY"), names(kp))[1]
  kp <- kp[kp$TOP_WINDOW %in% map$Marker, ]
  pk <- data.frame(chr = sub("^chr", "", kp$CHR), pos = map$pos[match(kp$TOP_WINDOW, map$Marker)] / 1e6,
                   ci_lo = kp$INTERVAL_START / 1e6, ci_hi = kp$INTERVAL_END / 1e6, found_by = "k-mer scan (step 26)",
                   Marker = kp$TOP_WINDOW, gridlmm_mlog10p = kp$TOP_MLOG10P, gridlmm_beta = kp[[bcol]], stringsAsFactors = FALSE)
  pk$lod <- sc$lod_additive[match(pk$Marker, sc$Marker)]
} else {
  peaks <- function(s, thr, lab) { if (is.na(thr)) return(NULL)
    p <- find_peaks(s, pmap, threshold = thr, peakdrop = PEAKDROP, drop = DROP)
    if (nrow(p)) p$found_by <- lab; p }
  pk <- rbind(peaks(add, thr_add, "additive"), peaks(int, thr_int, "per-taxon"))
  if (is.null(pk) || !nrow(pk)) { say("no QTL above the thresholds"); quit(save = "no") }
  pk$Marker <- mapply(function(ch, p) find_marker(pmap, ch, p), pk$chr, pk$pos)
  ## a per-taxon peak inside an additive peak's interval is the same QTL
  pk <- pk[order(pk$found_by != "additive", -pk$lod), ]
  keep <- rep(TRUE, nrow(pk))
  for (i in seq_len(nrow(pk))) for (j in seq_len(i - 1)) if (keep[j] && pk$chr[i] == pk$chr[j] &&
      pk$pos[i] >= pk$ci_lo[j] && pk$pos[i] <= pk$ci_hi[j]) keep[i] <- FALSE
  pk <- pk[keep, ]
  pk$gridlmm_mlog10p <- NA; pk$gridlmm_beta <- NA
}
Q <- pk$Marker

## ---------- more QTL inside the significant regions? conditional scans (forward selection) ----------
COND_MLOG10P <- Sys.getenv("COND_MLOG10P", ""); MAX_ADD <- envn("MAX_ADDED_QTL", 10); R2_MAX <- envn("MAX_R2_WITH_QTL", 0.5)
thr_cond <- if (COND_MLOG10P != "") qchisq(10^-as.numeric(COND_MLOG10P), 1, lower.tail = FALSE) / (2 * log(10)) else thr_add
GS <- Sys.getenv("GRIDLMM_SCAN", "")
gs <- if (GS != "" && file.exists(GS)) fread(GS, data.table = FALSE) else NULL
if (!is.na(thr_cond) && MAX_ADD > 0) {
  marg <- setNames(sc$lod_additive, sc$Marker)[colnames(G)]               # ordinary scan: only significant windows
  added <- list(); cs <- NULL
  for (rnd in seq_len(MAX_ADD)) {
    cs <- scan1(probs, y, K, addcovar = cbind(covar, G[, Q, drop = FALSE] / 2), cores = CORES)
    r2 <- suppressWarnings(apply(cor(G, G[, Q, drop = FALSE])^2, 1, max)); r2[is.na(r2)] <- 1
    ok <- rownames(cs)[marg[rownames(cs)] > thr_cond & r2[rownames(cs)] <= R2_MAX]
    if (!length(ok)) break
    best <- ok[which.max(cs[ok, 1])]
    if (cs[best, 1] <= thr_cond) break
    ch <- sub("^chr", "", map$chr[match(best, map$Marker)]); pos <- map$pos[match(best, map$Marker)] / 1e6
    cs2 <- cs; cs2[!rownames(cs2) %in% ok, 1] <- 0                         # support interval around this window
    ci <- tryCatch(lod_int(cs2, pmap, chr = ch, drop = DROP), error = function(e) NULL)
    lo <- if (is.null(ci)) pos else ci[1, "ci_lo"]; hi <- if (is.null(ci)) pos else ci[1, "ci_hi"]
    added[[rnd]] <- data.frame(trait = trait, round = rnd, qtl = best, chr = paste0("chr", ch), pos_Mb = pos,
                               lod_conditional = round(cs[best, 1], 2), mlog10p_conditional = round(lod2p(cs[best, 1], 1), 2),
                               max_r2_with_earlier_qtl = round(r2[best], 2))
    say("conditional scan %d: added %s (conditional LOD %.2f, -log10 p %.2f)", rnd, best, cs[best, 1], lod2p(cs[best, 1], 1))
    gi <- if (!is.null(gs)) match(best, gs$Marker) else NA
    pk <- rbindlist(list(pk, data.frame(chr = ch, pos = pos, ci_lo = lo, ci_hi = hi, found_by = "conditional scan",
                                        Marker = best, lod = sc$lod_additive[match(best, sc$Marker)],
                                        gridlmm_mlog10p = if (is.na(gi)) NA else -log10(gs$p_add[gi]),
                                        gridlmm_beta = if (is.na(gi)) NA else gs$beta_add[gi])), fill = TRUE)
    pk <- as.data.frame(pk); Q <- c(Q, best)
  }
  if (length(added)) fwrite(rbindlist(added), file.path(out_dir, sprintf("qtl2_%s_added.tsv", trait)), sep = "\t")
  if (!is.null(cs)) fwrite(data.frame(Marker = rownames(cs), lod_conditional = round(cs[, 1], 3)),
                           file.path(out_dir, sprintf("qtl2_%s_conditional.tsv", trait)), sep = "\t")
  say("conditional scans: %d QTL added (threshold LOD %.2f)", length(added), thr_cond)
}
pk <- pk[order(match(pk$chr, cn), pk$pos), ]
Q <- pk$Marker
say("QTL: %s", paste(sprintf("%s (%s)", Q, pk$found_by), collapse = ", "))

## ---------- joint model and per-taxon effects (fit1) ----------
ctr <- cbind(mu = c(1, 1), a = c(-1, 1))                      # a = effect per teosinte copy
qdose <- G[, Q, drop = FALSE] / 2                             # teosinte allele probability at each QTL
yv <- var(y[, 1])
rows <- list(); trows <- list()
for (i in seq_along(Q)) {
  q <- Q[i]; ch <- pk$chr[i]
  pr_q <- pull_genoprobpos(probs, map = pmap, chr = ch, pos = pk$pos[i])
  others <- qdose[, setdiff(Q, q), drop = FALSE]
  f <- fit1(pr_q, y, K[[ch]], addcovar = cbind(covar, others), contrasts = ctr)
  fs <- fit1(pr_q, y, K[[ch]], addcovar = covar, contrasts = ctr)           # alone, for PVE
  pve <- var(G[, q] * fs$coef["a"]) / yv
  li <- sc$lod_taxa_differ[sc$Marker == q]
  rows[[q]] <- data.frame(trait = trait, qtl = q, chr = paste0("chr", ch), pos_Mb = pk$pos[i],
                          interval_start_Mb = pk$ci_lo[i], interval_end_Mb = pk$ci_hi[i], found_by = pk$found_by[i],
                          carriers = sum(G[, q] > 0), gridlmm_mlog10p = round(pk$gridlmm_mlog10p[i], 2),
                          gridlmm_effect_per_copy = round(pk$gridlmm_beta[i], 3),
                          lod = round(sc$lod_additive[sc$Marker == q], 2), mlog10p = round(lod2p(sc$lod_additive[sc$Marker == q], 1), 2),
                          lod_joint = round(f$lod, 2), mlog10p_joint = round(lod2p(f$lod, 1), 2),
                          effect_per_copy_joint = round(f$coef["a"], 3), se = round(f$SE["a"], 3),
                          lod_taxa_differ = round(li, 2), mlog10p_taxa_differ = round(lod2p(li, length(taxa) - 1), 2),
                          pve_alone = round(pve, 4), threshold_lod = round(thr_add, 2))
  for (t in taxa) {
    s <- names(taxon)[taxon == t]; car <- sum(G[s, q] > 0)
    r <- data.frame(trait = trait, qtl = q, taxon = t, plots = length(s), carriers = car,
                    seg = sum(G[s, q] == 1), donor = sum(G[s, q] == 2), effect_per_copy = NA, se = NA, lod = NA)
    if (car >= MIN_TC) {
      oc <- qdose[s, setdiff(Q, q), drop = FALSE]                    # other QTL that vary in this taxon
      oc <- oc[, apply(oc, 2, function(x) length(unique(x)) > 1 && abs(suppressWarnings(cor(x, G[s, q]))) < 0.9), drop = FALSE]
      ft <- tryCatch(fit1(pr_q[s, , drop = FALSE], y[s, , drop = FALSE], K[[ch]][s, s],
                          addcovar = if (ncol(oc)) oc else NULL, contrasts = ctr),
                     error = function(e) NULL)
      if (!is.null(ft)) { r$effect_per_copy <- round(ft$coef["a"], 3); r$se <- round(ft$SE["a"], 3); r$lod <- round(ft$lod, 2) }
    }
    trows[[paste(q, t)]] <- r
  }
}
R <- rbindlist(rows); TB <- rbindlist(trows)
fwrite(R, file.path(out_dir, sprintf("qtl2_%s_qtl.tsv", trait)), sep = "\t")
fwrite(TB, file.path(out_dir, sprintf("qtl2_%s_by_taxon.tsv", trait)), sep = "\t")

## ---------- genes in the intervals ----------
if (GFF != "" && file.exists(GFF)) {
  gf <- fread(cmd = sprintf("grep -v '^#' %s | awk -F'\\t' '$3==\"gene\"'", shQuote(GFF)), header = FALSE, sep = "\t", quote = "")
  gf[, chr := sub("^chr", "", V1)]
  gl <- rbindlist(lapply(seq_len(nrow(R)), function(i) {
    x <- gf[chr == sub("^chr", "", R$chr[i]) & V5 >= R$interval_start_Mb[i] * 1e6 & V4 <= R$interval_end_Mb[i] * 1e6]
    if (!nrow(x)) return(NULL)
    data.table(trait = trait, qtl = R$qtl[i], gene = sub(".*(gene_id|ID)=(gene:)?([^;]+).*", "\\3", x$V9),
               chr = R$chr[i], gene_start = x$V4, gene_end = x$V5,
               name = ifelse(grepl("Name=", x$V9), sub(".*Name=([^;]+).*", "\\1", x$V9), ""),
               description = ifelse(grepl("description=", x$V9), sapply(sub(".*description=([^;]+).*", "\\1", x$V9), URLdecode), ""))
  }))
  if (nrow(gl)) fwrite(gl, file.path(out_dir, sprintf("qtl2_%s_genes.tsv", trait)), sep = "\t")
}

## ---------- plot ----------
png(file.path(out_dir, sprintf("qtl2_%s_manhattan.png", trait)), width = 2400, height = 1400, res = 200)
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
plot(add, pmap, col = "black", altcol = "grey50", main = sprintf("%s: additive (taxon + LOCO K), %d plots", trait, length(ids)))
if (!is.na(thr_add)) abline(h = thr_add, lty = 3)
plot(int, pmap, col = "black", altcol = "grey50", main = sprintf("%s: per-taxon QTL effects (intcovar = taxon)", trait))
if (!is.na(thr_int)) abline(h = thr_int, lty = 3)
invisible(dev.off())

print(as.data.frame(R[, .(qtl, interval_start_Mb, interval_end_Mb, found_by, lod, lod_joint, effect_per_copy_joint, se, lod_taxa_differ)]), row.names = FALSE)
print(as.data.frame(TB[, .(qtl, taxon, carriers, effect_per_copy, se, lod)]), row.names = FALSE)
say("%.1f min. WROTE %s/qtl2_%s_{qtl,by_taxon,genes,scan,perm}.tsv and _manhattan.png",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, trait)
