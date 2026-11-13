#!/usr/bin/env Rscript
# =============================================================================
# 14_call_introgression.R — call introgressions per pool from window read counts (base R)
#
#   Usage: Rscript 14_call_introgression.R <D13> <pools.tsv> <ref.fai> <outdir> <win_size> <switch_prob> <chroms,comma,separated>
#
#   Each pool = one BC2S3 family. Every window gets one of three states:
#     B73    fixed B73 in the family
#     SEG    segregating in the family (some plants carry donor)
#     DONOR  fixed donor in the family
#   In a window, ALT reads come from:  background (errors; B73 controls)  +  dose x q
#     dose = 0 (B73), 0.5 (SEG), 1 (DONOR)
#     q    = how often this donor carries the ALT allele at the window's sites
#   q is estimated per window and pool, iteratively (N_ROUNDS):
#     start : one value per taxon, from windows that clearly carry donor
#     round : call all pools; then for each window, q = ALT excess of the OTHER pools of the
#             same accession that were called DONOR there, blended with the taxon's value for
#             that window (weight N_PRIOR reads) so a few reads cannot swing it; call again
#   Reads are modeled as beta-binomial (allows extra spread: SEG pools vary in how many
#   plants carry the donor; donor haplotypes vary site to site). States along a chromosome
#   are linked by a hidden Markov model (switch probability per window), run forward-backward.
#
#   Outputs in <outdir>:
#     calls.tsv.gz         window x pool: reads, ALT, q, P(B73) P(SEG) P(DONOR), DOSAGE, CALL
#     dosage_matrix.tsv.gz windows x pools, expected donor dosage (0..1) — for GWAS / QTL
#     segments.tsv         merged runs of SEG / DONOR per pool
#     pool_summary.tsv     per pool: donor share of the genome, SEG / DONOR fractions, segments
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) stop("usage: Rscript 14_call_introgression.R <D13> <pools.tsv> <fai> <outdir> <win> <switch> <chroms>")
d13 <- args[1]; pools_file <- args[2]; fai_file <- args[3]; outdir <- args[4]
W <- as.numeric(args[5]); sw <- as.numeric(args[6]); chroms <- strsplit(args[7], ",")[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

STATES <- c("B73", "SEG", "DONOR"); DOSE <- c(0, 0.5, 1)
PRIOR  <- c(0.844, 0.0625, 0.094)       # BC2S3 family expectation (fixed B73 / segregating / fixed donor)
PHI    <- c(200, 8, 30)                  # beta-binomial concentration per state (smaller = more spread)
N_ROUNDS      <- 3                       # rounds of: call -> re-estimate q -> call
Q_MIN_FACTOR  <- 3                       # window uninformative if q < 3 x background

# ---- read window counts --------------------------------------------------------
pools <- read.delim(pools_file, stringsAsFactors = FALSE)
pools <- pools[pools$ROLE %in% c("control", "donor"), ]
fai <- read.delim(fai_file, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
names(fai) <- c("CHROM", "LEN")

wc <- do.call(rbind, lapply(chroms, function(c) {
  f <- file.path(d13, paste0("BZea_", c, ".win_counts.tsv.gz"))
  if (!file.exists(f)) stop("missing ", f, " (run 13)")
  read.delim(f, stringsAsFactors = FALSE)
}))
wc <- wc[wc$SAMPLE %in% pools$SAMPLE, ]
message("window rows: ", nrow(wc), " | pools: ", length(unique(wc$SAMPLE)))

# full window grid for every chromosome (windows without sites get 0 reads)
grid <- do.call(rbind, lapply(chroms, function(c) {
  len <- fai$LEN[fai$CHROM == c]; s <- seq(1, len, by = W)
  data.frame(CHROM = c, START = s, END = pmin(s + W - 1, len), stringsAsFactors = FALSE)
}))
grid$WID <- seq_len(nrow(grid))
samples <- pools$SAMPLE
P <- length(samples); NW <- nrow(grid)
key <- paste(grid$CHROM, grid$START)
wi <- match(paste(wc$CHROM, wc$START), key)
si <- match(wc$SAMPLE, samples)
N <- matrix(0, NW, P); K <- matrix(0, NW, P)
N[cbind(wi, si)] <- wc$REF + wc$ALT
K[cbind(wi, si)] <- wc$ALT
colnames(N) <- colnames(K) <- samples

is_ctrl <- pools$ROLE == "control"
grp <- pools$GROUP; acc <- pools$ACCESSION

# ---- background ALT rate (B73 controls, genome-wide) -----------------------------
bg <- max(sum(K[, is_ctrl]) / max(sum(N[, is_ctrl]), 1), 1e-4)
message(sprintf("background ALT rate (B73 controls): %.5f", bg))

# ---- HMM ----------------------------------------------------------------------------
lbb <- function(k, n, m, phi) {                     # beta-binomial log-lik without the choose term
  a <- m * phi; b <- (1 - m) * phi
  lbeta(k + a, n - k + b) - lbeta(a, b)
}
logT <- matrix(log(sw / 2), 3, 3); diag(logT) <- log(1 - sw)
lse3 <- function(x1, x2, x3) { m <- pmax(x1, x2, x3); m + log(exp(x1 - m) + exp(x2 - m) + exp(x3 - m)) }

run_hmm <- function(Q) {        # Q: NW x P matrix of donor ALT rate -> list of posteriors (NW x P each)
  E <- lapply(1:3, function(s) {
    m <- pmin(bg + DOSE[s] * Q, 0.999)
    e <- lbb(K, N, m, PHI[s])
    e[Q < Q_MIN_FACTOR * bg] <- 0                    # uninformative window: no evidence either way
    e
  })
  post <- lapply(1:3, function(s) matrix(0, NW, P))
  for (c in chroms) {
    idx <- which(grid$CHROM == c); L <- length(idx)
    fa <- lapply(1:3, function(s) matrix(0, L, P)); fb <- fa
    for (s in 1:3) fa[[s]][1, ] <- log(PRIOR[s]) + E[[s]][idx[1], ]
    if (L > 1) for (t in 2:L) for (s in 1:3)
      fa[[s]][t, ] <- lse3(fa[[1]][t - 1, ] + logT[1, s], fa[[2]][t - 1, ] + logT[2, s], fa[[3]][t - 1, ] + logT[3, s]) + E[[s]][idx[t], ]
    if (L > 1) for (t in (L - 1):1) for (s in 1:3)
      fb[[s]][t, ] <- lse3(logT[s, 1] + E[[1]][idx[t + 1], ] + fb[[1]][t + 1, ],
                           logT[s, 2] + E[[2]][idx[t + 1], ] + fb[[2]][t + 1, ],
                           logT[s, 3] + E[[3]][idx[t + 1], ] + fb[[3]][t + 1, ])
    z <- lse3(fa[[1]] + fb[[1]], fa[[2]] + fb[[2]], fa[[3]] + fb[[3]])
    for (s in 1:3) post[[s]][idx, ] <- exp(fa[[s]] + fb[[s]] - z)
  }
  post
}

# ---- donor ALT rate q: start value, then iterate -----------------------------------
# start: per taxon, one genome-wide value from windows that clearly carry donor
X <- K / pmax(N, 1) - bg                               # ALT excess per window x pool
donor_cols <- which(!is_ctrl)
taxa <- sort(unique(grp[donor_cols]))
q_init <- sapply(taxa, function(g) {
  cols <- which(grp == g & !is_ctrl)
  v <- X[, cols][N[, cols] >= 100]
  v <- v[v > 10 * bg]
  if (length(v) < 20) return(0.05)
  unname(quantile(v, 0.8))
})
Q <- matrix(0, NW, P)
for (g in taxa) Q[, grp == g & !is_ctrl] <- q_init[[g]]
Q[, is_ctrl] <- median(q_init)                          # controls are called like any pool
message("start q per taxon: ", paste(sprintf("%s=%.3f", taxa, q_init), collapse = " "))

# each round: call -> re-estimate q per taxon x window from pools called DONOR
#             -> per pool from the OTHER pools of its accession called DONOR, blended
#                with the taxon value (N_PRIOR reads of weight) -> call again
N_PRIOR <- 300
for (round in 1:N_ROUNDS) {
  post <- run_hmm(Q)
  dn <- post[[3]] > 0.8
  for (g in taxa) {
    cols <- which(grp == g & !is_ctrl)
    sk <- rowSums(K[, cols, drop = FALSE] * dn[, cols, drop = FALSE])
    sn <- rowSums(N[, cols, drop = FALSE] * dn[, cols, drop = FALSE])
    qg_all <- if (sum(sn) > 0) max(sum(sk) / sum(sn) - bg, 0) else q_init[[g]]
    qg <- (sk + N_PRIOR * (qg_all + bg)) / (sn + N_PRIOR) - bg          # taxon x window, shrunk to taxon mean
    for (a in unique(acc[cols])) {
      ac <- cols[acc[cols] == a]
      ka <- rowSums(K[, ac, drop = FALSE] * dn[, ac, drop = FALSE])
      na <- rowSums(N[, ac, drop = FALSE] * dn[, ac, drop = FALSE])
      for (j in ac) {
        k_o <- ka - K[, j] * dn[, j]; n_o <- na - N[, j] * dn[, j]      # leave this pool out
        Q[, j] <- pmax((k_o + N_PRIOR * (qg + bg)) / (n_o + N_PRIOR) - bg, 0)
      }
    }
    message(sprintf("round %d  %-3s pools %4d  q (taxon mean) %.3f  windows x pools called DONOR %.1f%%",
                    round, g, length(cols), qg_all, 100 * mean(dn[, cols])))
  }
}
post <- run_hmm(Q)
q_acc <- Q

# ---- outputs ---------------------------------------------------------------------
dosage <- 0.5 * post[[2]] + post[[3]]
callm <- matrix(STATES[max.col(cbind(c(post[[1]]), c(post[[2]]), c(post[[3]])), ties.method = "first")], NW, P)
colnames(dosage) <- colnames(callm) <- samples

long <- data.frame(grid[rep(seq_len(NW), P), c("CHROM", "START", "END")],
                   SAMPLE = rep(samples, each = NW), READS = c(N), ALT = c(K), Q = round(c(q_acc), 4),
                   P_B73 = round(c(post[[1]]), 4), P_SEG = round(c(post[[2]]), 4), P_DONOR = round(c(post[[3]]), 4),
                   DOSAGE = round(c(dosage), 4), CALL = c(callm), stringsAsFactors = FALSE)
gz <- gzfile(file.path(outdir, "calls.tsv.gz"), "w"); write.table(long, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

dm <- data.frame(grid[, c("CHROM", "START", "END")], round(dosage, 4), check.names = FALSE)
gz <- gzfile(file.path(outdir, "dosage_matrix.tsv.gz"), "w"); write.table(dm, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

seg <- do.call(rbind, lapply(seq_len(P), function(j) do.call(rbind, lapply(chroms, function(c) {
  idx <- which(grid$CHROM == c); r <- rle(callm[idx, j]); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  keep <- r$values != "B73"
  if (!any(keep)) return(NULL)
  data.frame(SAMPLE = samples[j], CHROM = c, START = grid$START[idx[s[keep]]], END = grid$END[idx[e[keep]]],
             STATE = r$values[keep], N_WINDOWS = r$lengths[keep],
             MEAN_DOSAGE = round(sapply(which(keep), function(i) mean(dosage[idx[s[i]:e[i]], j])), 3))
}))))
if (is.null(seg)) seg <- data.frame(SAMPLE = character(), CHROM = character(), START = numeric(), END = numeric(),
                                    STATE = character(), N_WINDOWS = numeric(), MEAN_DOSAGE = numeric())
write.table(seg, file.path(outdir, "segments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

wlen <- grid$END - grid$START + 1
ps <- data.frame(pools[, c("SAMPLE", "ROLE", "GROUP", "ACCESSION")],
                 READS = colSums(N),
                 DONOR_SHARE = round(colSums(dosage * wlen) / sum(wlen), 4),
                 FRAC_DONOR = round(colSums((callm == "DONOR") * wlen) / sum(wlen), 4),
                 FRAC_SEG = round(colSums((callm == "SEG") * wlen) / sum(wlen), 4),
                 N_SEGMENTS = as.integer(table(factor(seg$SAMPLE, levels = samples))))
write.table(ps, file.path(outdir, "pool_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDonor share of the genome per pool (expected ~0.125 for BC2-derived lines; ~0 for B73 controls):\n")
print(do.call(rbind, lapply(split(ps, ps$GROUP), function(d) data.frame(
  GROUP = d$GROUP[1], N = nrow(d), MEDIAN_DONOR_SHARE = median(d$DONOR_SHARE),
  MEDIAN_FRAC_DONOR = median(d$FRAC_DONOR), MEDIAN_FRAC_SEG = median(d$FRAC_SEG),
  MEDIAN_SEGMENTS = median(d$N_SEGMENTS)))), row.names = FALSE)
cat("\nWrote:", file.path(outdir, c("calls.tsv.gz", "dosage_matrix.tsv.gz", "segments.tsv", "pool_summary.tsv")), sep = "\n  ")
