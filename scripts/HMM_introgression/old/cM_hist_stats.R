

suppressPackageStartupMessages({
  library(data.table)
})

# ---------- USER SETTINGS ----------
state_dir <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/statepaths"
map_dir   <- "/rsstu/users/r/rrellan/sara/ref/NAM_genetic_map"

# Option B knob: bridge islands shorter than this (in cM)
cm_gap_RH_inside_donor <- 0.05   # try 0.02, 0.05, 0.10 based on your data
cm_gap_RR_inside_donor <- 0.05   # optional; can set 0 to disable

# Donor definition
donor_mode <- c("HH_only","HH_RH")  # compute both

out_stats_tsv <- "all_samples.state_stats.cmBridge.tsv"
out_pdf       <- "BC2S3_expected_vs_observed_hist.cmBridge.pdf"

# ---------- EXPECTED BC2S3 ----------
exp_HH <- 100 * (1/8) * (7/16)      # 5.46875
exp_RH <- 100 * (1/8) * (1/8)       # 1.5625
exp_HH_RH <- exp_HH + exp_RH        # 7.03125
exp_RR <- 100 - exp_HH_RH           # 92.96875

# ---------- MAP HELPERS ----------
read_map <- function(chr) {
  f <- file.path(map_dir, paste0(chr, ".map"))
  m <- fread(f, header=FALSE)
  # chrN.map is: bp 1.0 cM  (we want bp + cM)
  setnames(m, c("bp","col2","cM"))
  m <- m[order(bp)]
  m
}

interp_cm <- function(pos, m) {
  # linear interpolation with edge clamp
  bp <- m$bp; cm <- m$cM
  if (pos <= bp[1]) return(cm[1])
  if (pos >= bp[length(bp)]) return(cm[length(cm)])
  i <- findInterval(pos, bp)
  x0 <- bp[i]; x1 <- bp[i+1]
  y0 <- cm[i]; y1 <- cm[i+1]
  if (x1 == x0) return(y0)
  y0 + (pos-x0) * (y1-y0) / (x1-x0)
}

# ---------- RUN-LENGTH / BRIDGING ----------
rle_segments <- function(states, weights) {
  # returns data.table of runs with run_state and run_weight
  r <- rle(states)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  run_w <- vapply(seq_along(starts), function(i) sum(weights[starts[i]:ends[i]]), numeric(1))
  data.table(run_id=seq_along(r$values), state=r$values, w=run_w, start=starts, end=ends)
}

bridge_islands <- function(states, weights, donor_set, cm_gap_island, island_states) {
  # bridges islands of specified state(s) inside donor, if island run weight < cm_gap
  # only bridges if island is flanked by donor on BOTH sides
  seg <- rle_segments(states, weights)
  new_states <- states

  for (i in seq_len(nrow(seg))) {
    st <- seg$state[i]
    if (!(st %in% island_states)) next
    if (seg$w[i] >= cm_gap_island) next
    if (i == 1 || i == nrow(seg)) next

    left_ok  <- seg$state[i-1] %in% donor_set
    right_ok <- seg$state[i+1] %in% donor_set
    if (left_ok && right_ok) {
      # relabel this entire island to DONOR (pick a representative donor state)
      idx <- seg$start[i]:seg$end[i]
      # If donor_set is HH only -> relabel to HH; if HH+RH -> relabel to HH (safer)
      new_states[idx] <- "HH"
    }
  }
  new_states
}

# ---------- PER SAMPLE COMPUTATION ----------
compute_sample <- function(f_statepath) {
  sample <- sub("\\.statepath\\.tsv\\.gz$", "", basename(f_statepath))

  dt <- fread(cmd=paste("zcat", shQuote(f_statepath)))
  # Expect columns: CHROM POS REF ALT DP STATE ...
  dt <- dt[, .(CHROM, POS, STATE)]
  dt <- dt[CHROM %in% paste0("chr",1:10)]
  setorder(dt, CHROM, POS)

  # load maps once per chr used
  chrs <- unique(dt$CHROM)
  maps <- lapply(chrs, read_map); names(maps) <- chrs

  # interpolate cM per SNP
  dt[, cM := mapply(function(chr,pos) interp_cm(pos, maps[[chr]]), CHROM, POS)]

  # build interval weights between SNP i and i+1
  # weight is delta cM to next SNP, within same chr; last SNP on chr gets 0
  dt[, next_cM := shift(cM, type="lead"), by=CHROM]
  dt[, w_cM := pmax(0, next_cM - cM)]
  dt[is.na(w_cM), w_cM := 0]

  total_cM <- sum(dt$w_cM)

  # baseline (no bridge) percentages by cM
  pct_by_state <- function(states_vec) {
    sRR <- sum(dt$w_cM[states_vec == "RR"])
    sRH <- sum(dt$w_cM[states_vec == "RH"])
    sHH <- sum(dt$w_cM[states_vec == "HH"])
    c(pct_RR=100*sRR/total_cM, pct_RH=100*sRH/total_cM, pct_HH=100*sHH/total_cM)
  }

  base_pct <- pct_by_state(dt$STATE)

  out <- list(
    sample=sample,
    total_cM=total_cM,
    pct_RR_cM=base_pct["pct_RR"],
    pct_RH_cM=base_pct["pct_RH"],
    pct_HH_cM=base_pct["pct_HH"]
  )

  # donor definitions + bridging
  for (mode in donor_mode) {
    donor_set <- if (mode=="HH_only") c("HH") else c("HH","RH")

    st <- dt$STATE

    # bridge RH islands inside donor
    st2 <- bridge_islands(
      states=st, weights=dt$w_cM,
      donor_set=donor_set,
      cm_gap_island=cm_gap_RH_inside_donor,
      island_states=c("RH")
    )

    # optionally bridge RR islands inside donor (often helpful for occasional RR noise)
    if (cm_gap_RR_inside_donor > 0) {
      st2 <- bridge_islands(
        states=st2, weights=dt$w_cM,
        donor_set=c("HH","RH"),  # after first bridge, donor-like is HH/RH
        cm_gap_island=cm_gap_RR_inside_donor,
        island_states=c("RR")
      )
    }

    donor_cM <- sum(dt$w_cM[st2 %in% donor_set])
    donor_pct <- 100 * donor_cM / total_cM

    nm <- if (mode=="HH_only") "pct_donor_cM_HH" else "pct_donor_cM_HH_RH"
    out[[nm]] <- donor_pct
  }

  as.data.table(out)
}

# ---------- RUN ALL ----------
files <- sort(list.files(state_dir, pattern="\\.statepath\\.tsv\\.gz$", full.names=TRUE))
if (length(files)==0) stop("No statepath files found")

res <- rbindlist(lapply(files, compute_sample), fill=TRUE)
fwrite(res, out_stats_tsv, sep="\t")

summ <- function(x) c(mean=mean(x,na.rm=TRUE), median=median(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE),
                      q05=as.numeric(quantile(x,0.05,na.rm=TRUE)),
                      q95=as.numeric(quantile(x,0.95,na.rm=TRUE)))

cat("Expected BC2S3 (% genome):\n")
cat(sprintf("RR ~ %.4f, RH ~ %.4f, HH ~ %.4f, (HH+RH) ~ %.4f\n\n", exp_RR, exp_RH, exp_HH, exp_HH_RH))

cat("Observed (% genome, cM-based, with cM-bridging):\n")
print(rbind(
  RR    = summ(res$pct_RR_cM),
  RH    = summ(res$pct_RH_cM),
  HH    = summ(res$pct_donor_cM_HH),
  HH_RH = summ(res$pct_donor_cM_HH_RH)
))

pdf(out_pdf, width=10, height=7)

hist(res$pct_donor_cM_HH, breaks=50,
     main=paste0("Donor % (HH only) with cM-bridging (gap=",cm_gap_RH_inside_donor," cM)"),
     xlab="% genome donor (cM)", ylab="Number of samples")
abline(v=exp_HH, lwd=3)

hist(res$pct_donor_cM_HH_RH, breaks=50,
     main=paste0("Donor % (HH+RH) with cM-bridging (gap=",cm_gap_RH_inside_donor," cM)"),
     xlab="% genome donor (cM)", ylab="Number of samples")
abline(v=exp_HH_RH, lwd=3)

hist(res$pct_RH_cM, breaks=50,
     main="RH % (cM-based) after bridging",
     xlab="% genome RH (cM)", ylab="Number of samples")
abline(v=exp_RH, lwd=3)

dev.off()

cat("\nWrote:\n")
cat("  ", out_stats_tsv, "\n")
cat("  ", out_pdf, "\n")
