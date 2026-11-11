suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# TEST RIGIDITY ON SINGLE SAMPLE
# ============================================================================

input_file <- "/Users/nirwantandukar/Documents/Research/data/BZea/angsd_genotyping/GL/PN17_SID1627.statepath.tsv.gz"
out_dir <- "/Users/nirwantandukar/Documents/Github/BZea_genotyping/scripts/HMM_introgression"

cat("=== RIGIDITY TEST FOR PN17_SID1627 ===\n\n")

# Read data
cat("Reading data...\n")
dt <- fread(cmd = paste("gunzip -c", shQuote(input_file)), showProgress = FALSE)
cat("Total SNPs:", nrow(dt), "\n\n")

# Original state counts
cat("=== ORIGINAL STATE COUNTS ===\n")
print(table(dt$STATE))
cat("\n")

# Chromosome sizes
chr_len <- c(308452471, 243675191, 238017767, 250330460, 226353449,
             181357234, 185808916, 182411202, 163004744, 152435371)
names(chr_len) <- paste0("chr", 1:10)

# Helper: get posterior for called state
get_state_prob <- function(state, p_rr, p_rh, p_hh) {
  ifelse(state == "RR", p_rr, ifelse(state == "RH", p_rh, p_hh))
}

# Smoothing function
smooth_states <- function(states, probs_matrix, min_markers) {
  if (length(states) == 0 || min_markers <= 1) return(states)

  changed <- TRUE
  iter <- 0

  while (changed && iter < 50) {
    changed <- FALSE
    iter <- iter + 1

    r <- rle(states)
    n_runs <- length(r$values)
    if (n_runs <= 1) break

    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L

    for (i in seq_len(n_runs)) {
      if (r$lengths[i] < min_markers) {
        seg_idx <- starts[i]:ends[i]

        log_probs <- log(probs_matrix[seg_idx, , drop = FALSE] + 1e-10)
        state_scores <- colSums(log_probs)

        left_state <- if (i > 1) r$values[i - 1] else NULL
        right_state <- if (i < n_runs) r$values[i + 1] else NULL
        bonus <- 0.5 * r$lengths[i]

        if (!is.null(left_state)) {
          state_scores[left_state] <- state_scores[left_state] + bonus
        }
        if (!is.null(right_state)) {
          state_scores[right_state] <- state_scores[right_state] + bonus
        }

        best <- names(which.max(state_scores))

        if (best != r$values[i]) {
          states[seg_idx] <- best
          changed <- TRUE
        }
      }
    }
  }

  return(states)
}

# Analyze segments
analyze_segments <- function(states, chrom) {
  # Count segments by state
  r <- rle(states)
  seg_dt <- data.table(
    state = r$values,
    length = r$lengths,
    chr = chrom[cumsum(r$lengths)]
  )

  # Segment size distribution
  list(
    n_segments = nrow(seg_dt),
    by_state = table(seg_dt$state),
    size_dist = summary(seg_dt$length),
    short_segments = sum(seg_dt$length < 5),
    very_short = sum(seg_dt$length < 3)
  )
}

# Test different rigidity values
rigidity_values <- c(1, 3, 5, 10, 15, 20)

results <- list()

for (min_markers in rigidity_values) {
  cat("Testing min_markers =", min_markers, "...\n")

  # Process each chromosome
  all_states_smooth <- character(0)
  all_states_orig <- character(0)
  all_chrom <- character(0)

  for (chr in paste0("chr", 1:10)) {
    dchr <- dt[CHROM == chr]
    if (nrow(dchr) == 0) next

    orig_states <- dchr$STATE
    probs <- as.matrix(dchr[, .(P_RR, P_RH, P_HH)])
    colnames(probs) <- c("RR", "RH", "HH")

    smooth_st <- smooth_states(orig_states, probs, min_markers)

    all_states_orig <- c(all_states_orig, orig_states)
    all_states_smooth <- c(all_states_smooth, smooth_st)
    all_chrom <- c(all_chrom, rep(chr, length(orig_states)))
  }

  # Compare
  orig_counts <- table(all_states_orig)
  smooth_counts <- table(all_states_smooth)

  orig_seg <- analyze_segments(all_states_orig, all_chrom)
  smooth_seg <- analyze_segments(all_states_smooth, all_chrom)

  results[[as.character(min_markers)]] <- list(
    min_markers = min_markers,
    orig_counts = orig_counts,
    smooth_counts = smooth_counts,
    orig_segments = orig_seg$n_segments,
    smooth_segments = smooth_seg$n_segments,
    segments_removed = orig_seg$n_segments - smooth_seg$n_segments,
    orig_short = orig_seg$short_segments,
    smooth_short = smooth_seg$short_segments
  )
}

# ============================================================================
# DISPLAY RESULTS
# ============================================================================

cat("\n")
cat("=====================================================================\n")
cat("                    RIGIDITY COMPARISON RESULTS                      \n")
cat("=====================================================================\n\n")

# Create summary table
summary_dt <- rbindlist(lapply(results, function(r) {
  data.table(
    min_markers = r$min_markers,
    RR_snps = as.numeric(r$smooth_counts["RR"]),
    RH_snps = as.numeric(r$smooth_counts["RH"]),
    HH_snps = as.numeric(r$smooth_counts["HH"]),
    n_segments = r$smooth_segments,
    segments_removed = r$segments_removed,
    short_segs_remaining = r$smooth_short
  )
}))

# Add percentages
summary_dt[, RH_pct := round(100 * RH_snps / (RR_snps + RH_snps + HH_snps), 2)]
summary_dt[, HH_pct := round(100 * HH_snps / (RR_snps + RH_snps + HH_snps), 2)]
summary_dt[, donor_pct := round(RH_pct + HH_pct, 2)]

cat("SNP-level state counts after smoothing:\n")
print(summary_dt[, .(min_markers, RR_snps, RH_snps, HH_snps,
                     RH_pct, HH_pct, donor_pct)])

cat("\nSegment-level changes:\n")
print(summary_dt[, .(min_markers, n_segments, segments_removed,
                     short_segs_remaining)])

# ============================================================================
# DETAILED TRACT OUTPUT FOR BEST RIGIDITY
# ============================================================================

cat("\n")
cat("=====================================================================\n")
cat("             INTROGRESSION TRACTS (min_markers = 5)                  \n")
cat("=====================================================================\n\n")

# Use min_markers = 5 for tract output
min_markers <- 5

tracts_list <- list()

for (chr in paste0("chr", 1:10)) {
  dchr <- dt[CHROM == chr]
  if (nrow(dchr) == 0) next

  orig_states <- dchr$STATE
  probs <- as.matrix(dchr[, .(P_RR, P_RH, P_HH)])
  colnames(probs) <- c("RR", "RH", "HH")

  smooth_st <- smooth_states(orig_states, probs, min_markers)

  # Get smoothed tract info
  r <- rle(smooth_st)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  for (i in seq_along(r$values)) {
    if (r$values[i] %in% c("HH", "RH")) {
      seg_idx <- starts[i]:ends[i]

      tracts_list[[length(tracts_list) + 1]] <- data.table(
        Chr = chr,
        Start = dchr$POS[starts[i]],
        End = dchr$POS[ends[i]],
        State = r$values[i],
        n_SNPs = r$lengths[i],
        len_bp = dchr$POS[ends[i]] - dchr$POS[starts[i]] + 1,
        mean_posterior = mean(get_state_prob(
          smooth_st[seg_idx],
          dchr$P_RR[seg_idx],
          dchr$P_RH[seg_idx],
          dchr$P_HH[seg_idx]
        ))
      )
    }
  }
}

tracts <- rbindlist(tracts_list)
tracts[, len_Mb := round(len_bp / 1e6, 3)]

cat("Total donor tracts:", nrow(tracts), "\n")
cat("Total HH tracts:", sum(tracts$State == "HH"), "\n")
cat("Total RH tracts:", sum(tracts$State == "RH"), "\n\n")

cat("Tract size summary (Mb):\n")
print(summary(tracts$len_Mb))

cat("\nTop 20 largest introgression tracts:\n")
print(tracts[order(-len_bp)][1:20,
  .(Chr, Start, End, State, n_SNPs, len_Mb,
    mean_post = round(mean_posterior, 4))])

# Save tracts
out_tracts <- file.path(out_dir, "PN17_SID1627.tracts.rigidity5.tsv")
fwrite(tracts, out_tracts, sep = "\t")
cat("\nTracts saved to:", out_tracts, "\n")

# ============================================================================
# CHROMOSOME SUMMARY
# ============================================================================

cat("\n")
cat("=====================================================================\n")
cat("                    CHROMOSOME-LEVEL SUMMARY                         \n")
cat("=====================================================================\n\n")

chr_summary <- tracts[, .(
  n_tracts = .N,
  total_bp = sum(len_bp),
  total_Mb = round(sum(len_bp) / 1e6, 2),
  HH_Mb = round(sum(len_bp[State == "HH"]) / 1e6, 2),
  RH_Mb = round(sum(len_bp[State == "RH"]) / 1e6, 2),
  mean_tract_Mb = round(mean(len_Mb), 3),
  max_tract_Mb = round(max(len_Mb), 3)
), by = Chr]

chr_summary[, chr_len_Mb := chr_len[Chr] / 1e6]
chr_summary[, donor_pct := round(100 * total_Mb / chr_len_Mb, 2)]

setorder(chr_summary, Chr)
print(chr_summary)

cat("\n=== GENOME TOTALS ===\n")
cat("Total donor (HH+RH):", round(sum(tracts$len_bp) / 1e6, 2), "Mb\n")
cat("Total HH:", round(sum(tracts$len_bp[tracts$State == "HH"]) / 1e6, 2), "Mb\n")
cat("Total RH:", round(sum(tracts$len_bp[tracts$State == "RH"]) / 1e6, 2), "Mb\n")
cat("Genome-wide donor %:",
    round(100 * sum(tracts$len_bp) / sum(chr_len), 2), "%\n")

cat("\n=== DONE ===\n")
