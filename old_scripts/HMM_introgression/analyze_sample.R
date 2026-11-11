#!/usr/bin/env Rscript
# ============================================================================
# ANALYZE SINGLE SAMPLE WITH RIGIDITY-BASED INTROGRESSION CALLING
# ============================================================================
# Usage: Rscript analyze_sample.R <input.statepath.tsv.gz> [min_markers]
#
# Example:
#   Rscript analyze_sample.R PN12_SID1120.statepath.tsv.gz 5
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript analyze_sample.R <input.statepath.tsv.gz> [min_markers]")
}

input_file <- args[1]
min_markers <- if (length(args) >= 2) as.integer(args[2]) else 5

if (!file.exists(input_file)) {
  stop("File not found: ", input_file)
}

sample_name <- gsub("\\.statepath\\.tsv\\.gz$", "", basename(input_file))
out_dir <- dirname(input_file)

cat("=== INTROGRESSION ANALYSIS ===\n")
cat("Sample:", sample_name, "\n")
cat("Rigidity (min_markers):", min_markers, "\n\n")

# Chromosome sizes
chr_len <- c(
  chr1 = 308452471, chr2 = 243675191, chr3 = 238017767,
  chr4 = 250330460, chr5 = 226353449, chr6 = 181357234,
  chr7 = 185808916, chr8 = 182411202, chr9 = 163004744,
  chr10 = 152435371
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

get_state_prob <- function(state, p_rr, p_rh, p_hh) {
  ifelse(state == "RR", p_rr, ifelse(state == "RH", p_rh, p_hh))
}

smooth_states <- function(states, probs_matrix, min_mk) {
  if (length(states) == 0 || min_mk <= 1) return(states)

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
      if (r$lengths[i] < min_mk) {
        seg_idx <- starts[i]:ends[i]

        log_probs <- log(probs_matrix[seg_idx, , drop = FALSE] + 1e-10)
        state_scores <- colSums(log_probs)

        left_st <- if (i > 1) r$values[i - 1] else NULL
        right_st <- if (i < n_runs) r$values[i + 1] else NULL
        bonus <- 0.5 * r$lengths[i]

        if (!is.null(left_st)) {
          state_scores[left_st] <- state_scores[left_st] + bonus
        }
        if (!is.null(right_st)) {
          state_scores[right_st] <- state_scores[right_st] + bonus
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

# ============================================================================
# READ AND PROCESS DATA
# ============================================================================

cat("Reading data...\n")
dt <- fread(cmd = paste("gunzip -c", shQuote(input_file)), showProgress = FALSE)
cat("Total SNPs:", format(nrow(dt), big.mark = ","), "\n\n")

cat("=== ORIGINAL STATE COUNTS ===\n")
orig_table <- table(dt$STATE)
print(orig_table)
cat("\n")

# Process each chromosome
tracts_list <- list()
smoothed_counts <- c(RR = 0, RH = 0, HH = 0)

for (chr in paste0("chr", 1:10)) {
  dchr <- dt[CHROM == chr]
  if (nrow(dchr) == 0) next

  orig_states <- dchr$STATE
  probs <- as.matrix(dchr[, .(P_RR, P_RH, P_HH)])
  colnames(probs) <- c("RR", "RH", "HH")

  smooth_st <- smooth_states(orig_states, probs, min_markers)

  # Count smoothed states
  smooth_tab <- table(smooth_st)
  for (s in names(smooth_tab)) {
    smoothed_counts[s] <- smoothed_counts[s] + smooth_tab[s]
  }

  # Extract tracts
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
tracts[, len_Mb := round(len_bp / 1e6, 4)]

# ============================================================================
# OUTPUT RESULTS
# ============================================================================

cat("=== SMOOTHED STATE COUNTS ===\n")
print(smoothed_counts)
cat("\n")

cat("=== INTROGRESSION SUMMARY ===\n")
cat("Total donor tracts:", nrow(tracts), "\n")
cat("  HH tracts:", sum(tracts$State == "HH"), "\n")
cat("  RH tracts:", sum(tracts$State == "RH"), "\n\n")

if (nrow(tracts) > 0) {
  cat("Tract size summary (Mb):\n")
  print(summary(tracts$len_Mb))
  cat("\n")

  cat("Top 15 largest introgression tracts:\n")
  print(tracts[order(-len_bp)][1:min(15, nrow(tracts)),
    .(Chr, Start, End, State, n_SNPs, len_Mb,
      mean_post = round(mean_posterior, 4))])

  # Chromosome summary
  cat("\n=== CHROMOSOME SUMMARY ===\n")
  chr_summary <- tracts[, .(
    n_tracts = .N,
    total_Mb = round(sum(len_bp) / 1e6, 2),
    HH_Mb = round(sum(len_bp[State == "HH"]) / 1e6, 2),
    RH_Mb = round(sum(len_bp[State == "RH"]) / 1e6, 2)
  ), by = Chr]

  chr_summary[, chr_len_Mb := chr_len[Chr] / 1e6]
  chr_summary[, donor_pct := round(100 * total_Mb / chr_len_Mb, 2)]
  setorder(chr_summary, Chr)
  print(chr_summary)

  cat("\n=== GENOME TOTALS ===\n")
  total_donor <- sum(tracts$len_bp)
  cat("Total donor (HH+RH):", round(total_donor / 1e6, 2), "Mb\n")
  cat("Total HH:", round(sum(tracts$len_bp[tracts$State == "HH"]) / 1e6, 2),
      "Mb\n")
  cat("Total RH:", round(sum(tracts$len_bp[tracts$State == "RH"]) / 1e6, 2),
      "Mb\n")
  cat("Genome-wide donor %:", round(100 * total_donor / sum(chr_len), 2), "%\n")

  # Save tracts
  out_file <- file.path(
    out_dir,
    paste0(sample_name, ".tracts.rigidity", min_markers, ".tsv")
  )
  fwrite(tracts, out_file, sep = "\t")
  cat("\nTracts saved to:", out_file, "\n")
} else {
  cat("No donor tracts found.\n")
}

cat("\n=== DONE ===\n")
