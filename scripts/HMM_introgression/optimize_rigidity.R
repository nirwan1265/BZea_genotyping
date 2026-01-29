suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# OPTIMIZE RIGIDITY (min_markers) FOR HMM INTROGRESSION CALLING
# ============================================================================
# This script finds the optimal rigidity value by analyzing:
# 1. Posterior probability confidence of called segments
# 2. Trade-off between false positives (noise) and false negatives (lost signal)
#
# Analogous to RTIGER's optimize_R() but adapted for your HMM output
# ============================================================================

# ---------------- SETTINGS ----------------
state_dir <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/statepaths"
out_dir <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/rigidity_optimization"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Number of samples to use for optimization (use subset for speed)
n_samples_to_test <- 20  # Set to NULL to use all samples

# Range of rigidity values to test
rigidity_values <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 50)

# Chromosome sizes
chr_len <- c(308452471, 243675191, 238017767, 250330460, 226353449,
             181357234, 185808916, 182411202, 163004744, 152435371)
names(chr_len) <- paste0("chr", 1:10)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Get posterior probability for the called state
get_state_prob <- function(state, P_RR, P_RH, P_HH) {
  ifelse(state == "RR", P_RR,
         ifelse(state == "RH", P_RH, P_HH))
}

# Apply smoothing with given min_markers
smooth_states <- function(states, probs_matrix, min_markers) {
  if (length(states) == 0 || min_markers <= 1) return(states)

  changed <- TRUE
  max_iter <- 100
  iter <- 0

  while (changed && iter < max_iter) {
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

        # Sum log-probs for each state
        log_probs <- log(probs_matrix[seg_idx, , drop = FALSE] + 1e-10)
        state_scores <- colSums(log_probs)

        # Add neighbor bonus
        left_state <- if (i > 1) r$values[i - 1] else NULL
        right_state <- if (i < n_runs) r$values[i + 1] else NULL
        neighbor_bonus <- 0.5 * r$lengths[i]

        if (!is.null(left_state)) state_scores[left_state] <- state_scores[left_state] + neighbor_bonus
        if (!is.null(right_state)) state_scores[right_state] <- state_scores[right_state] + neighbor_bonus

        best_state <- names(which.max(state_scores))

        if (best_state != r$values[i]) {
          states[seg_idx] <- best_state
          changed <- TRUE
        }
      }
    }
  }

  return(states)
}

# Calculate metrics for a given rigidity value
calculate_metrics <- function(dt, min_markers) {
  if (nrow(dt) == 0) {
    return(list(
      n_segments_original = 0,
      n_segments_after = 0,
      segments_removed = 0,
      mean_posterior_original = NA,
      mean_posterior_after = NA,
      low_confidence_removed = 0,
      high_confidence_removed = 0,
      donor_bp_original = 0,
      donor_bp_after = 0
    ))
  }

  # Original state analysis
  original_states <- dt$STATE
  probs_matrix <- as.matrix(dt[, .(P_RR, P_RH, P_HH)])
  colnames(probs_matrix) <- c("RR", "RH", "HH")

  # Get original segment info
  r_orig <- rle(original_states)
  n_orig <- length(r_orig$values)

  # Calculate posterior for each SNP's called state
  dt$state_prob <- get_state_prob(dt$STATE, dt$P_RR, dt$P_RH, dt$P_HH)
  mean_post_orig <- mean(dt$state_prob, na.rm = TRUE)

  # Identify short segments that would be removed
  ends_orig <- cumsum(r_orig$lengths)
  starts_orig <- ends_orig - r_orig$lengths + 1L

  short_segments <- which(r_orig$lengths < min_markers)

  # Calculate confidence of segments to be removed
  low_conf_removed <- 0
  high_conf_removed <- 0

  for (s in short_segments) {
    seg_idx <- starts_orig[s]:ends_orig[s]
    seg_mean_prob <- mean(dt$state_prob[seg_idx], na.rm = TRUE)
    if (seg_mean_prob < 0.8) {
      low_conf_removed <- low_conf_removed + 1
    } else {
      high_conf_removed <- high_conf_removed + 1
    }
  }

  # Original donor bp
  donor_idx_orig <- which(original_states %in% c("HH", "RH"))
  donor_bp_orig <- length(donor_idx_orig)

  # Apply smoothing
  smoothed_states <- smooth_states(original_states, probs_matrix, min_markers)

  # Smoothed segment analysis
  r_smooth <- rle(smoothed_states)
  n_smooth <- length(r_smooth$values)

  # Recalculate posterior for smoothed states
  smoothed_prob <- get_state_prob(smoothed_states, dt$P_RR, dt$P_RH, dt$P_HH)
  mean_post_smooth <- mean(smoothed_prob, na.rm = TRUE)

  # Smoothed donor bp
  donor_idx_smooth <- which(smoothed_states %in% c("HH", "RH"))
  donor_bp_smooth <- length(donor_idx_smooth)

  return(list(
    n_segments_original = n_orig,
    n_segments_after = n_smooth,
    segments_removed = n_orig - n_smooth,
    mean_posterior_original = mean_post_orig,
    mean_posterior_after = mean_post_smooth,
    low_confidence_removed = low_conf_removed,
    high_confidence_removed = high_conf_removed,
    donor_snps_original = donor_bp_orig,
    donor_snps_after = donor_bp_smooth
  ))
}

# ============================================================================
# MAIN OPTIMIZATION
# ============================================================================

cat("=== RIGIDITY OPTIMIZATION FOR HMM INTROGRESSION ===\n\n")

# Get sample files
files <- sort(list.files(state_dir, pattern = "\\.statepath\\.tsv\\.gz$", full.names = TRUE))
if (length(files) == 0) stop("No *.statepath.tsv.gz in state_dir")

cat("Found", length(files), "samples\n")

# Subset samples if requested
if (!is.null(n_samples_to_test) && n_samples_to_test < length(files)) {
  set.seed(42)
  files <- sample(files, n_samples_to_test)
  cat("Using random subset of", n_samples_to_test, "samples for optimization\n")
}

cat("Testing rigidity values:", paste(rigidity_values, collapse = ", "), "\n\n")

# Initialize results storage
results_list <- vector("list", length(rigidity_values))

for (r_idx in seq_along(rigidity_values)) {
  min_markers <- rigidity_values[r_idx]
  cat("Testing min_markers =", min_markers, "... ")

  # Aggregate metrics across all samples
  all_metrics <- list(
    n_segments_original = 0,
    n_segments_after = 0,
    segments_removed = 0,
    mean_posterior_original = c(),
    mean_posterior_after = c(),
    low_confidence_removed = 0,
    high_confidence_removed = 0,
    donor_snps_original = 0,
    donor_snps_after = 0
  )

  n_processed <- 0

  for (f in files) {
    dt <- tryCatch(
      fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE),
      error = function(e) NULL
    )

    if (is.null(dt) || nrow(dt) == 0) next

    # Standardize column names
    if ("CHROM" %in% names(dt)) {
      # ok
    } else if ("Chr" %in% names(dt)) {
      setnames(dt, "Chr", "CHROM")
    } else next

    if (!"STATE" %in% names(dt)) {
      if ("State" %in% names(dt)) setnames(dt, "State", "STATE")
      else next
    }

    # Check for probability columns
    if (!all(c("P_RR", "P_RH", "P_HH") %in% names(dt))) next

    dt <- dt[CHROM %in% names(chr_len)]

    # Process each chromosome separately
    for (chr in unique(dt$CHROM)) {
      dchr <- dt[CHROM == chr]
      if (nrow(dchr) < 10) next

      metrics <- calculate_metrics(dchr, min_markers)

      all_metrics$n_segments_original <- all_metrics$n_segments_original + metrics$n_segments_original
      all_metrics$n_segments_after <- all_metrics$n_segments_after + metrics$n_segments_after
      all_metrics$segments_removed <- all_metrics$segments_removed + metrics$segments_removed
      all_metrics$mean_posterior_original <- c(all_metrics$mean_posterior_original, metrics$mean_posterior_original)
      all_metrics$mean_posterior_after <- c(all_metrics$mean_posterior_after, metrics$mean_posterior_after)
      all_metrics$low_confidence_removed <- all_metrics$low_confidence_removed + metrics$low_confidence_removed
      all_metrics$high_confidence_removed <- all_metrics$high_confidence_removed + metrics$high_confidence_removed
      all_metrics$donor_snps_original <- all_metrics$donor_snps_original + metrics$donor_snps_original
      all_metrics$donor_snps_after <- all_metrics$donor_snps_after + metrics$donor_snps_after
    }

    n_processed <- n_processed + 1
  }

  cat("processed", n_processed, "samples\n")

  # Summarize results
  results_list[[r_idx]] <- data.table(
    min_markers = min_markers,
    n_samples = n_processed,
    total_segments_original = all_metrics$n_segments_original,
    total_segments_after = all_metrics$n_segments_after,
    segments_removed = all_metrics$segments_removed,
    pct_segments_removed = 100 * all_metrics$segments_removed / max(1, all_metrics$n_segments_original),
    mean_posterior_original = mean(all_metrics$mean_posterior_original, na.rm = TRUE),
    mean_posterior_after = mean(all_metrics$mean_posterior_after, na.rm = TRUE),
    posterior_improvement = mean(all_metrics$mean_posterior_after, na.rm = TRUE) -
                            mean(all_metrics$mean_posterior_original, na.rm = TRUE),
    low_conf_segments_removed = all_metrics$low_confidence_removed,
    high_conf_segments_removed = all_metrics$high_confidence_removed,
    # Ratio of good removals to bad removals (higher is better)
    removal_quality_ratio = all_metrics$low_confidence_removed /
                            max(1, all_metrics$high_confidence_removed),
    donor_snps_original = all_metrics$donor_snps_original,
    donor_snps_after = all_metrics$donor_snps_after,
    donor_snps_lost = all_metrics$donor_snps_original - all_metrics$donor_snps_after,
    pct_donor_lost = 100 * (all_metrics$donor_snps_original - all_metrics$donor_snps_after) /
                     max(1, all_metrics$donor_snps_original)
  )
}

# Combine results
results <- rbindlist(results_list)

# ============================================================================
# SCORING AND RECOMMENDATION
# ============================================================================

# Calculate a composite score:
# - Maximize posterior improvement (removing low-confidence calls)
# - Maximize removal quality ratio (remove bad segments, keep good ones)
# - Minimize loss of donor signal

# Normalize each metric to 0-1 scale
normalize <- function(x, higher_is_better = TRUE) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0.5, length(x)))
  norm <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  if (!higher_is_better) norm <- 1 - norm
  return(norm)
}

results[, score_posterior := normalize(posterior_improvement, higher_is_better = TRUE)]
results[, score_quality := normalize(removal_quality_ratio, higher_is_better = TRUE)]
results[, score_donor_preserved := normalize(pct_donor_lost, higher_is_better = FALSE)]

# Composite score (weighted average)
# Weight: posterior improvement (30%), removal quality (40%), donor preservation (30%)
results[, composite_score := 0.30 * score_posterior +
                             0.40 * score_quality +
                             0.30 * score_donor_preserved]

# Find optimal
optimal_idx <- which.max(results$composite_score)
optimal_rigidity <- results$min_markers[optimal_idx]

# ============================================================================
# OUTPUT RESULTS
# ============================================================================

cat("\n=== OPTIMIZATION RESULTS ===\n\n")

print(results[, .(min_markers,
                  pct_segments_removed = round(pct_segments_removed, 1),
                  posterior_improvement = round(posterior_improvement * 100, 2),
                  removal_quality_ratio = round(removal_quality_ratio, 2),
                  pct_donor_lost = round(pct_donor_lost, 2),
                  composite_score = round(composite_score, 3))])

cat("\n=== RECOMMENDATION ===\n")
cat("Optimal min_markers (rigidity):", optimal_rigidity, "\n\n")

cat("Interpretation:\n")
cat("- posterior_improvement: Higher = better (removing low-confidence noise)\n")
cat("- removal_quality_ratio: Higher = better (removing bad segments, not good ones)\n")
cat("- pct_donor_lost: Lower = better (preserving real introgression signal)\n")
cat("- composite_score: Weighted combination (higher = better)\n\n")

# Save results
out_file <- file.path(out_dir, "rigidity_optimization_results.tsv")
fwrite(results, out_file, sep = "\t")
cat("Results saved to:", out_file, "\n")

# ============================================================================
# PLOT RESULTS
# ============================================================================

# Save plot
plot_file <- file.path(out_dir, "rigidity_optimization_plot.pdf")
pdf(plot_file, width = 10, height = 8)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Plot 1: Posterior improvement
plot(results$min_markers, results$posterior_improvement * 100,
     type = "b", pch = 19, col = "blue",
     xlab = "min_markers (rigidity)",
     ylab = "Posterior Improvement (%)",
     main = "A) Confidence Improvement",
     log = "x")
abline(v = optimal_rigidity, lty = 2, col = "red")

# Plot 2: Removal quality ratio
plot(results$min_markers, results$removal_quality_ratio,
     type = "b", pch = 19, col = "darkgreen",
     xlab = "min_markers (rigidity)",
     ylab = "Removal Quality Ratio\n(low-conf / high-conf removed)",
     main = "B) Quality of Segment Removal",
     log = "x")
abline(v = optimal_rigidity, lty = 2, col = "red")

# Plot 3: Donor signal loss
plot(results$min_markers, results$pct_donor_lost,
     type = "b", pch = 19, col = "orange",
     xlab = "min_markers (rigidity)",
     ylab = "% Donor SNPs Lost",
     main = "C) Donor Signal Loss",
     log = "x")
abline(v = optimal_rigidity, lty = 2, col = "red")

# Plot 4: Composite score
plot(results$min_markers, results$composite_score,
     type = "b", pch = 19, col = "purple", lwd = 2,
     xlab = "min_markers (rigidity)",
     ylab = "Composite Score",
     main = paste0("D) Composite Score (Optimal: ", optimal_rigidity, ")"),
     log = "x")
abline(v = optimal_rigidity, lty = 2, col = "red")
points(optimal_rigidity, results$composite_score[optimal_idx],
       pch = 19, col = "red", cex = 2)

dev.off()

cat("Plot saved to:", plot_file, "\n")

cat("\n=== DONE ===\n")
cat("Use min_markers =", optimal_rigidity, "in introgression_stats_rigidity.R\n")
