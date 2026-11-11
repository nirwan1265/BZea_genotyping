suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# INTROGRESSION STATS WITH MARKER-COUNT RIGIDITY (RTIGER-style)
# ============================================================================
# This version uses marker-count based smoothing instead of cM-based bridging.
# Two approaches implemented:
#   1. Simple rigidity: flip short segments (< min_markers) to neighbor state
#   2. Probability-based: use posterior probs to refine boundaries
# ============================================================================

# ---------------- SETTINGS ----------------
state_dir <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/statepaths"
map_dir   <- "/rsstu/users/r/rrellan/sara/ref/NAM_genetic_map"

out_dir   <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/introgression_rigidity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============== KEY PARAMETERS ==============
# Minimum markers to keep a segment (RTIGER-style rigidity)
min_markers <- 5  # Try 3, 5, 10 - segments with fewer markers get smoothed

# Method for smoothing: "simple" or "probability"
# - "simple": flip short segments to most common neighbor state
# - "probability": use posterior probs to decide flips
smoothing_method <- "probability"

# Define donor states
donor_set <- c("HH", "RH")

# Output files
out_summary <- file.path(out_dir, "all_samples.ref_alt_het_perc.rigidity.tsv")
out_all_tracts <- file.path(out_dir, "all_samples.alt_het.tracts.rigidity.tsv.gz")

# Chromosome sizes (bp)
chr_len <- c(308452471, 243675191, 238017767, 250330460, 226353449,
             181357234, 185808916, 182411202, 163004744, 152435371)
names(chr_len) <- paste0("chr", 1:10)
genome_bp <- sum(chr_len)

# ---------------- MAPS (LOAD ONCE) ----------------
maps <- lapply(paste0("chr", 1:10), function(chr) {
  f <- file.path(map_dir, paste0(chr, ".map"))
  m <- fread(f, header = FALSE)
  setnames(m, c("bp", "col2", "cM"))
  setorder(m, bp)
  m
})
names(maps) <- paste0("chr", 1:10)

interp_cm_vec <- function(chr, pos) {
  m <- maps[[chr]]
  approx(x = m$bp, y = m$cM, xout = pos, rule = 2)$y
}

# ============================================================================
# METHOD 1: SIMPLE MARKER-COUNT SMOOTHING
# ============================================================================
# If a segment has < min_markers SNPs, flip it to the dominant neighbor state
# This mimics RTIGER's rigidity concept applied post-hoc

smooth_by_marker_count <- function(states, min_markers) {
  if (length(states) == 0) return(states)

  # Find runs
  r <- rle(states)
  n_runs <- length(r$values)

  if (n_runs <= 1) return(states)

  # Iteratively smooth until no short segments remain
  changed <- TRUE
  max_iter <- 100
  iter <- 0

  while (changed && iter < max_iter) {
    changed <- FALSE
    iter <- iter + 1

    r <- rle(states)
    n_runs <- length(r$values)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L

    for (i in seq_len(n_runs)) {
      if (r$lengths[i] < min_markers) {
        # This segment is too short - flip it
        # Decide what to flip to: use neighbor states
        left_state <- if (i > 1) r$values[i - 1] else NA
        right_state <- if (i < n_runs) r$values[i + 1] else NA

        # If both neighbors same, use that
        # Otherwise, use left neighbor (or right if left is NA)
        new_state <- if (!is.na(left_state) && !is.na(right_state) && left_state == right_state) {
          left_state
        } else if (!is.na(left_state)) {
          left_state
        } else if (!is.na(right_state)) {
          right_state
        } else {
          r$values[i]  # Keep original if no neighbors
        }

        if (new_state != r$values[i]) {
          states[starts[i]:ends[i]] <- new_state
          changed <- TRUE
        }
      }
    }
  }

  return(states)
}

# ============================================================================
# METHOD 2: PROBABILITY-BASED SMOOTHING (RTIGER-style findsplit)
# ============================================================================
# Uses posterior probabilities to refine boundaries and smooth short segments

smooth_by_probability <- function(dt, min_markers) {
  # dt must have: STATE, P_RR, P_RH, P_HH
  if (nrow(dt) == 0) return(dt)

  # Check if probability columns exist
  has_probs <- all(c("P_RR", "P_RH", "P_HH") %in% names(dt))
  if (!has_probs) {
    # Fall back to simple method
    dt$STATE <- smooth_by_marker_count(dt$STATE, min_markers)
    return(dt)
  }

  states <- dt$STATE
  probs <- as.matrix(dt[, .(P_RR, P_RH, P_HH)])
  colnames(probs) <- c("RR", "RH", "HH")

  # Convert to log-probs for numerical stability (add small epsilon)
  log_probs <- log(probs + 1e-10)

  # Iteratively smooth short segments
  changed <- TRUE
  max_iter <- 100
  iter <- 0

  while (changed && iter < max_iter) {
    changed <- FALSE
    iter <- iter + 1

    r <- rle(states)
    n_runs <- length(r$values)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L

    for (i in seq_len(n_runs)) {
      if (r$lengths[i] < min_markers && n_runs > 1) {
        # Short segment - evaluate using probabilities
        seg_idx <- starts[i]:ends[i]

        # Sum of log-probs for each possible state across this segment
        state_scores <- colSums(log_probs[seg_idx, , drop = FALSE])

        # Also consider context from neighbors (RTIGER's findsplit idea)
        # Weight current evidence + neighbor preference
        left_state <- if (i > 1) r$values[i - 1] else NULL
        right_state <- if (i < n_runs) r$values[i + 1] else NULL

        # Bonus for matching neighbors (encourages continuity)
        neighbor_bonus <- 0.5 * r$lengths[i]  # Scale with segment length

        if (!is.null(left_state)) {
          state_scores[left_state] <- state_scores[left_state] + neighbor_bonus
        }
        if (!is.null(right_state)) {
          state_scores[right_state] <- state_scores[right_state] + neighbor_bonus
        }

        # Pick best state
        best_state <- names(which.max(state_scores))

        if (best_state != r$values[i]) {
          states[seg_idx] <- best_state
          changed <- TRUE
        }
      }
    }
  }

  dt$STATE <- states
  return(dt)
}

# ============================================================================
# FINDSPLIT: RTIGER-style boundary refinement
# ============================================================================
# For each boundary between segments, find optimal split point using probs

refine_boundaries <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  if (!all(c("P_RR", "P_RH", "P_HH") %in% names(dt))) return(dt)

  states <- dt$STATE
  probs <- as.matrix(dt[, .(P_RR, P_RH, P_HH)])
  colnames(probs) <- c("RR", "RH", "HH")
  log_probs <- log(probs + 1e-10)

  r <- rle(states)
  n_runs <- length(r$values)

  if (n_runs <= 1) return(dt)

  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  # For each boundary, find optimal split
  for (i in 1:(n_runs - 1)) {
    left_state <- r$values[i]
    right_state <- r$values[i + 1]

    if (left_state == right_state) next

    # Window around boundary: 5 SNPs on each side (or available)
    window_size <- 5
    boundary_pos <- ends[i]

    win_start <- max(starts[i], boundary_pos - window_size + 1)
    win_end <- min(ends[i + 1], boundary_pos + window_size)

    if (win_end - win_start < 2) next

    win_idx <- win_start:win_end

    # RTIGER's findsplit: find position that maximizes
    # cumsum(left_state probs from left) + cumsum(right_state probs from right)
    left_probs <- log_probs[win_idx, left_state]
    right_probs <- log_probs[win_idx, right_state]

    left_cumsum <- cumsum(left_probs)
    right_cumsum <- rev(cumsum(rev(right_probs)))

    # Total score for each possible split point
    # Split at position k means: left state for 1:k, right state for (k+1):end
    scores <- c(0, left_cumsum[-length(left_cumsum)]) +
              c(right_cumsum[-1], 0)

    best_split <- which.max(scores)
    new_boundary <- win_start + best_split - 1

    # Update states based on new boundary
    if (new_boundary != boundary_pos && new_boundary >= starts[i] && new_boundary < ends[i + 1]) {
      states[starts[i]:new_boundary] <- left_state
      states[(new_boundary + 1):ends[i + 1]] <- right_state
    }
  }

  dt$STATE <- states
  return(dt)
}

# ---------------- RUN SEGMENTS ----------------
make_runs <- function(state, w_cM, pos, chr) {
  r <- rle(state)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  seg <- data.table(
    run_id = seq_along(r$values),
    state = r$values,
    start_i = starts,
    end_i = ends,
    n_sites = r$lengths
  )

  seg[, Chr := chr[start_i]]
  seg[, Start := pos[start_i]]
  seg[, End := pos[end_i]]

  seg[, len_cM := vapply(seq_len(.N), function(i) {
    sum(w_cM[start_i[i]:end_i[i]], na.rm = TRUE)
  }, numeric(1))]

  seg[, len_bp := as.numeric(End - Start + 1L)]
  seg[]
}

# ---------------- EXTRACT DONOR TRACTS ----------------
extract_donor_tracts <- function(seg, donor_set) {
  if (nrow(seg) == 0) return(data.table())

  # Filter to donor segments only, then merge adjacent
  donor_seg <- seg[state %in% donor_set]

  if (nrow(donor_seg) == 0) return(data.table())

  # Merge adjacent donor segments (they may have different states like HH vs RH)
  out <- list()
  k <- 0L
  i <- 1L

  while (i <= nrow(donor_seg)) {
    chr <- donor_seg$Chr[i]
    j <- i

    # Extend through consecutive donor segments on same chromosome
    while (j + 1L <= nrow(donor_seg) &&
           donor_seg$Chr[j + 1L] == chr &&
           donor_seg$run_id[j + 1L] == donor_seg$run_id[j] + 1L) {
      j <- j + 1L
    }

    k <- k + 1L
    tract_runs <- donor_seg[i:j]
    region <- if (any(tract_runs$state == "HH")) "Alt" else "Het"

    out[[k]] <- data.table(
      Chr = chr,
      Start = min(tract_runs$Start),
      End = max(tract_runs$End),
      Region = region,
      n_sites = sum(tract_runs$n_sites),
      len_bp = as.numeric(max(tract_runs$End) - min(tract_runs$Start) + 1L),
      len_cM = sum(tract_runs$len_cM)
    )

    i <- j + 1L
  }

  rbindlist(out, fill = TRUE)
}

# gzip writer
write_tsv_gz <- function(dt, out_gz) {
  tmp <- tempfile(fileext = ".tsv")
  fwrite(dt, tmp, sep = "\t")
  system2("gzip", c("-c", tmp), stdout = out_gz)
  unlink(tmp)
}

# ---------------- PER SAMPLE ----------------
files <- sort(list.files(state_dir, pattern = "\\.statepath\\.tsv\\.gz$", full.names = TRUE))
if (length(files) == 0) stop("No *.statepath.tsv.gz in state_dir")

all_tracts <- vector("list", length(files))
all_stats <- vector("list", length(files))

cat("Using smoothing method:", smoothing_method, "\n")
cat("Minimum markers per segment:", min_markers, "\n\n")

for (ii in seq_along(files)) {
  f <- files[ii]
  sample <- sub("\\.statepath\\.tsv\\.gz$", "", basename(f))

  dt <- tryCatch(
    fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE),
    error = function(e) NULL
  )

  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) == 0) {
    message("SKIP empty statepath: ", basename(f))

    all_stats[[ii]] <- data.table(
      Sample = sample,
      Ref_perc = 100, Alt_perc = 0, Het_perc = 0,
      Ref_bp = genome_bp, Alt_bp = 0, Het_bp = 0
    )

    all_tracts[[ii]] <- data.table(
      Sample = sample,
      Chr = character(), Start = integer(), End = integer(), Region = character(),
      n_sites = integer(), len_bp = numeric(), len_cM = numeric()
    )

    next
  }

  # Standardize column names
  nm <- names(dt)

  if ("CHROM" %in% nm) {
    # ok
  } else if ("Chr" %in% nm) {
    setnames(dt, "Chr", "CHROM")
  } else {
    message("SKIP missing CHROM col: ", basename(f))
    next
  }

  if (!"POS" %in% names(dt)) {
    if ("Pos" %in% names(dt)) setnames(dt, "Pos", "POS")
    else { message("SKIP missing POS col: ", basename(f)); next }
  }

  if (!"STATE" %in% names(dt)) {
    if ("State" %in% names(dt)) setnames(dt, "State", "STATE")
    else { message("SKIP missing STATE col: ", basename(f)); next }
  }

  dt <- dt[CHROM %in% names(chr_len)]
  setorder(dt, CHROM, POS)

  # ============== APPLY SMOOTHING ==============
  tr_list <- list()

  for (chr in names(chr_len)) {
    dchr <- dt[CHROM == chr]
    if (nrow(dchr) == 0) next

    # Apply smoothing based on method
    if (smoothing_method == "probability") {
      dchr <- smooth_by_probability(dchr, min_markers)
      dchr <- refine_boundaries(dchr)
    } else {
      dchr$STATE <- smooth_by_marker_count(dchr$STATE, min_markers)
    }

    # Interpolate cM
    dchr[, cM := interp_cm_vec(chr, POS)]
    dchr[, next_cM := shift(cM, type = "lead")]
    dchr[, w_cM := pmax(0, next_cM - cM)]
    dchr[is.na(w_cM), w_cM := 0]

    # Make runs from smoothed states
    seg <- make_runs(dchr$STATE, dchr$w_cM, dchr$POS, dchr$CHROM)

    # Extract donor tracts
    tr_chr <- extract_donor_tracts(seg, donor_set)

    if (nrow(tr_chr) > 0) tr_list[[chr]] <- tr_chr
  }

  tr <- rbindlist(tr_list, fill = TRUE)
  if (nrow(tr) == 0) {
    tr <- data.table(Chr = character(), Start = integer(), End = integer(), Region = character(),
                     n_sites = integer(), len_bp = numeric(), len_cM = numeric())
  }

  tr[, Sample := sample]
  setcolorder(tr, c("Sample", "Chr", "Start", "End", "Region", "n_sites", "len_bp", "len_cM"))

  # Per-sample output
  out_f <- file.path(out_dir, paste0(sample, ".alt_het.tracts.rigidity.tsv.gz"))
  write_tsv_gz(tr, out_f)

  # Percentages
  Alt_bp <- sum(tr$len_bp[tr$Region == "Alt"], na.rm = TRUE)
  Het_bp <- sum(tr$len_bp[tr$Region == "Het"], na.rm = TRUE)
  Ref_bp <- genome_bp - Alt_bp - Het_bp
  if (Ref_bp < 0) Ref_bp <- 0

  st <- data.table(
    Sample = sample,
    Ref_perc = 100 * Ref_bp / genome_bp,
    Alt_perc = 100 * Alt_bp / genome_bp,
    Het_perc = 100 * Het_bp / genome_bp,
    Ref_bp = Ref_bp,
    Alt_bp = Alt_bp,
    Het_bp = Het_bp
  )

  all_tracts[[ii]] <- tr
  all_stats[[ii]] <- st

  if (ii %% 50 == 0) message("Done ", ii, " / ", length(files))
}

tract_all <- rbindlist(all_tracts, fill = TRUE)
stats_all <- rbindlist(all_stats, fill = TRUE)

# TOTAL row
nS <- nrow(stats_all)
tot_alt <- sum(stats_all$Alt_bp)
tot_het <- sum(stats_all$Het_bp)
tot_ref <- (genome_bp * nS) - tot_alt - tot_het
if (tot_ref < 0) tot_ref <- 0

total_row <- data.table(
  Sample = "TOTAL",
  Ref_perc = 100 * tot_ref / (genome_bp * nS),
  Alt_perc = 100 * tot_alt / (genome_bp * nS),
  Het_perc = 100 * tot_het / (genome_bp * nS),
  Ref_bp = tot_ref, Alt_bp = tot_alt, Het_bp = tot_het
)

stats_out <- rbind(stats_all, total_row, fill = TRUE)
fwrite(stats_out, out_summary, sep = "\t")

write_tsv_gz(tract_all, out_all_tracts)

cat("\n========================================\n")
cat("Wrote:\n")
cat("  Tracts per sample: ", out_dir, "\n", sep = "")
cat("  Combined tracts:   ", out_all_tracts, "\n", sep = "")
cat("  Summary perc:      ", out_summary, "\n", sep = "")
cat("========================================\n")
cat("\nSettings used:\n")
cat("  smoothing_method: ", smoothing_method, "\n")
cat("  min_markers:      ", min_markers, "\n")
cat("  donor_set:        ", paste(donor_set, collapse = ", "), "\n")
