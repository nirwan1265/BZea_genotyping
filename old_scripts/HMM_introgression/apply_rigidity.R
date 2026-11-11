#!/usr/bin/env Rscript
# ============================================================================
# APPLY RIGIDITY TO EXISTING STATEPATH FILES
# ============================================================================
# This script applies RTIGER-style rigidity to already-computed statepath files.
# It uses a sliding window to smooth states based on posterior probabilities.
#
# Usage: Rscript apply_rigidity.R <input.statepath.tsv.gz> <rigidity> [min_block_kb]
#
# Arguments:
#   input: statepath.tsv.gz file (from your HMM)
#   rigidity: Window size in SNPs (try 50-200 for BC2S3)
#   min_block_kb: Minimum block size in kb (default: 100)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript apply_rigidity.R <input.statepath.tsv.gz> <rigidity> [min_block_kb]\n")
  cat("\nArguments:\n")
  cat("  input:       statepath.tsv.gz file\n")
  cat("  rigidity:    Window size in SNPs (50-200 for BC2S3)\n")
  cat("  min_block_kb: Min block size to keep (default: 100 kb)\n")
  stop()
}

input_file <- args[1]
rigidity <- as.integer(args[2])
min_block_bp <- if (length(args) >= 3) as.numeric(args[3]) * 1000 else 100000

if (!file.exists(input_file)) stop("File not found: ", input_file)

sample_name <- gsub("\\.statepath\\.tsv\\.gz$", "", basename(input_file))
out_dir <- dirname(input_file)

cat("=== APPLY RIGIDITY ===\n")
cat("Sample:", sample_name, "\n")
cat("Rigidity (R):", rigidity, "SNPs\n")
cat("Min block size:", min_block_bp / 1000, "kb\n\n")

chr_len <- c(
  chr1 = 308452471, chr2 = 243675191, chr3 = 238017767,
  chr4 = 250330460, chr5 = 226353449, chr6 = 181357234,
  chr7 = 185808916, chr8 = 182411202, chr9 = 163004744,
  chr10 = 152435371
)

state_to_geno <- c(RR = "AA", RH = "AB", HH = "BB")
geno_colors <- c(AA = "#3366CC", AB = "#663399", BB = "#CC3333")

# ============================================================================
# RIGIDITY SMOOTHING FUNCTION (VECTORIZED - FAST)
# ============================================================================
# Uses data.table frollsum for fast sliding window smoothing

apply_rigidity_smooth_fast <- function(dt, R) {
  # Use frollsum for vectorized rolling sums (much faster than for-loop)
  # align="center" for centered window

  dt[, roll_RR := frollsum(P_RR, n = R, align = "center", na.rm = TRUE, fill = NA)]
  dt[, roll_RH := frollsum(P_RH, n = R, align = "center", na.rm = TRUE, fill = NA)]
  dt[, roll_HH := frollsum(P_HH, n = R, align = "center", na.rm = TRUE, fill = NA)]

  # Handle NA at edges (use smaller window)
  half_R <- R %/% 2
  n <- nrow(dt)

  # Fill edges with cumsum for start
  for (i in 1:min(half_R, n)) {
    window_end <- min(i + half_R, n)
    dt$roll_RR[i] <- sum(dt$P_RR[1:window_end], na.rm = TRUE)
    dt$roll_RH[i] <- sum(dt$P_RH[1:window_end], na.rm = TRUE)
    dt$roll_HH[i] <- sum(dt$P_HH[1:window_end], na.rm = TRUE)
  }

  # Fill edges at end
  for (i in max(n - half_R + 1, 1):n) {
    window_start <- max(i - half_R, 1)
    dt$roll_RR[i] <- sum(dt$P_RR[window_start:n], na.rm = TRUE)
    dt$roll_RH[i] <- sum(dt$P_RH[window_start:n], na.rm = TRUE)
    dt$roll_HH[i] <- sum(dt$P_HH[window_start:n], na.rm = TRUE)
  }

  # Assign state with highest rolling sum
  smoothed_state <- ifelse(
    dt$roll_RR >= dt$roll_RH & dt$roll_RR >= dt$roll_HH, "RR",
    ifelse(dt$roll_RH >= dt$roll_HH, "RH", "HH")
  )

  # Clean up temp columns
  dt[, c("roll_RR", "roll_RH", "roll_HH") := NULL]

  return(smoothed_state)
}

# ============================================================================
# READ DATA
# ============================================================================

cat("Reading data...\n")
dt <- fread(cmd = paste("gunzip -c", shQuote(input_file)), showProgress = FALSE)
cat("Total SNPs:", format(nrow(dt), big.mark = ","), "\n")

# Original state counts
cat("\nOriginal states:\n")
print(table(dt$STATE))

# Original segment count
orig_segs <- sum(diff(as.numeric(factor(dt$STATE))) != 0) + 1
cat("Original segments:", orig_segs, "\n")

# ============================================================================
# APPLY RIGIDITY PER CHROMOSOME
# ============================================================================

cat("\nApplying rigidity (R =", rigidity, ")...\n")

dt[, SMOOTH_STATE := "RR"]  # Initialize

for (chr in paste0("chr", 1:10)) {
  idx <- which(dt$CHROM == chr)
  if (length(idx) == 0) next

  chr_data <- dt[idx]
  smoothed <- apply_rigidity_smooth_fast(chr_data, rigidity)
  dt$SMOOTH_STATE[idx] <- smoothed
}

# Smoothed state counts
cat("\nSmoothed states:\n")
print(table(dt$SMOOTH_STATE))

# Smoothed segment count
smooth_segs <- sum(diff(as.numeric(factor(dt$SMOOTH_STATE))) != 0) + 1
cat("Smoothed segments:", smooth_segs, "\n")
cat("Reduction:", round(100 * (1 - smooth_segs / orig_segs), 1), "%\n")

# ============================================================================
# CREATE COMPLETE BLOCKS
# ============================================================================

cat("\nCreating complete blocks...\n")

all_blocks <- list()

for (chr in paste0("chr", 1:10)) {
  chr_data <- dt[CHROM == chr]
  if (nrow(chr_data) == 0) {
    # No data - fill with AA
    all_blocks[[chr]] <- data.table(
      chrom = chr,
      start = 1,
      end = chr_len[chr],
      genotype = "AA",
      n_snps = 0
    )
    next
  }

  # Get runs of smoothed state
  r <- rle(chr_data$SMOOTH_STATE)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  chr_blocks <- list()

  for (i in seq_along(r$values)) {
    geno <- state_to_geno[r$values[i]]
    start_pos <- chr_data$POS[starts[i]]
    end_pos <- chr_data$POS[ends[i]]

    chr_blocks[[i]] <- data.table(
      chrom = chr,
      start = start_pos,
      end = end_pos,
      genotype = geno,
      n_snps = r$lengths[i]
    )
  }

  blocks <- rbindlist(chr_blocks)
  blocks[, len_bp := end - start]

  # Remove NA genotypes (shouldn't happen but just in case)
  blocks <- blocks[!is.na(genotype)]

  # Filter small blocks (keep all non-AA blocks regardless of size)
  blocks <- blocks[len_bp >= min_block_bp | genotype != "AA"]

  # Fill gaps: extend to chromosome boundaries
  if (nrow(blocks) > 0) {
    setorder(blocks, start)

    # Extend first block to start of chromosome
    blocks$start[1] <- 1

    # Extend last block to end of chromosome
    blocks$end[nrow(blocks)] <- chr_len[chr]

    # Fill gaps between blocks
    if (nrow(blocks) > 1) {
      for (j in 2:nrow(blocks)) {
        if (blocks$start[j] > blocks$end[j - 1] + 1) {
          # Gap exists - extend previous block
          blocks$end[j - 1] <- blocks$start[j] - 1
        }
      }
    }

    # Merge adjacent same-state blocks
    if (nrow(blocks) == 1) {
      # Only one block, nothing to merge
    } else {
      merged <- list()
      current <- blocks[1]

      for (j in 2:nrow(blocks)) {
        g_current <- current$genotype
        g_next <- blocks$genotype[j]
        # Safe comparison handling NAs
        if (!is.na(g_current) && !is.na(g_next) && g_next == g_current) {
          current$end <- blocks$end[j]
          current$n_snps <- current$n_snps + blocks$n_snps[j]
        } else {
          merged[[length(merged) + 1]] <- current
          current <- blocks[j]
        }
      }
      merged[[length(merged) + 1]] <- current

      blocks <- rbindlist(merged)
    }
  } else {
    blocks <- data.table(
      chrom = chr,
      start = 1,
      end = chr_len[chr],
      genotype = "AA",
      n_snps = 0
    )
  }

  blocks[, len_bp := end - start]
  all_blocks[[chr]] <- blocks
}

final_blocks <- rbindlist(all_blocks)
final_blocks[, len_bp := end - start]

cat("Final blocks:", nrow(final_blocks), "\n")

# ============================================================================
# OUTPUT
# ============================================================================

# BED file
bed_file <- file.path(out_dir, paste0(sample_name, ".rigidity", rigidity, ".bed"))
bed_out <- final_blocks[, .(chrom, start, end, genotype)]
fwrite(bed_out, bed_file, sep = "\t", col.names = FALSE)
cat("\nBED file:", bed_file, "\n")

# Show blocks
cat("\nComplete blocks:\n")
print(bed_out)

# PNG visualization
png_file <- file.path(out_dir, paste0(sample_name, ".rigidity", rigidity, ".png"))

n_chr <- 10
plot_height <- 150 + n_chr * 80

png(png_file, width = 1400, height = plot_height, res = 100)

par(mar = c(4, 8, 4, 2))

plot(NULL, xlim = c(0, max(chr_len)), ylim = c(0, n_chr + 1),
     xlab = "Position (Mb)", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n",
     main = paste0(sample_name, " - Rigidity R=", rigidity))

axis_pos <- seq(0, 300e6, by = 50e6)
axis(1, at = axis_pos, labels = axis_pos / 1e6)

for (i in 1:n_chr) {
  chr <- paste0("chr", i)
  y_pos <- n_chr - i + 1

  rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35,
       col = "gray90", border = NA)

  chr_blocks <- final_blocks[chrom == chr]

  for (j in seq_len(nrow(chr_blocks))) {
    blk <- chr_blocks[j]
    rect(blk$start, y_pos - 0.35, blk$end, y_pos + 0.35,
         col = geno_colors[blk$genotype], border = NA)
  }

  rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35,
       col = NA, border = "gray40", lwd = 0.5)

  text(-5e6, y_pos, chr, adj = 1, cex = 0.9, font = 2)
}

legend("topright",
       legend = c("AA (B73/B73)", "AB (Het)", "BB (Teo/Teo)"),
       fill = geno_colors[c("AA", "AB", "BB")],
       border = "gray40", cex = 0.9, bg = "white")

dev.off()

cat("PNG file:", png_file, "\n")

# Summary
cat("\n=== SUMMARY ===\n")
summary_dt <- final_blocks[, .(
  n_blocks = .N,
  total_Mb = round(sum(len_bp) / 1e6, 2)
), by = genotype]
print(summary_dt)

cat("\nChromosome donor %:\n")
chr_summary <- final_blocks[, .(
  donor_Mb = round(sum(len_bp[genotype != "AA"]) / 1e6, 2)
), by = chrom]
chr_summary[, chr_len_Mb := chr_len[chrom] / 1e6]
chr_summary[, donor_pct := round(100 * donor_Mb / chr_len_Mb, 1)]
setorder(chr_summary, chrom)
print(chr_summary)

cat("\n=== DONE ===\n")
