#!/usr/bin/env Rscript
# ============================================================================
# PLOT INTROGRESSION - BED OUTPUT + VISUALIZATION
# ============================================================================
# Outputs:
#   1. BED file with continuous segments (AA/AB/BB format)
#   2. PNG visualization like RTIGER genotype plot
#
# Usage: Rscript plot_introgression.R <input.statepath.tsv.gz> [min_markers]
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_introgression.R <input.statepath.tsv.gz> [min_markers]")
}

input_file <- args[1]
min_markers <- if (length(args) >= 2) as.integer(args[2]) else 5

if (!file.exists(input_file)) {
  stop("File not found: ", input_file)
}

sample_name <- gsub("\\.statepath\\.tsv\\.gz$", "", basename(input_file))
out_dir <- dirname(input_file)

cat("=== INTROGRESSION PLOT ===\n")
cat("Sample:", sample_name, "\n")
cat("Rigidity:", min_markers, "\n\n")

# Chromosome sizes
chr_len <- c(
  chr1 = 308452471, chr2 = 243675191, chr3 = 238017767,
  chr4 = 250330460, chr5 = 226353449, chr6 = 181357234,
  chr7 = 185808916, chr8 = 182411202, chr9 = 163004744,
  chr10 = 152435371
)

# State to genotype mapping (RTIGER style)
state_to_geno <- c(RR = "AA", RH = "AB", HH = "BB")

# Colors (RTIGER style)
geno_colors <- c(AA = "#3366CC", AB = "#663399", BB = "#CC3333")

# ============================================================================
# SMOOTHING FUNCTION
# ============================================================================

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
cat("Total SNPs:", format(nrow(dt), big.mark = ","), "\n")

# Process each chromosome and build continuous segments
all_segments <- list()

for (chr in paste0("chr", 1:10)) {
  dchr <- dt[CHROM == chr]
  if (nrow(dchr) == 0) next

  # Apply smoothing
  orig_states <- dchr$STATE
  probs <- as.matrix(dchr[, .(P_RR, P_RH, P_HH)])
  colnames(probs) <- c("RR", "RH", "HH")

  smooth_st <- smooth_states(orig_states, probs, min_markers)

  # Get continuous segments using rle
  r <- rle(smooth_st)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  for (i in seq_along(r$values)) {
    # Get physical positions
    start_pos <- dchr$POS[starts[i]]
    end_pos <- dchr$POS[ends[i]]

    # Convert state to genotype
    geno <- state_to_geno[r$values[i]]

    all_segments[[length(all_segments) + 1]] <- data.table(
      chrom = chr,
      start = start_pos,
      end = end_pos,
      genotype = geno,
      state = r$values[i],
      n_snps = r$lengths[i]
    )
  }
}

segments <- rbindlist(all_segments)

# ============================================================================
# OUTPUT BED FILE
# ============================================================================

bed_file <- file.path(out_dir, paste0(sample_name, ".genotypes.bed"))

# BED format: chrom start end genotype
bed_out <- segments[, .(chrom, start, end, genotype)]
fwrite(bed_out, bed_file, sep = "\t", col.names = FALSE)

cat("\nBED file saved:", bed_file, "\n")
cat("Preview:\n")
print(head(bed_out, 20))

# Also save detailed version with all info
detailed_file <- file.path(out_dir, paste0(sample_name, ".segments.tsv"))
fwrite(segments, detailed_file, sep = "\t")
cat("\nDetailed segments:", detailed_file, "\n")

# ============================================================================
# CREATE VISUALIZATION (RTIGER-style)
# ============================================================================

cat("\nGenerating plot...\n")

# PNG output
png_file <- file.path(out_dir, paste0(sample_name, ".genotype_plot.png"))

# Calculate plot dimensions
n_chr <- 10
plot_height <- 150 + n_chr * 80  # pixels per chromosome

png(png_file, width = 1400, height = plot_height, res = 100)

# Set up layout
par(mar = c(4, 8, 4, 2), family = "sans")

# Create empty plot
plot(NULL, xlim = c(0, max(chr_len)), ylim = c(0, n_chr + 1),
     xlab = "Position (Mb)", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n",
     main = paste0(sample_name, " - Introgression Map"))

# X-axis in Mb
axis_pos <- seq(0, 300e6, by = 50e6)
axis(1, at = axis_pos, labels = axis_pos / 1e6)

# Plot each chromosome
for (i in 1:n_chr) {
  chr <- paste0("chr", i)
  y_pos <- n_chr - i + 1

  # Chromosome background (gray)
  rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35,
       col = "gray90", border = NA)

  # Get segments for this chromosome
  chr_segs <- segments[chrom == chr]

  # Plot each segment
  for (j in seq_len(nrow(chr_segs))) {
    seg <- chr_segs[j]
    rect(seg$start, y_pos - 0.35, seg$end, y_pos + 0.35,
         col = geno_colors[seg$genotype], border = NA)
  }

  # Chromosome outline
  rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35,
       col = NA, border = "gray40", lwd = 0.5)

  # Chromosome label
  text(-5e6, y_pos, chr, adj = 1, cex = 0.9, font = 2)
}

# Legend
legend("topright",
       legend = c("AA (B73/B73)", "AB (Het)", "BB (Teo/Teo)"),
       fill = geno_colors[c("AA", "AB", "BB")],
       border = "gray40", cex = 0.9, bg = "white")

dev.off()

cat("Plot saved:", png_file, "\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== GENOME SUMMARY ===\n")

# Calculate genome coverage by genotype
genome_total <- sum(chr_len)

geno_summary <- segments[, .(
  total_bp = sum(end - start + 1),
  n_segments = .N
), by = genotype]

geno_summary[, pct := round(100 * total_bp / genome_total, 2)]
geno_summary[, total_Mb := round(total_bp / 1e6, 2)]

print(geno_summary[order(genotype)])

# Chromosome-level summary
cat("\n=== CHROMOSOME SUMMARY ===\n")
chr_summary <- segments[genotype != "AA", .(
  AB_Mb = round(sum((end - start + 1)[genotype == "AB"]) / 1e6, 2),
  BB_Mb = round(sum((end - start + 1)[genotype == "BB"]) / 1e6, 2),
  donor_Mb = round(sum(end - start + 1) / 1e6, 2)
), by = chrom]

chr_summary[, chr_len_Mb := chr_len[chrom] / 1e6]
chr_summary[, donor_pct := round(100 * donor_Mb / chr_len_Mb, 2)]
setorder(chr_summary, chrom)
print(chr_summary)

cat("\n=== DONE ===\n")
cat("Files created:\n")
cat("  BED:", bed_file, "\n")
cat("  TSV:", detailed_file, "\n")
cat("  PNG:", png_file, "\n")
