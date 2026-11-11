#!/usr/bin/env Rscript
# ============================================================================
# PLOT INTROGRESSION - COMPLETE BLOCKS (RTIGER-style)
# ============================================================================
# This version creates CONTINUOUS BLOCKS like RTIGER by:
#   1. Merging adjacent same-state segments
#   2. Bridging small gaps between same-state blocks
#   3. Filtering out tiny blocks (< min_block_size)
#
# Usage: Rscript plot_introgression_blocks.R <input> [max_gap_kb] [min_block_kb]
#
# Arguments:
#   input: statepath.tsv.gz file
#   max_gap_kb: Maximum gap (kb) to bridge between same-state blocks (default: 500)
#   min_block_kb: Minimum block size (kb) to keep (default: 10)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript plot_introgression_blocks.R <input> [max_gap_kb] [min_block_kb]\n")
  cat("\nArguments:\n")
  cat("  input:        statepath.tsv.gz file\n")
  cat("  max_gap_kb:   Max gap to bridge (default: 500 kb)\n")
  cat("  min_block_kb: Min block size to keep (default: 10 kb)\n")
  stop()
}

input_file <- args[1]
max_gap_bp <- if (length(args) >= 2) as.numeric(args[2]) * 1000 else 500000
min_block_bp <- if (length(args) >= 3) as.numeric(args[3]) * 1000 else 10000
# 10000 kb 5000 kb good
if (!file.exists(input_file)) stop("File not found: ", input_file)

sample_name <- gsub("\\.statepath\\.tsv\\.gz$", "", basename(input_file))
out_dir <- dirname(input_file)

cat("=== COMPLETE BLOCK ANALYSIS ===\n")
cat("Sample:", sample_name, "\n")
cat("Max gap to bridge:", max_gap_bp / 1000, "kb\n")
cat("Min block size:", min_block_bp / 1000, "kb\n\n")

# Chromosome sizes
chr_len <- c(
  chr1 = 308452471, chr2 = 243675191, chr3 = 238017767,
  chr4 = 250330460, chr5 = 226353449, chr6 = 181357234,
  chr7 = 185808916, chr8 = 182411202, chr9 = 163004744,
  chr10 = 152435371
)

state_to_geno <- c(RR = "AA", RH = "AB", HH = "BB")
geno_colors <- c(AA = "#3366CC", AB = "#663399", BB = "#CC3333")

# ============================================================================
# MERGE BLOCKS FUNCTION - This is the key!
# ============================================================================
# Merges adjacent blocks of same state, bridging gaps up to max_gap_bp

merge_blocks <- function(segments, max_gap, min_size) {
  if (nrow(segments) == 0) return(segments)

  # Sort by position
  setorder(segments, start)

  merged <- list()
  current <- segments[1]

  for (i in 2:nrow(segments)) {
    next_seg <- segments[i]
    gap <- next_seg$start - current$end

    # Merge if same state AND gap is small enough
    if (next_seg$genotype == current$genotype && gap <= max_gap) {
      # Extend current block
      current$end <- next_seg$end
      current$n_snps <- current$n_snps + next_seg$n_snps
    } else {
      # Save current and start new
      merged[[length(merged) + 1]] <- current
      current <- next_seg
    }
  }
  # Don't forget last block
  merged[[length(merged) + 1]] <- current

  result <- rbindlist(merged)

  # Update length
  result[, len_bp := end - start]

  # Filter by minimum size
  result <- result[len_bp >= min_size]

  return(result)
}

# ============================================================================
# FILL CHROMOSOME GAPS
# ============================================================================
# Fills gaps at chromosome start/end with flanking state

fill_chromosome <- function(segments, chr_length, chr_name) {
  if (nrow(segments) == 0) {
    # No segments - fill entire chromosome with AA (reference)
    return(data.table(
      chrom = chr_name,
      start = 1,
      end = chr_length,
      genotype = "AA",
      n_snps = 0,
      len_bp = chr_length
    ))
  }

  setorder(segments, start)

  filled <- list()

  # Fill from chromosome start to first segment
  if (segments$start[1] > 1) {
    filled[[1]] <- data.table(
      chrom = chr_name,
      start = 1,
      end = segments$start[1] - 1,
      genotype = segments$genotype[1],  # Use first segment's state
      n_snps = 0,
      len_bp = segments$start[1] - 1
    )
  }

  # Add all existing segments and fill gaps between them
  for (i in 1:nrow(segments)) {
    filled[[length(filled) + 1]] <- segments[i]

    # If there's a gap to next segment, fill it
    if (i < nrow(segments)) {
      gap_start <- segments$end[i] + 1
      gap_end <- segments$start[i + 1] - 1

      if (gap_end > gap_start) {
        # Decide state for gap: use state of longer flanking block
        left_len <- segments$len_bp[i]
        right_len <- segments$len_bp[i + 1]
        gap_state <- if (left_len >= right_len) {
          segments$genotype[i]
        } else {
          segments$genotype[i + 1]
        }

        filled[[length(filled) + 1]] <- data.table(
          chrom = chr_name,
          start = gap_start,
          end = gap_end,
          genotype = gap_state,
          n_snps = 0,
          len_bp = gap_end - gap_start
        )
      }
    }
  }

  # Fill from last segment to chromosome end
  last_end <- segments$end[nrow(segments)]
  if (last_end < chr_length) {
    filled[[length(filled) + 1]] <- data.table(
      chrom = chr_name,
      start = last_end + 1,
      end = chr_length,
      genotype = segments$genotype[nrow(segments)],  # Use last segment's state
      n_snps = 0,
      len_bp = chr_length - last_end
    )
  }

  result <- rbindlist(filled)
  setorder(result, start)

  # Final merge of adjacent same-state blocks
  if (nrow(result) > 1) {
    final_merged <- list()
    current <- result[1]

    for (i in 2:nrow(result)) {
      if (result$genotype[i] == current$genotype) {
        current$end <- result$end[i]
        current$n_snps <- current$n_snps + result$n_snps[i]
        current$len_bp <- current$end - current$start
      } else {
        final_merged[[length(final_merged) + 1]] <- current
        current <- result[i]
      }
    }
    final_merged[[length(final_merged) + 1]] <- current
    result <- rbindlist(final_merged)
  }

  return(result)
}

# ============================================================================
# READ AND PROCESS DATA
# ============================================================================

cat("Reading data...\n")
dt <- fread(cmd = paste("gunzip -c", shQuote(input_file)), showProgress = FALSE)
cat("Total SNPs:", format(nrow(dt), big.mark = ","), "\n")

# Create initial segments from consecutive same-state SNPs
all_segments <- list()

for (chr in paste0("chr", 1:10)) {
  dchr <- dt[CHROM == chr]
  if (nrow(dchr) == 0) next

  setorder(dchr, POS)

  # Get runs of same state
  r <- rle(dchr$STATE)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  for (i in seq_along(r$values)) {
    geno <- state_to_geno[r$values[i]]

    all_segments[[length(all_segments) + 1]] <- data.table(
      chrom = chr,
      start = dchr$POS[starts[i]],
      end = dchr$POS[ends[i]],
      genotype = geno,
      n_snps = r$lengths[i],
      len_bp = dchr$POS[ends[i]] - dchr$POS[starts[i]]
    )
  }
}

segments <- rbindlist(all_segments)

cat("Initial segments:", nrow(segments), "\n")

# ============================================================================
# MERGE AND FILL BLOCKS
# ============================================================================

cat("\nMerging blocks (max gap:", max_gap_bp/1000, "kb)...\n")

complete_blocks <- list()

for (chr in paste0("chr", 1:10)) {
  chr_segs <- segments[chrom == chr]

  if (nrow(chr_segs) == 0) {
    # No data - fill with AA
    complete_blocks[[chr]] <- data.table(
      chrom = chr,
      start = 1,
      end = chr_len[chr],
      genotype = "AA",
      n_snps = 0,
      len_bp = chr_len[chr]
    )
    next
  }

  # Step 1: Merge adjacent same-state blocks
  merged <- merge_blocks(chr_segs, max_gap_bp, min_block_bp)

  # Step 2: Fill chromosome gaps
  filled <- fill_chromosome(merged, chr_len[chr], chr)

  complete_blocks[[chr]] <- filled
}

blocks <- rbindlist(complete_blocks)
blocks[, len_bp := end - start]

cat("Complete blocks:", nrow(blocks), "\n\n")

# ============================================================================
# OUTPUT BED FILE
# ============================================================================

bed_file <- file.path(out_dir, paste0(sample_name, ".complete_blocks.bed"))
bed_out <- blocks[, .(chrom, start, end, genotype)]
fwrite(bed_out, bed_file, sep = "\t", col.names = FALSE)

cat("BED file saved:", bed_file, "\n")
cat("\nComplete blocks:\n")
print(bed_out)

# ============================================================================
# CREATE VISUALIZATION
# ============================================================================

cat("\nGenerating plot...\n")

png_file <- file.path(out_dir, paste0(sample_name, ".complete_blocks.png"))

n_chr <- 10
plot_height <- 150 + n_chr * 80

png(png_file, width = 1400, height = plot_height, res = 100)

par(mar = c(4, 8, 4, 2))

plot(NULL, xlim = c(0, max(chr_len)), ylim = c(0, n_chr + 1),
     xlab = "Position (Mb)", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n",
     main = paste0(sample_name, " - Complete Introgression Blocks"))

axis_pos <- seq(0, 300e6, by = 50e6)
axis(1, at = axis_pos, labels = axis_pos / 1e6)

for (i in 1:n_chr) {
  chr <- paste0("chr", i)
  y_pos <- n_chr - i + 1

  rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35,
       col = "gray90", border = NA)

  chr_blocks <- blocks[chrom == chr]

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

cat("Plot saved:", png_file, "\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== BLOCK SUMMARY ===\n")

block_summary <- blocks[, .(
  n_blocks = .N,
  total_Mb = round(sum(len_bp) / 1e6, 2)
), by = genotype]

print(block_summary)

cat("\n=== CHROMOSOME SUMMARY ===\n")

chr_summary <- blocks[, .(
  n_blocks = .N,
  AA_Mb = round(sum(len_bp[genotype == "AA"]) / 1e6, 2),
  AB_Mb = round(sum(len_bp[genotype == "AB"]) / 1e6, 2),
  BB_Mb = round(sum(len_bp[genotype == "BB"]) / 1e6, 2)
), by = chrom]

chr_summary[, chr_len_Mb := chr_len[chrom] / 1e6]
chr_summary[, donor_pct := round(100 * (AB_Mb + BB_Mb) / chr_len_Mb, 1)]
setorder(chr_summary, chrom)
print(chr_summary)

cat("\n=== DONE ===\n")
cat("Files:\n")
cat("  BED:", bed_file, "\n")
cat("  PNG:", png_file, "\n")




#Usage:                                                                                                                                        
  # Process a folder                                                                                                                            
 # Rscript batch_introgression_analysis.R <input_folder> [max_gap_kb] [min_block_kb] [output_folder]                                             
                                                                                                                                                
  # Example with your settings                                                                                                                  
  #Rscript batch_introgression_analysis.R /path/to/folder 10000 5000                