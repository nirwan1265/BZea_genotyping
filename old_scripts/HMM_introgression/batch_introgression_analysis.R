#!/usr/bin/env Rscript
# ============================================================================
# BATCH INTROGRESSION ANALYSIS
# ============================================================================
# Process multiple samples from a folder and generate summary files.
#
# Usage: Rscript batch_introgression_analysis.R <input_folder> [max_gap_kb] [min_block_kb] [output_folder]
#
# Arguments:
#   input_folder:  Folder containing .statepath.tsv.gz files
#   max_gap_kb:    Maximum gap (kb) to bridge between same-state blocks (default: 500)
#   min_block_kb:  Minimum block size (kb) to keep (default: 10)
#   output_folder: Where to save outputs (default: input_folder)
#
# Outputs:
#   Per sample:
#     - <sample>.complete_blocks.bed
#     - <sample>.complete_blocks.png
#     - <sample>.block_summary.tsv
#     - <sample>.chr_summary.tsv
#
#   Aggregate (all samples):
#     - introgression_summary_total.tsv (Ref%, Het%, Teo% per sample)
#     - introgression_summary_per_chr_donor.tsv (donor % per chr per sample)
#     - introgression_summary_per_chr_het.tsv (Het % per chr)
#     - introgression_summary_per_chr_teo.tsv (Teo % per chr)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript batch_introgression_analysis.R <input_folder> [max_gap_kb] [min_block_kb] [output_folder]\n")
  cat("\nArguments:\n")
  cat("  input_folder:  Folder with .statepath.tsv.gz files\n")
  cat("  max_gap_kb:    Max gap to bridge (default: 500 kb)\n")
  cat("  min_block_kb:  Min block size (default: 10 kb)\n")
  cat("  output_folder: Output location (default: input_folder)\n")
  cat("\nExample:\n")
  cat("  Rscript batch_introgression_analysis.R ./GL_files/ 10000 5000\n")
  stop()
}

input_folder <- args[1]
max_gap_bp <- if (length(args) >= 2) as.numeric(args[2]) * 1000 else 500000
min_block_bp <- if (length(args) >= 3) as.numeric(args[3]) * 1000 else 10000
output_folder <- if (length(args) >= 4) args[4] else input_folder

if (!dir.exists(input_folder)) stop("Folder not found: ", input_folder)
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Find all statepath files
statepath_files <- list.files(input_folder, pattern = "\\.statepath\\.tsv\\.gz$", full.names = TRUE)

if (length(statepath_files) == 0) {
  stop("No .statepath.tsv.gz files found in: ", input_folder)
}

cat("=== BATCH INTROGRESSION ANALYSIS ===\n")
cat("Input folder:", input_folder, "\n")
cat("Output folder:", output_folder, "\n")
cat("Max gap:", max_gap_bp / 1000, "kb\n")
cat("Min block:", min_block_bp / 1000, "kb\n")
cat("Files found:", length(statepath_files), "\n\n")

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
# HELPER FUNCTIONS (from plot_introgression_blocks.R)
# ============================================================================

merge_blocks <- function(segments, max_gap, min_size) {
  if (nrow(segments) == 0) return(segments)

  # Remove any rows with NA values in key columns
  segments <- segments[!is.na(start) & !is.na(end) & !is.na(genotype)]
  if (nrow(segments) == 0) return(segments)

  setorder(segments, start)

  merged <- list()
  current <- segments[1]

  for (i in 2:nrow(segments)) {
    next_seg <- segments[i]
    gap <- next_seg$start - current$end

    # Safe comparison: only merge if gap is valid and genotypes match
    if (!is.na(gap) && isTRUE(next_seg$genotype == current$genotype) && gap <= max_gap) {
      current$end <- next_seg$end
      current$n_snps <- current$n_snps + next_seg$n_snps
    } else {
      merged[[length(merged) + 1]] <- current
      current <- next_seg
    }
  }
  merged[[length(merged) + 1]] <- current

  result <- rbindlist(merged)
  result[, len_bp := end - start]
  result <- result[!is.na(len_bp) & len_bp >= min_size]

  return(result)
}

fill_chromosome <- function(segments, chr_length, chr_name) {
  if (nrow(segments) == 0) {
    return(data.table(
      chrom = chr_name, start = 1, end = chr_length,
      genotype = "AA", n_snps = 0, len_bp = chr_length
    ))
  }

  setorder(segments, start)
  filled <- list()

  if (segments$start[1] > 1) {
    filled[[1]] <- data.table(
      chrom = chr_name, start = 1, end = segments$start[1] - 1,
      genotype = segments$genotype[1], n_snps = 0, len_bp = segments$start[1] - 1
    )
  }

  for (i in 1:nrow(segments)) {
    filled[[length(filled) + 1]] <- segments[i]

    if (i < nrow(segments)) {
      gap_start <- segments$end[i] + 1
      gap_end <- segments$start[i + 1] - 1

      if (gap_end > gap_start) {
        left_len <- segments$len_bp[i]
        right_len <- segments$len_bp[i + 1]
        gap_state <- if (left_len >= right_len) segments$genotype[i] else segments$genotype[i + 1]

        filled[[length(filled) + 1]] <- data.table(
          chrom = chr_name, start = gap_start, end = gap_end,
          genotype = gap_state, n_snps = 0, len_bp = gap_end - gap_start
        )
      }
    }
  }

  last_end <- segments$end[nrow(segments)]
  if (last_end < chr_length) {
    filled[[length(filled) + 1]] <- data.table(
      chrom = chr_name, start = last_end + 1, end = chr_length,
      genotype = segments$genotype[nrow(segments)], n_snps = 0,
      len_bp = chr_length - last_end
    )
  }

  result <- rbindlist(filled)
  setorder(result, start)

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
# PROCESS SINGLE SAMPLE
# ============================================================================

process_sample <- function(input_file, out_dir, max_gap_bp, min_block_bp) {
  sample_name <- gsub("\\.statepath\\.tsv\\.gz$", "", basename(input_file))

  # Read data
  dt <- fread(cmd = paste("gunzip -c", shQuote(input_file)), showProgress = FALSE)

  # Check if file was empty or malformed
  if (nrow(dt) == 0 || !"STATE" %in% names(dt)) {
    warning(paste("Empty or malformed file:", input_file))
    return(NULL)
  }

  # Filter to valid states only (RR, RH, HH)
  valid_states <- c("RR", "RH", "HH")
  dt <- dt[STATE %in% valid_states]
  if (nrow(dt) == 0) {
    warning(paste("No valid states in file:", input_file))
    return(NULL)
  }

  # Create segments
  all_segments <- list()

  for (chr in paste0("chr", 1:10)) {
    dchr <- dt[CHROM == chr]
    if (nrow(dchr) == 0) next

    setorder(dchr, POS)
    r <- rle(dchr$STATE)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L

    for (i in seq_along(r$values)) {
      state_val <- r$values[i]
      # Skip if state is NA or not in mapping
      if (is.na(state_val) || !state_val %in% names(state_to_geno)) next

      geno <- state_to_geno[state_val]

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

  if (length(all_segments) == 0) {
    warning(paste("No segments created for:", input_file))
    return(NULL)
  }

  segments <- rbindlist(all_segments)

  # Merge and fill blocks
  complete_blocks <- list()

  for (chr in paste0("chr", 1:10)) {
    chr_segs <- segments[chrom == chr]

    if (nrow(chr_segs) == 0) {
      complete_blocks[[chr]] <- data.table(
        chrom = chr, start = 1, end = chr_len[chr],
        genotype = "AA", n_snps = 0, len_bp = chr_len[chr]
      )
      next
    }

    merged <- merge_blocks(chr_segs, max_gap_bp, min_block_bp)
    filled <- fill_chromosome(merged, chr_len[chr], chr)
    complete_blocks[[chr]] <- filled
  }

  blocks <- rbindlist(complete_blocks)
  blocks[, len_bp := end - start]

  # === OUTPUT BED FILE ===
  bed_file <- file.path(out_dir, paste0(sample_name, ".complete_blocks.bed"))
  bed_out <- blocks[, .(chrom, start, end, genotype)]
  fwrite(bed_out, bed_file, sep = "\t", col.names = FALSE)

  # === BLOCK SUMMARY ===
  block_summary <- blocks[, .(
    n_blocks = .N,
    total_Mb = round(sum(len_bp) / 1e6, 2)
  ), by = genotype]

  summary_file <- file.path(out_dir, paste0(sample_name, ".block_summary.tsv"))
  fwrite(block_summary, summary_file, sep = "\t")

  # === CHROMOSOME SUMMARY ===
  chr_summary <- blocks[, .(
    n_blocks = .N,
    AA_Mb = round(sum(len_bp[genotype == "AA"]) / 1e6, 2),
    AB_Mb = round(sum(len_bp[genotype == "AB"]) / 1e6, 2),
    BB_Mb = round(sum(len_bp[genotype == "BB"]) / 1e6, 2)
  ), by = chrom]

  chr_summary[, chr_len_Mb := chr_len[chrom] / 1e6]
  chr_summary[, donor_pct := round(100 * (AB_Mb + BB_Mb) / chr_len_Mb, 1)]
  chr_summary[, het_pct := round(100 * AB_Mb / chr_len_Mb, 1)]
  chr_summary[, teo_pct := round(100 * BB_Mb / chr_len_Mb, 1)]
  setorder(chr_summary, chrom)

  chr_file <- file.path(out_dir, paste0(sample_name, ".chr_summary.tsv"))
  fwrite(chr_summary, chr_file, sep = "\t")

  # === PNG VISUALIZATION ===
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

    rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35, col = "gray90", border = NA)

    chr_blocks <- blocks[chrom == chr]
    for (j in seq_len(nrow(chr_blocks))) {
      blk <- chr_blocks[j]
      rect(blk$start, y_pos - 0.35, blk$end, y_pos + 0.35,
           col = geno_colors[blk$genotype], border = NA)
    }

    rect(0, y_pos - 0.35, chr_len[chr], y_pos + 0.35, col = NA, border = "gray40", lwd = 0.5)
    text(-5e6, y_pos, chr, adj = 1, cex = 0.9, font = 2)
  }

  legend("topright",
         legend = c("AA (B73/B73)", "AB (Het)", "BB (Teo/Teo)"),
         fill = geno_colors[c("AA", "AB", "BB")],
         border = "gray40", cex = 0.9, bg = "white")

  dev.off()

  # === RETURN SUMMARY FOR AGGREGATION ===
  total_genome <- sum(chr_len)

  total_AA <- blocks[genotype == "AA", sum(len_bp)]
  total_AB <- blocks[genotype == "AB", sum(len_bp)]
  total_BB <- blocks[genotype == "BB", sum(len_bp)]

  sample_summary <- data.table(
    sample = sample_name,
    Ref_pct = round(100 * total_AA / total_genome, 2),
    Het_pct = round(100 * total_AB / total_genome, 2),
    Teo_pct = round(100 * total_BB / total_genome, 2),
    Donor_pct = round(100 * (total_AB + total_BB) / total_genome, 2)
  )

  # Per-chromosome stats
  chr_stats <- chr_summary[, .(
    sample = sample_name,
    chrom,
    donor_pct,
    het_pct,
    teo_pct
  )]

  return(list(
    sample_summary = sample_summary,
    chr_stats = chr_stats,
    n_blocks = nrow(blocks)
  ))
}

# ============================================================================
# PROCESS ALL SAMPLES
# ============================================================================

all_sample_summaries <- list()
all_chr_stats <- list()

for (i in seq_along(statepath_files)) {
  f <- statepath_files[i]
  sample_name <- gsub("\\.statepath\\.tsv\\.gz$", "", basename(f))

  cat(sprintf("[%d/%d] Processing %s...\n", i, length(statepath_files), sample_name))

  result <- tryCatch({
    process_sample(f, output_folder, max_gap_bp, min_block_bp)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (!is.null(result)) {
    all_sample_summaries[[i]] <- result$sample_summary
    all_chr_stats[[i]] <- result$chr_stats
    cat(sprintf("  -> %d blocks\n", result$n_blocks))
  }
}

# ============================================================================
# AGGREGATE SUMMARIES
# ============================================================================

cat("\n=== CREATING AGGREGATE SUMMARIES ===\n")

# Filter out NULL entries
all_sample_summaries <- all_sample_summaries[!sapply(all_sample_summaries, is.null)]
all_chr_stats <- all_chr_stats[!sapply(all_chr_stats, is.null)]

if (length(all_sample_summaries) == 0) {
  cat("WARNING: No samples processed successfully!\n")
  # Create empty summary files so downstream scripts don't fail
  total_file <- file.path(output_folder, "introgression_summary_total.tsv")
  fwrite(data.table(sample = character(), Ref_pct = numeric(), Het_pct = numeric(),
                    Teo_pct = numeric(), Donor_pct = numeric()),
         total_file, sep = "\t")
  cat("Empty summary:", total_file, "\n")
  stop("No samples processed successfully")
}

# 1. Total introgression summary
total_summary <- rbindlist(all_sample_summaries, fill = TRUE)
total_file <- file.path(output_folder, "introgression_summary_total.tsv")
fwrite(total_summary, total_file, sep = "\t")
cat("Total summary:", total_file, "\n")

# 2. Per-chromosome summaries (wide format)
if (length(all_chr_stats) > 0) {
  all_chr <- rbindlist(all_chr_stats, fill = TRUE)

  # Donor % per chromosome
  tryCatch({
    donor_wide <- dcast(all_chr, sample ~ chrom, value.var = "donor_pct")
    # Only reorder columns that exist
    cols_present <- intersect(c("sample", paste0("chr", 1:10)), names(donor_wide))
    setcolorder(donor_wide, cols_present)
    donor_file <- file.path(output_folder, "introgression_summary_per_chr_donor.tsv")
    fwrite(donor_wide, donor_file, sep = "\t")
    cat("Per-chr donor %:", donor_file, "\n")
  }, error = function(e) cat("  Warning: Could not create donor summary:", conditionMessage(e), "\n"))

  # Het % per chromosome
  tryCatch({
    het_wide <- dcast(all_chr, sample ~ chrom, value.var = "het_pct")
    cols_present <- intersect(c("sample", paste0("chr", 1:10)), names(het_wide))
    setcolorder(het_wide, cols_present)
    het_file <- file.path(output_folder, "introgression_summary_per_chr_het.tsv")
    fwrite(het_wide, het_file, sep = "\t")
    cat("Per-chr het %:", het_file, "\n")
  }, error = function(e) cat("  Warning: Could not create het summary:", conditionMessage(e), "\n"))

  # Teo % per chromosome
  tryCatch({
    teo_wide <- dcast(all_chr, sample ~ chrom, value.var = "teo_pct")
    cols_present <- intersect(c("sample", paste0("chr", 1:10)), names(teo_wide))
    setcolorder(teo_wide, cols_present)
    teo_file <- file.path(output_folder, "introgression_summary_per_chr_teo.tsv")
    fwrite(teo_wide, teo_file, sep = "\t")
    cat("Per-chr teo %:", teo_file, "\n")
  }, error = function(e) cat("  Warning: Could not create teo summary:", conditionMessage(e), "\n"))
}

# ============================================================================
# PRINT SUMMARIES
# ============================================================================

cat("\n=== TOTAL INTROGRESSION SUMMARY ===\n")
print(total_summary)

cat("\n=== PER-CHROMOSOME DONOR % ===\n")
print(donor_wide)

cat("\n=== PER-CHROMOSOME HET % ===\n")
print(het_wide)

cat("\n=== PER-CHROMOSOME TEO % ===\n")
print(teo_wide)

cat("\n=== DONE ===\n")
cat("Processed", length(statepath_files), "samples\n")
cat("\nOutput files in:", output_folder, "\n")
