#!/usr/bin/env Rscript
# =============================================================================
# 10_counts_qc.R — per-pool read-count QC (base R only).
#   Usage: Rscript 10_counts_qc.R <D09 dir> <name_change_final.txt> <output dir>
#   (normally run through 10_counts_qc.sh, which takes the paths from config.sh)
#
#   Outputs:
#     sample_qc.tsv      one row per pool: reads, coverage, ALT and other-base fractions, flags
#     taxon_summary.tsv  medians per donor taxon (and B73 / PurpleCheck controls)
#     site_qc.tsv        per-chromosome site summary
#     qc_plots.png       reads, ALT fraction and other-base fraction by taxon
#
#   What to look for:
#     - B73 controls: ALT_FRAC is the background (errors, paralogs, B73 seed-lot
#       differences). Donor pools should sit clearly above it.
#     - OTHER_FRAC (third-allele reads) estimates the per-pool error rate.
#     - LOW_COVERAGE pools (< 20% of the median read count) need checking.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: Rscript 10_counts_qc.R <D09> <SAMPLE_META> <OUTDIR>")
d09 <- args[1]; meta_file <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- per-pool totals over all chromosomes ------------------------------------
qc_files <- list.files(d09, pattern = "\\.sampleqc\\.tsv$", full.names = TRUE)
if (length(qc_files) == 0) stop("no *.sampleqc.tsv in ", d09)
message("chromosome files: ", length(qc_files))
qc <- do.call(rbind, lapply(qc_files, read.delim, colClasses = c("character", rep("numeric", 4))))
q <- aggregate(cbind(SITES_COVERED, REF_READS, ALT_READS, OTHER_READS) ~ SAMPLE, data = qc, FUN = sum)

site <- do.call(rbind, lapply(list.files(d09, "\\.siteqc\\.tsv$", full.names = TRUE), read.delim))
n_sites <- sum(site$N_SITES)
write.table(site, file.path(outdir, "site_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- metadata (file has Windows line endings; one ID is listed twice) ---------
meta <- read.delim(meta_file, colClasses = "character")
meta[] <- lapply(meta, function(x) trimws(sub("\r$", "", x)))
names(meta)[1:2] <- c("SAMPLE", "GENOTYPE")
dup_ids <- unique(meta$SAMPLE[duplicated(meta$SAMPLE)])
if (length(dup_ids)) message("WARNING: IDs listed more than once in metadata: ", paste(dup_ids, collapse = ", "))
meta <- aggregate(GENOTYPE ~ SAMPLE, data = meta, FUN = function(x) paste(unique(x), collapse = "|"))
is_donor <- grepl("^Z[a-z]\\.", meta$GENOTYPE)
meta$TAXON <- ifelse(is_donor, sub("\\..*$", "", meta$GENOTYPE), meta$GENOTYPE)
meta$ACCESSION <- ifelse(is_donor, sub("_.*$", "", meta$GENOTYPE), meta$GENOTYPE)
meta$TAXON[grepl("\\|", meta$GENOTYPE)] <- "conflicting_metadata"

q <- merge(q, meta, by = "SAMPLE", all.x = TRUE)
q$TAXON[is.na(q$TAXON)] <- "no_metadata"

# ---- per-pool metrics ---------------------------------------------------------
q$READS <- q$REF_READS + q$ALT_READS
q$FRAC_SITES_COVERED <- q$SITES_COVERED / n_sites
q$MEAN_DEPTH_AT_SNPS <- q$READS / n_sites
q$ALT_FRAC <- q$ALT_READS / pmax(q$READS, 1)
q$OTHER_FRAC <- q$OTHER_READS / pmax(q$READS + q$OTHER_READS, 1)
q$LOW_COVERAGE <- q$READS < 0.2 * median(q$READS)
q <- q[order(q$TAXON, q$SAMPLE), ]
write.table(q, file.path(outdir, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- per-taxon summary --------------------------------------------------------
tx <- do.call(rbind, lapply(split(q, q$TAXON), function(d) data.frame(
  TAXON = d$TAXON[1], N_POOLS = nrow(d),
  MEDIAN_READS = median(d$READS),
  MEDIAN_DEPTH_AT_SNPS = round(median(d$MEAN_DEPTH_AT_SNPS), 3),
  MEDIAN_FRAC_SITES_COVERED = round(median(d$FRAC_SITES_COVERED), 3),
  MEDIAN_ALT_FRAC = signif(median(d$ALT_FRAC), 3),
  MEDIAN_OTHER_FRAC = signif(median(d$OTHER_FRAC), 3),
  N_LOW_COVERAGE = sum(d$LOW_COVERAGE))))
write.table(tx, file.path(outdir, "taxon_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- plots --------------------------------------------------------------------
lev <- names(sort(tapply(q$ALT_FRAC, q$TAXON, median)))
q$TAXON_F <- factor(q$TAXON, levels = lev)
png(file.path(outdir, "qc_plots.png"), width = 1800, height = 600, res = 120)
par(mfrow = c(1, 3), mar = c(8, 5, 3, 1))
boxplot(log10(READS + 1) ~ TAXON_F, data = q, las = 2, xlab = "", ylab = "log10(B73+ALT reads at SNPs)",
        main = "Reads per pool")
boxplot(ALT_FRAC ~ TAXON_F, data = q, las = 2, xlab = "", ylab = "ALT reads / (B73+ALT)",
        main = "ALT fraction (B73 = background)")
if ("B73" %in% lev) abline(h = median(q$ALT_FRAC[q$TAXON == "B73"]), lty = 2, col = "red")
boxplot(OTHER_FRAC ~ TAXON_F, data = q, las = 2, xlab = "", ylab = "third-allele reads / all reads",
        main = "Error proxy")
dev.off()

# ---- console summary ----------------------------------------------------------
cat("\nSNP sites with read counts:", n_sites, "\n")
cat("Pools:", nrow(q), " | low coverage:", sum(q$LOW_COVERAGE),
    " | no metadata:", sum(q$TAXON == "no_metadata"), "\n\n")
print(tx, row.names = FALSE)
cat("\nWrote:", file.path(outdir, c("sample_qc.tsv", "taxon_summary.tsv", "site_qc.tsv", "qc_plots.png")), sep = "\n  ")
