################################################################################
### LIBRARIES
################################################################################
library(tidyverse)
library(vroom)
library(data.table)
library(purrr)
library(CMplot)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)


################################################################################
### Getting the GWAS results
################################################################################

# Initialize empty list to hold per-chromosome results
gwas_list <- list()

# Loop through chr1 to chr10
for (chr in 1:10) {
  # Construct the file path
  file_path <- paste0("/Users/nirwantandukar/Documents/Github/BZea_genotyping/GAPIT_FarmCPU_DTA_chr", chr, ".rds")
  
  # Read the RDS file
  x <- readRDS(file_path)
  
  # Extract the GWAS results table for TN_maize
  gwas_chr <- x[["GWAS"]]
  
  # Append to the list
  gwas_list[[chr]] <- gwas_chr
}

# Combine all into a single data frame
gwas_raw <- do.call(rbind, gwas_list)
head(gwas_raw)


################################################################################
### LOAD THE REF FILES
################################################################################

ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]



library(GenomicRanges)
library(IRanges)
library(dplyr)
library(tibble)

window_bp <- 25000  # total width; +/- 12.5 kb around SNP

# ---- 1) Build SNP GRanges from gwas_raw ----
snps <- GRanges(
  seqnames = Rle(paste0("chr", as.character(gwas_raw$Chromosome))),  # "chr1" style to match genes_only
  ranges   = IRanges(start = as.integer(gwas_raw$Position),
                     end   = as.integer(gwas_raw$Position))
)

mcols(snps)$SNP    <- as.character(gwas_raw$SNP)
mcols(snps)$pvalue <- suppressWarnings(as.numeric(gwas_raw$P.value))
mcols(snps)$pos0   <- suppressWarnings(as.integer(gwas_raw$Position))  # keep original SNP position

# NA-safe filtering (keep everything with valid p in [0,1])
pv   <- mcols(snps)$pvalue
keep <- !is.na(pv) & pv <= 1 & pv >= 0
snp_f <- snps[keep]

# ---- 2) Expand to window and overlap genes ----
snp_e <- resize(snp_f, width = window_bp + 1, fix = "center")

ol <- findOverlaps(genes_only, snp_e, ignore.strand = TRUE)

# ---- 3) Make annotated table ----
if (length(ol) == 0) {
  annotated_df <- tibble()
} else {
  hg <- queryHits(ol)    # genes
  hs <- subjectHits(ol)  # snps
  
  snp_pos <- mcols(snp_e)$pos0[hs]
  
  annotated_df <- tibble(
    Model    = "MLM",  # or whatever you want to label it
    GeneID   = mcols(genes_only)$ID[hg],
    SNP      = mcols(snp_e)$SNP[hs],
    Chr      = as.character(seqnames(genes_only)[hg]),
    Pos      = snp_pos,
    P.value  = mcols(snp_e)$pvalue[hs],
    Relation = case_when(
      snp_pos >= start(genes_only)[hg] & snp_pos <= end(genes_only)[hg] ~ "within",
      snp_pos <  start(genes_only)[hg]                                   ~ "upstream",
      TRUE                                                               ~ "downstream"
    )
  ) %>%
    distinct(Model, GeneID, SNP, Chr, Pos, P.value, Relation)
}

# check
annotated_df

# Save the annotated results to a CSV file
write.csv(annotated_df, "/Users/nirwantandukar/Documents/Github/BZea_genotyping/results/GAPIT_FarmCPU_FloweringTime_DTA_annotated_results.csv", row.names = FALSE)
