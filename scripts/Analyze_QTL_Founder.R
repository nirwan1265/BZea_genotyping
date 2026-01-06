# Load required libraries
library(GridLMM)          # for mixed-model GWAS/QTL mapping
library(data.table)       # for fast data frame operations (optional)
library(dplyr)            # for data manipulation (optional)

# 1. Read in phenotype data (BLUE values for flowering time)
pheno <- read.csv("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv")  # assuming a CSV with columns: Line, Family, Flowering_TIME_BLUE (for example)


# Add a column called Family with the first two strings before the "." in the new_genotype column
pheno$Family <- sub("\\..*$", "", pheno$new_genotype)
colnames(pheno)[1] <- c("Line")
pheno <- pheno %>% select(-c(Sequencing_ID,DTA))
colnames(pheno)[2] <- c("Flowering_TIME_BLUE")

# Ensure the phenotype data has an identifier for each line that matches genotype data
# E.g., columns: Line (unique line ID), Family, Trait (BLUE)
pheno <- pheno %>% mutate(Family = factor(Family))  # ensure Family is a factor
head(pheno)
table(pheno$Family)


# 2. Read or prepare genotype state probability data (RTIGER output)
# Assume RTIGER 2-state posterior probabilities are in a matrix or file where rows = lines, cols = marker positions
# For example, a CSV where first column is Line ID, then columns for each marker with probability of Teosinte allele
geno_probs <- fread("RTiger_binned_dosage_genome_step1e+05.csv")
str(geno_probs)

# Convert to matrix and set rownames as line IDs:
geno_mat <- as.matrix(geno_probs[ , -1])                     # all columns except the first (which might be line ID)
rownames(geno_mat) <- geno_probs[[1]]                        # first column as rownames (line IDs)
colnames(geno_mat)
# The values in geno_mat are probabilities of having the donor (teosinte) allele at each marker.
# (If it's a 3-state output, you would combine states: e.g., P(donor allele) = P(het) + P(homo_teo)*1,
#  effectively the expected dosage of teosinte allele at that locus.)

# 3. (Optional) Construct window-based introgression dosage matrix.
# Define window size (e.g., 10 Mb physical or a fixed number of markers)
window_size <- 100000  # 10 Mb for example

# Suppose we have marker positions in a vector `marker_pos` of length = ncol(geno_mat),
# and corresponding chromosomes in `marker_chr`.
parse_geno_colnames <- function(mk) {
  mk <- as.character(mk)
  chr <- sub(":.*$", "", mk)
  pos <- as.integer(sub("^.*?:", "", mk))
  Chrom <- as.integer(sub("^chr", "", chr))
  data.frame(Marker = mk, chr = chr, Chrom = Chrom, pos = pos, stringsAsFactors = FALSE)
}

marker_map <- parse_geno_colnames(colnames(geno_mat))

head(marker_map)
table(marker_map$chr)

marker_chr <- marker_map$chr
marker_pos <- marker_map$pos


# Create windows: for simplicity, let's do it per chromosome
window_markers <- list()
window_names <- c()
for(chr in unique(marker_chr)) {
  chr_mask <- marker_chr == chr
  chr_positions <- marker_pos[chr_mask]
  chr_markers <- which(chr_mask)
  # break the chromosome into windows of specified size
  breaks <- seq(min(chr_positions), max(chr_positions), by = window_size)
  for(i in seq_along(breaks)) {
    start <- breaks[i]
    end <- start + window_size - 1
    # get markers in this window
    markers_in_window <- chr_markers[chr_positions >= start & chr_positions < end]
    if(length(markers_in_window) == 0) next
    window_name <- paste0("Chr", chr, "_", start, "-", end)
    window_names <- c(window_names, window_name)
    window_markers[[window_name]] <- markers_in_window
  }
}

# Now compute window dosage matrix: average donor probability in each window for each line
window_mat <- sapply(window_markers, function(idx){
  if(length(idx)==1) {
    # if only one marker in window, just use its probability
    geno_mat[, idx]
  } else {
    # average probability across markers in the window
    rowMeans(geno_mat[, idx, drop=FALSE])
  }
})
rownames(window_mat) <- rownames(geno_mat)
colnames(window_mat) <- window_names
# window_mat now has entries from 0 to 1 indicating the fraction of each window that is teosinte in each line.
# Many will be ~0 or ~1 for fixed segments; intermediate values could occur if a window partially overlaps a crossover or heterozygous region.

# 4. Read kinship matrices (LOCO) prepared from PLINK or other tools.
# Assume we have one kinship matrix per chromosome in, e.g., CSV or binary format.
# Here, let's assume they are .csv with row/col names as line IDs:
read_K_matrix_txt <- function(f) {
  Kdf <- data.table::fread(f, data.table = FALSE)
  
  # first column should be sample IDs (header name "IID")
  id_col <- names(Kdf)[1]
  ids <- as.character(Kdf[[id_col]])
  
  # remove id column, convert rest to numeric matrix
  K <- as.matrix(Kdf[, -1, drop = FALSE])
  storage.mode(K) <- "double"
  
  rownames(K) <- ids
  
  # column names: if first column name is "IID", remaining names are the column IDs
  colnames(K) <- names(Kdf)[-1]
  
  # sanity: square?
  if (nrow(K) != ncol(K)) {
    warning("K is not square: ", nrow(K), " x ", ncol(K),
            " file=", basename(f),
            "\nThis usually means the header/format isn't a proper square matrix.")
  }
  
  # enforce symmetry if minor numeric noise
  # K <- (K + t(K))/2
  
  K
}

#K <- read_K_matrix_txt("/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K/K_matrix_chr1.txt")
#dim(K)
#head(rownames(K))
#head(colnames(K))  
# Now kinship_list is a named list of kinship matrices for chr1...chrN (with N=10 for maize).
# Each K excludes markers of that chromosome.

# 5. Prepare the analysis data frame, ensuring it contains the columns for line ID and family.
data <- pheno  # using the phenotype data as base
# Make sure the data frame has a column that identifies lines matching the kinship matrix rownames and geno_mat rownames.

# Assuming pheno$Line is that identifier:
identical(as.character(pheno$Line), rownames(geno_mat))  # ideally TRUE or ensure to match order if not.
length(pheno$Line)
length(rownames(geno_mat))


# IDs present in genotype
ids_g <- rownames(geno_mat)

# subset phenotype to genotyped lines
pheno_sub <- pheno[pheno$Line %in% ids_g, , drop = FALSE]

# reorder phenotype to match geno_mat row order
pheno_sub <- pheno_sub[match(ids_g, pheno_sub$Line), , drop = FALSE]

stopifnot(identical(pheno_sub$Line, ids_g))
nrow(pheno_sub)

pheno <- pheno_sub





# geno_mat: rows = Line IDs (100)
# marker_map: Marker/chr/pos from colnames(geno_mat)
# pheno: your full 906-row phenotype

stripB <- function(z) sub("\\.B$", "", as.character(z))

# make an analysis phenotype that matches geno_mat rows
pheno$Line <- stripB(pheno$Line)
ids <- stripB(rownames(geno_mat))
rownames(geno_mat) <- ids

pheno_sub <- pheno[pheno$Line %in% ids, , drop=FALSE]
pheno_sub <- pheno_sub[match(ids, pheno_sub$Line), , drop=FALSE]
stopifnot(identical(pheno_sub$Line, ids))

# if not already in pheno
pheno_sub$Family <- sub("\\..*$", "", pheno_sub$Line)
pheno_sub$Family <- factor(pheno_sub$Family)



library(data.table)

read_K_matrix_txt <- function(f) {
  Kdf <- fread(f, data.table = FALSE)
  ids <- as.character(Kdf[[1]])
  K <- as.matrix(Kdf[, -1, drop=FALSE])
  storage.mode(K) <- "double"
  rownames(K) <- ids
  colnames(K) <- names(Kdf)[-1]
  
  # if colnames are not IDs but still square, force colnames=rownames
  if (nrow(K) == ncol(K) && !all(colnames(K) %in% rownames(K))) {
    colnames(K) <- rownames(K)
  }
  K
}

kinship_list <- list()
for (chr in unique(marker_map$chr)) {
  f <- paste0("/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K/K_matrix_", chr, ".txt")
  K <- read_K_matrix_txt(f)
  
  # normalize IDs
  rownames(K) <- stripB(rownames(K))
  colnames(K) <- stripB(colnames(K))
  
  # subset & order to ids
  keep <- ids[ids %in% rownames(K)]
  if (length(keep) < 20) next
  
  K <- K[keep, keep, drop=FALSE]
  
  # IMPORTANT: if some ids missing in K, drop those from geno/pheno too for that chr
  kinship_list[[chr]] <- K
}

impute_colmean <- function(X) {
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (anyNA(v)) {
      mu <- mean(v, na.rm=TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      X[, j] <- v
    }
  }
  X
}

# 6. Combined analysis: Scan each chromosome using GridLMM
results_combined <- list()

for (chr in unique(marker_map$chr)) {
  
  markers_on_chr <- which(marker_map$chr == chr)
  X_chr <- geno_mat[, markers_on_chr, drop = FALSE]
  
  K_chr <- kinship_list[[chr]]
  if (is.null(K_chr)) next
  
  keep <- intersect(ids, rownames(K_chr))
  if (length(keep) < 20) next
  
  # phenotype aligned to keep
  ydat <- data.frame(
    Line   = keep,
    y      = pheno_sub$Flowering_TIME_BLUE[match(keep, pheno_sub$Line)],
    Family = pheno_sub$Family[match(keep, pheno_sub$Line)],
    stringsAsFactors = FALSE
  )
  
  # subset & order X and K to keep
  X_use <- X_chr[keep, , drop = FALSE]
  K_use <- K_chr[keep, keep, drop = FALSE]
  
  # impute NAs in X
  X_use <- impute_colmean(X_use)
  
  # GWAS
  res <- GridLMM_GWAS(
    formula         = y ~ Family + (1|Line),  # <-- required random effect term
    test_formula    = ~ 1,
    reduced_formula = ~ 0,
    data            = ydat,
    X               = X_use,
    X_ID            = "Line",
    relmat          = list(Line = K_use),
    method          = "ML",
    centerX         = FALSE,
    scaleX          = FALSE,
    fillNAX         = FALSE,
    verbose         = FALSE,
    mc.cores        = 1
  )
  
  out <- as.data.frame(res$results)
  
  # attach marker IDs in the same order as X columns
  # GridLMM_GWAS returns one row per marker in X, same order.
  out$Marker <- colnames(X_use)
  out$chr <- chr
  out$pos <- as.integer(sub("^.*?:", "", out$Marker))
  
  results_combined[[chr]] <- out
}

# combine
res_all <- do.call(rbind, results_combined)
head(res_all)

# Standardize the results
standardize_gridlmm <- function(out) {
  out <- as.data.frame(out)
  if (!"p_value_ML" %in% names(out)) stop("No p_value_ML column found.")
  out$p <- as.numeric(out$p_value_ML)
  out$p <- pmin(pmax(out$p, 1e-300), 1)
  out$mlog10p <- -log10(out$p)
  out
}
res_all <- standardize_gridlmm(res_all)
res_all <- res_all[order(res_all$p), ]
head(res_all[, c("Marker","chr","pos","p","mlog10p")], 10)


# 7. Family-wise analysis: loop over each family and perform the scan within that subset.
run_gridlmm_family <- function(fam, pheno_sub, geno_mat, marker_map, kinship_list) {
  
  ids_fam <- pheno_sub$Line[pheno_sub$Family == fam]
  ids_fam <- ids_fam[ids_fam %in% rownames(geno_mat)]
  
  if (length(ids_fam) < 20) {
    message("Skipping ", fam, " (n<20): ", length(ids_fam))
    return(NULL)
  }
  
  out_list <- list()
  
  for (chr in unique(marker_map$chr)) {
    
    K_chr <- kinship_list[[chr]]
    if (is.null(K_chr)) next
    
    keep <- intersect(ids_fam, rownames(K_chr))
    if (length(keep) < 20) next
    
    markers_on_chr <- which(marker_map$chr == chr)
    X_chr <- geno_mat[, markers_on_chr, drop = FALSE]
    X_use <- X_chr[keep, , drop = FALSE]
    K_use <- K_chr[keep, keep, drop = FALSE]
    
    # phenotype for this family, aligned to keep
    ydat <- data.frame(
      Line = keep,
      y    = pheno_sub$Flowering_TIME_BLUE[match(keep, pheno_sub$Line)],
      stringsAsFactors = FALSE
    )
    
    # impute X NAs (same function you already have)
    X_use <- impute_colmean(X_use)
    
    stopifnot(identical(ydat$Line, rownames(K_use)))
    stopifnot(identical(ydat$Line, rownames(X_use)))
    
    res <- GridLMM_GWAS(
      formula         = y ~ 1 + (1|Line),
      test_formula    = ~ 1,
      reduced_formula = ~ 0,
      data            = ydat,
      X               = X_use,
      X_ID            = "Line",
      relmat          = list(Line = K_use),
      method          = "ML",
      centerX         = FALSE,
      scaleX          = FALSE,
      fillNAX         = FALSE,
      verbose         = FALSE,
      mc.cores        = 1
    )
    
    out <- as.data.frame(res$results)
    out$Marker <- colnames(X_use)
    out$chr <- chr
    out$pos <- as.integer(sub("^.*?:", "", out$Marker))
    out$Family <- fam
    
    out_list[[chr]] <- out
  }
  
  fam_res <- do.call(rbind, out_list)
  fam_res <- standardize_gridlmm(fam_res)
  fam_res <- fam_res[order(fam_res$p), ]
  fam_res
}

families <- levels(pheno_sub$Family)

res_by_family <- lapply(families, function(fam) {
  run_gridlmm_family(fam, pheno_sub, geno_mat, marker_map, kinship_list)
})
names(res_by_family) <- families
res_by_family <- res_by_family[!vapply(res_by_family, is.null, logical(1))]

# quick peek: top hit per family
lapply(res_by_family, function(df) df[1, c("Family","Marker","chr","pos","p","mlog10p")])

# 8. Peak calling
call_peaks_chr <- function(df_chr, drop = 1) {
  df_chr <- df_chr[order(df_chr$pos), ]
  ipeak <- which.max(df_chr$mlog10p)
  peak <- df_chr[ipeak, ]
  
  thr <- peak$mlog10p - drop
  keep <- df_chr$mlog10p >= thr
  
  # contiguous region around peak
  left <- max(which(keep[1:ipeak]), na.rm = TRUE)
  right <- ipeak - 1 + min(which(rev(keep[ipeak:length(keep)]) == FALSE), na.rm = TRUE) - 1
  if (!is.finite(left)) left <- ipeak
  if (!is.finite(right)) right <- ipeak
  
  ci_start <- df_chr$pos[left]
  ci_end   <- df_chr$pos[right]
  
  data.frame(
    chr = peak$chr,
    peak_pos = peak$pos,
    peak_marker = peak$Marker,
    peak_p = peak$p,
    peak_mlog10p = peak$mlog10p,
    ci_start = ci_start,
    ci_end = ci_end
  )
}
# peaks for combined
peaks_combined <- do.call(rbind, lapply(split(res_all, res_all$chr), call_peaks_chr))
peaks_combined

peaks_by_family <- lapply(res_by_family, function(df) {
  do.call(rbind, lapply(split(df, df$chr), call_peaks_chr))
})
peaks_by_family




# 9. Annotation:

## ============================================================
## Annotate GridLMM / RTIGER-GridLMM results (res_all) with genes
##   - Uses GFF3/GTF gene coordinates + +/- window around Marker pos
##   - Returns:
##       $long    : one row per (marker x gene) overlap
##       $summary : one row per marker, collapsed gene list
## ============================================================

## Fix: dplyr `.data$P.value` error happens because dplyr WILL ERROR
## if you reference a column that doesn't exist (even inside coalesce()).
## Your res_all already has `p`, so we should just USE `p` (or p_value_ML fallback)
## and NEVER touch `P.value`.

## ============================================================
## Robust annotation for GridLMM introgression-bin GWAS results
## (This avoids the seqlevels/seqnames mismatch that broke your earlier function.)
##
## Key fix:
##   NEVER do seqnames(gr) <- "chr1" if seqlevels(gr) are "1","2",...
##   Instead rename SEQLEVELS first (GenomeInfoDb::renameSeqlevels()).
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

## ---------------------------
## helpers
## ---------------------------
.fix_chr <- function(x, chr_style = c("chr", "nochr")) {
  chr_style <- match.arg(chr_style)
  x <- as.character(x)
  x <- ifelse(is.na(x), NA_character_, x)
  x <- gsub("^Chr", "chr", x)
  x <- gsub("^CHR", "chr", x)
  
  # extract number from things like "chr1", "1", "chr01"
  num <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  out <- ifelse(!is.na(num),
                if (chr_style == "chr") paste0("chr", num) else as.character(num),
                NA_character_)
  out
}

.pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

.parse_chr_pos_from_marker <- function(x) {
  x <- as.character(x)
  chr <- str_extract(x, "^[^:]+")
  pos <- suppressWarnings(as.integer(str_extract(x, "(?<=:)\\d+$")))
  list(chr = chr, pos = pos)
}

## ---------------------------
## 1) Harmonize/standardize your GWAS bins table
##    Works with:
##      - your res_all (Marker/chr/pos/p/mlog10p)
##      - your familywise gwas_bins (bin_start/bin_end/family, etc.)
## ---------------------------
prep_gwas_bins <- function(res_all,
                           chr_style = c("chr","nochr"),
                           sig_thr_mlog10 = NULL) {
  chr_style <- match.arg(chr_style)
  df <- as.data.frame(res_all)
  
  # marker column
  marker_col <- .pick_col(df, c("marker","Marker","X_ID"))
  if (is.null(marker_col)) stop("No marker column found. Need one of: marker / Marker / X_ID")
  
  # chr + pos columns (or parse from marker)
  chr_col <- .pick_col(df, c("chr","CHR","chrom","chromosome","seqnames","chr.x","chr.y"))
  pos_col <- .pick_col(df, c("pos","POS","position","bp","Pos","pos.x","pos.y"))
  
  if (is.null(chr_col) || is.null(pos_col)) {
    parsed <- .parse_chr_pos_from_marker(df[[marker_col]])
    if (is.null(chr_col)) df$chr <- parsed$chr else df$chr <- df[[chr_col]]
    if (is.null(pos_col)) df$pos <- parsed$pos else df$pos <- df[[pos_col]]
  } else {
    df$chr <- df[[chr_col]]
    df$pos <- df[[pos_col]]
  }
  
  # bin start/end (if absent, treat as point bins at pos)
  bs_col <- .pick_col(df, c("bin_start","start","Start","bin_start.x","bin_start.y"))
  be_col <- .pick_col(df, c("bin_end","end","End","bin_end.x","bin_end.y"))
  if (is.null(bs_col) || is.null(be_col)) {
    df$bin_start <- as.integer(df$pos)
    df$bin_end   <- as.integer(df$pos)
  } else {
    df$bin_start <- as.integer(df[[bs_col]])
    df$bin_end   <- as.integer(df[[be_col]])
  }
  
  # family optional
  fam_col <- .pick_col(df, c("family","Family","FAMILY"))
  if (!is.null(fam_col)) df$family <- as.character(df[[fam_col]]) else df$family <- NA_character_
  
  # p / mlog10p
  p_col <- .pick_col(df, c("p","p_value_ML","P","pvalue","pval","p.value"))
  if (is.null(p_col)) stop("No p-value column found. Need one of: p / p_value_ML / P / pvalue / pval")
  
  mlog_col <- .pick_col(df, c("mlog10p","-log10p","neglog10p"))
  if (is.null(mlog_col)) {
    ptmp <- suppressWarnings(as.numeric(df[[p_col]]))
    ptmp <- pmin(pmax(ptmp, 1e-300), 1)
    df$mlog10p <- -log10(ptmp)
  } else {
    df$mlog10p <- suppressWarnings(as.numeric(df[[mlog_col]]))
  }
  df$p <- suppressWarnings(as.numeric(df[[p_col]]))
  df$p <- pmin(pmax(df$p, 1e-300), 1)
  
  out <- df %>%
    transmute(
      marker = as.character(.data[[marker_col]]),
      family = .data$family,
      chr    = .fix_chr(.data$chr, chr_style = chr_style),
      pos    = as.integer(.data$pos),
      bin_start = as.integer(.data$bin_start),
      bin_end   = as.integer(.data$bin_end),
      mlog10p = as.numeric(.data$mlog10p),
      p = as.numeric(.data$p)
    ) %>%
    filter(!is.na(chr), !is.na(bin_start), !is.na(bin_end)) %>%
    arrange(p)
  
  if (!is.null(sig_thr_mlog10)) out <- out %>% filter(mlog10p >= sig_thr_mlog10)
  out
}

## ---------------------------
## 2) Load genes from GFF3 and FORCE chr style safely
##    (This is the part that avoids the seqlevels/seqnames error.)
## ---------------------------
load_genes_from_gff <- function(gff_file,
                                chr_style = c("chr","nochr"),
                                keep_chr = 1:10) {
  chr_style <- match.arg(chr_style)
  stopifnot(file.exists(gff_file))
  
  ref_gr <- rtracklayer::import(gff_file)
  genes  <- ref_gr[mcols(ref_gr)$type == "gene"]
  
  # robust gene id/name
  get_gene_id <- function(gr) {
    m <- mcols(gr)
    if ("ID" %in% names(m)) return(as.character(m$ID))
    if ("gene_id" %in% names(m)) return(as.character(m$gene_id))
    if ("Name" %in% names(m)) return(as.character(m$Name))
    rep(NA_character_, length(gr))
  }
  get_gene_name <- function(gr) {
    m <- mcols(gr)
    if ("Name" %in% names(m)) return(as.character(m$Name))
    if ("gene" %in% names(m)) return(as.character(m$gene))
    rep(NA_character_, length(gr))
  }
  
  mcols(genes)$gene_id   <- sub("^gene:", "", get_gene_id(genes))
  mcols(genes)$gene_name <- get_gene_name(genes)
  
  # ---- SAFE renaming of seqlevels (critical) ----
  old <- seqlevels(genes)
  old_num <- suppressWarnings(as.integer(str_extract(old, "\\d+")))
  keep_set <- paste0("chr", keep_chr)
  
  # build mapping old_level -> new_level
  new_levels <- ifelse(!is.na(old_num),
                       if (chr_style == "chr") paste0("chr", old_num) else as.character(old_num),
                       old)
  
  map <- setNames(new_levels, old)  # names=old, values=new
  genes <- GenomeInfoDb::renameSeqlevels(genes, map)
  
  if (chr_style == "chr") {
    genes <- keepSeqlevels(genes, keep_set, pruning.mode = "coarse")
  } else {
    genes <- keepSeqlevels(genes, as.character(keep_chr), pruning.mode = "coarse")
  }
  
  genes
}

## ---------------------------
## 3) Annotate bins -> genes within +/- window_bp
## ---------------------------
annotate_bins_to_genes <- function(gwas_bins, genes_gr, window_bp = 25000L) {
  gwas_bins <- as.data.frame(gwas_bins)
  
  bins_gr <- GRanges(
    seqnames = gwas_bins$chr,
    ranges   = IRanges(start = gwas_bins$bin_start, end = gwas_bins$bin_end)
  )
  mcols(bins_gr)$marker  <- gwas_bins$marker
  mcols(bins_gr)$family  <- gwas_bins$family
  mcols(bins_gr)$mlog10p <- gwas_bins$mlog10p
  mcols(bins_gr)$p       <- gwas_bins$p
  mcols(bins_gr)$pos     <- gwas_bins$pos
  
  bins_exp <- GRanges(
    seqnames = seqnames(bins_gr),
    ranges   = IRanges(
      start = pmax(1L, start(bins_gr) - as.integer(window_bp)),
      end   = end(bins_gr) + as.integer(window_bp)
    )
  )
  mcols(bins_exp) <- mcols(bins_gr)
  
  hits <- findOverlaps(bins_exp, genes_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(tibble())
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  # distance from ORIGINAL bin to gene (0 if overlaps)
  a1 <- start(bins_gr)[q]; a2 <- end(bins_gr)[q]
  b1 <- start(genes_gr)[s]; b2 <- end(genes_gr)[s]
  dist_bp <- pmax(b1 - a2, a1 - b2, 0L)
  
  tibble(
    family   = as.character(mcols(bins_gr)$family[q]),
    marker   = as.character(mcols(bins_gr)$marker[q]),
    chr      = as.character(seqnames(bins_gr)[q]),
    bin_start= as.integer(a1),
    bin_end  = as.integer(a2),
    pos      = as.integer(mcols(bins_gr)$pos[q]),
    mlog10p  = as.numeric(mcols(bins_gr)$mlog10p[q]),
    p        = as.numeric(mcols(bins_gr)$p[q]),
    gene_id  = as.character(mcols(genes_gr)$gene_id[s]),
    gene_name= as.character(mcols(genes_gr)$gene_name[s]),
    gene_start = as.integer(b1),
    gene_end   = as.integer(b2),
    gene_strand= as.character(strand(genes_gr)[s]),
    dist_bp    = as.integer(dist_bp),
    in_gene    = dist_bp == 0L,
    in_window  = dist_bp <= as.integer(window_bp)
  ) %>%
    arrange(p, dist_bp, gene_id)
}

## ---------------------------
## 4) Summary table per marker (collapse multiple genes)
## ---------------------------
summarize_annotations <- function(ann_long, keep_all_overlaps = TRUE) {
  if (nrow(ann_long) == 0) return(tibble())
  
  df <- ann_long
  if (!keep_all_overlaps) {
    df <- df %>% group_by(marker) %>% slice_min(order_by = dist_bp, n = 1, with_ties = FALSE) %>% ungroup()
  }
  
  df %>%
    group_by(marker, chr, pos, bin_start, bin_end, family) %>%
    summarise(
      p       = min(p, na.rm = TRUE),
      mlog10p = max(mlog10p, na.rm = TRUE),
      n_genes = n_distinct(gene_id),
      genes   = paste(unique(gene_id), collapse = ","),
      gene_names = paste(unique(na.omit(gene_name)), collapse = ","),
      nearest_gene = gene_id[which.min(dist_bp)][1],
      nearest_gene_name = gene_name[which.min(dist_bp)][1],
      nearest_dist_bp = min(dist_bp, na.rm = TRUE),
      any_in_gene = any(in_gene, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(p)
}

## ---------------------------
## ONE CALL: do everything
## ---------------------------
annotate_introgression_gwas <- function(res_all,
                                        gff_file,
                                        window_bp = 25000L,
                                        chr_style = c("chr","nochr"),
                                        sig_thr_mlog10 = NULL,
                                        keep_all_overlaps = TRUE,
                                        keep_chr = 1:10) {
  chr_style <- match.arg(chr_style)
  
  gwas_bins <- prep_gwas_bins(res_all, chr_style = chr_style, sig_thr_mlog10 = sig_thr_mlog10)
  genes_gr  <- load_genes_from_gff(gff_file, chr_style = chr_style, keep_chr = keep_chr)
  
  ann_long <- annotate_bins_to_genes(gwas_bins, genes_gr, window_bp = window_bp)
  ann_sum  <- summarize_annotations(ann_long, keep_all_overlaps = keep_all_overlaps)
  
  list(
    gwas_bins = gwas_bins,
    genes_gr  = genes_gr,
    long      = ann_long,
    summary   = ann_sum
  )
}

## ============================================================
## YOUR EXACT RUN (copy/paste)
## ============================================================
gff_file <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

window_bp <- 50000L
sig_thr_mlog10 <- 3

ann <- annotate_introgression_gwas(
  res_all = res_all,
  gff_file = gff_file,
  window_bp = window_bp,
  chr_style = "chr",              # your GWAS uses chr1..chr10
  sig_thr_mlog10 = sig_thr_mlog10,
  keep_all_overlaps = TRUE,
  keep_chr = 1:10
)
str(res_all)

# results
ann$summary %>% dplyr::arrange(p) %>% head(20)

# save the results
write.csv(ann$summary, file = "GridLMM_Combined_DTS_GWAS_annotation_summary_100kwindow.csv", row.names = FALSE)


# Annotation family
ann_Zd <- annotate_introgression_gwas(
  res_all = res_by_family$Zd,
  gff_file = gff_file,
  window_bp = window_bp,
  chr_style = "chr",              # your GWAS uses chr1..chr10
  sig_thr_mlog10 = sig_thr_mlog10,
  keep_all_overlaps = TRUE,
  keep_chr = 1:10
)

# save the results
write.csv(ann_Zd$summary, file = "GridLMM_Zd_DTS_GWAS_annotation_summary_100kwindow.csv", row.names = FALSE)


ann_Zl <- annotate_introgression_gwas(
  res_all = res_by_family$Zl,
  gff_file = gff_file,
  window_bp = window_bp,
  chr_style = "chr",              # your GWAS uses chr1..chr10
  sig_thr_mlog10 = sig_thr_mlog10,
  keep_all_overlaps = TRUE,
  keep_chr = 1:10
)

# save the results
write.csv(ann_Zl$summary, file = "GridLMM_Zl_DTS_GWAS_annotation_summary_100kwindow.csv", row.names = FALSE)

ann_Zv <- annotate_introgression_gwas(
  res_all = res_by_family$Zv,
  gff_file = gff_file,
  window_bp = window_bp,
  chr_style = "chr",              # your GWAS uses chr1..chr10
  sig_thr_mlog10 = sig_thr_mlog10,
  keep_all_overlaps = TRUE,
  keep_chr = 1:10
)

# save the results
write.csv(ann_Zv$summary, file = "GridLMM_Zv_DTS_GWAS_annotation_summary_100kwindow.csv", row.names = FALSE)



## ============================================================
## Count SNPs table
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

## ------------------------------------------------------------
## 1) Clean / standardize GWAS table from res_all
##    Output columns: Marker, chr, CHR, BP, p, logp
## ------------------------------------------------------------
gwas <- res_all %>%
  mutate(
    Marker = dplyr::coalesce(.data$Marker, .data$X_ID),
    
    chr = dplyr::coalesce(
      .data$chr,
      str_extract(.data$Marker, "^[^:]+")
    ),
    
    pos = dplyr::coalesce(
      .data$pos,
      suppressWarnings(as.integer(str_extract(.data$Marker, "(?<=:)\\d+")))
    ),
    
    # numeric chromosome
    CHR = suppressWarnings(as.integer(str_extract(.data$chr, "\\d+"))),
    
    # GWAS position column name people expect
    BP = as.integer(.data$pos),
    
    # p-value (you already have both p and p_value_ML)
    p = dplyr::coalesce(.data$p, .data$p_value_ML),
    p = suppressWarnings(as.numeric(p)),
    p = pmin(pmax(p, 1e-300), 1),
    
    logp = -log10(p)
  ) %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(p)) %>%
  arrange(CHR, BP)

# quick sanity
dplyr::glimpse(gwas[, c("Marker","chr","CHR","BP","p","logp")])

## ------------------------------------------------------------
## 2) Count bins per chromosome (TOTAL)
## ------------------------------------------------------------
gwas %>%
  count(CHR, name = "n_total_bins") %>%
  arrange(CHR)

## ------------------------------------------------------------
## 3) Count bins per chromosome ABOVE thresholds
## ------------------------------------------------------------
gwas %>%
  group_by(CHR) %>%
  summarise(
    n_total = n(),
    n_ge_3  = sum(logp >= 3, na.rm = TRUE),
    n_ge_5  = sum(logp >= 5, na.rm = TRUE),
    n_ge_7  = sum(logp >= 7, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(CHR)

## ------------------------------------------------------------
## 4) Check your "binning step" per chr (should be ~ step_bp)
##    (median diff of consecutive BP within each chromosome)
## ------------------------------------------------------------
gwas %>%
  group_by(CHR) %>%
  summarise(
    median_step_bp = suppressWarnings(as.integer(median(diff(sort(unique(BP))), na.rm = TRUE))),
    n_unique_bp = n_distinct(BP),
    .groups = "drop"
  ) %>%
  arrange(CHR)




## ============================================================
## ============================================================
## ============================================================
## Manhattan plot for GridLMM introgression-bin GWAS (res_all)
## Optional: gene labels using ann$summary from annotate_introgression_gwas()
## ============================================================
## ============================================================
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

## ---------------------------
## 1) Prep GWAS table from res_all
## ---------------------------
prep_gridlmm_gwas <- function(res_all) {
  df <- as.data.frame(res_all)
  
  # pick marker col
  marker_col <- dplyr::coalesce(
    if ("Marker" %in% names(df)) "Marker" else NA_character_,
    if ("marker" %in% names(df)) "marker" else NA_character_,
    if ("X_ID"   %in% names(df)) "X_ID"   else NA_character_
  )
  if (is.na(marker_col)) stop("Need Marker/marker/X_ID in res_all")
  
  df <- df %>%
    mutate(
      Marker = .data[[marker_col]],
      chr = dplyr::coalesce(
        if ("chr" %in% names(df)) as.character(.data$chr) else NA_character_,
        str_extract(as.character(Marker), "^[^:]+")
      ),
      BP = dplyr::coalesce(
        if ("pos" %in% names(df)) suppressWarnings(as.integer(.data$pos)) else NA_integer_,
        suppressWarnings(as.integer(str_extract(as.character(Marker), "(?<=:)\\d+")))
      ),
      CHR = suppressWarnings(as.integer(str_extract(chr, "\\d+"))),
      p = dplyr::coalesce(
        if ("p" %in% names(df)) suppressWarnings(as.numeric(.data$p)) else NA_real_,
        if ("p_value_ML" %in% names(df)) suppressWarnings(as.numeric(.data$p_value_ML)) else NA_real_
      ),
      p = pmin(pmax(p, 1e-300), 1),
      logp = -log10(p)
    ) %>%
    filter(!is.na(CHR), !is.na(BP), !is.na(p)) %>%
    arrange(CHR, BP)
  
  df
}

## ---------------------------
## 2) Build chromosome offsets (from chr_len if provided; else from data)
## ---------------------------
make_chr_offsets <- function(gwas, chr_len = NULL) {
  if (!is.null(chr_len)) {
    # chr_len can be named (chr1..chr10) or numeric vector length 10
    if (!is.null(names(chr_len))) {
      # accept chr1..chr10 or 1..10
      nms <- names(chr_len)
      chr_num <- suppressWarnings(as.integer(str_extract(nms, "\\d+")))
      chr_tbl <- tibble(
        CHR = chr_num,
        chr_len = as.numeric(chr_len)
      ) %>%
        filter(!is.na(CHR)) %>%
        group_by(CHR) %>% summarise(chr_len = max(chr_len), .groups="drop") %>%
        arrange(CHR)
    } else {
      chr_tbl <- tibble(CHR = seq_along(chr_len), chr_len = as.numeric(chr_len)) %>%
        arrange(CHR)
    }
  } else {
    chr_tbl <- gwas %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP, na.rm = TRUE), .groups="drop") %>%
      arrange(CHR)
  }
  
  chr_tbl %>%
    mutate(
      offset = dplyr::lag(cumsum(chr_len), default = 0),
      center = offset + chr_len/2
    )
}

## ---------------------------
## 3) Manhattan plot (optional annotation labels)
## ---------------------------
plot_introgression_manhattan <- function(res_all,
                                         ann_summary = NULL,
                                         chr_len = NULL,
                                         thr_main = 7,
                                         thr_sugg = 5,
                                         label_thr = 7,
                                         label_top_n = 30,
                                         label_mode = c("nearest_gene", "genes", "marker"),
                                         title = "Introgression-bin GWAS (GridLMM) Manhattan plot",
                                         out_png = NULL,
                                         width = 14, height = 6, dpi = 300) {
  
  label_mode <- match.arg(label_mode)
  
  gwas <- prep_gridlmm_gwas(res_all)
  chr_tbl <- make_chr_offsets(gwas, chr_len = chr_len)
  
  gwas <- gwas %>%
    left_join(chr_tbl %>% select(CHR, offset), by = "CHR") %>%
    mutate(pos_cum = BP + offset)
  
  # alternate greys by chromosome parity
  gwas <- gwas %>% mutate(chr_parity = factor(CHR %% 2))
  
  # if annotation provided, join by Marker
  lab_df <- gwas
  if (!is.null(ann_summary) && nrow(ann_summary) > 0) {
    anns <- as.data.frame(ann_summary)
    # require marker column in annotation
    if (!("marker" %in% names(anns))) stop("ann_summary must have column 'marker' (use ann$summary)")
    anns$marker <- as.character(anns$marker)
    
    lab_df <- lab_df %>%
      left_join(
        anns %>% transmute(
          Marker = as.character(marker),
          nearest_gene = dplyr::coalesce(nearest_gene, NA_character_),
          nearest_gene_name = dplyr::coalesce(nearest_gene_name, NA_character_),
          genes = dplyr::coalesce(genes, NA_character_)
        ),
        by = "Marker"
      )
  } else {
    lab_df$nearest_gene <- NA_character_
    lab_df$nearest_gene_name <- NA_character_
    lab_df$genes <- NA_character_
  }
  
  # label text
  lab_df <- lab_df %>%
    mutate(
      label = dplyr::case_when(
        label_mode == "nearest_gene" & !is.na(nearest_gene) & nearest_gene != "" ~ nearest_gene,
        label_mode == "genes" & !is.na(genes) & genes != "" ~ genes,
        TRUE ~ as.character(Marker)
      )
    )
  
  # choose which points to label
  to_label <- lab_df %>%
    filter(logp >= label_thr) %>%
    arrange(p) %>%
    slice_head(n = label_top_n)
  
  p <- ggplot(gwas, aes(x = pos_cum, y = logp)) +
    geom_point(aes(color = chr_parity), alpha = 0.35, size = 1.2) +
    scale_color_manual(values = c("0" = "grey65", "1" = "grey35"), guide = "none") +
    
    geom_hline(yintercept = thr_sugg, linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = thr_main, linetype = "solid",  linewidth = 0.7) +
    
    scale_x_continuous(
      breaks = chr_tbl$center,
      labels = chr_tbl$CHR,
      expand = c(0.01, 0.01)
    ) +
    labs(
      title = title,
      x = "Chromosome",
      y = expression(-log[10](p))
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  if (nrow(to_label) > 0) {
    p <- p +
      geom_point(data = to_label, alpha = 0.9, size = 2.2) +
      ggrepel::geom_text_repel(
        data = to_label,
        aes(label = label),
        size = 4.2,
        min.segment.length = 0,
        box.padding = 0.35,
        point.padding = 0.25,
        seed = 42,
        max.overlaps = Inf
      )
  }
  
  if (!is.null(out_png)) {
    ggsave(out_png, p, width = width, height = height, dpi = dpi)
  }
  
  return(list(plot = p, gwas = gwas, chr_tbl = chr_tbl, labeled = to_label))
}

## ============================================================
## YOUR RUN
## ============================================================

# Option A: no gene labels (just bins)
out0 <- plot_introgression_manhattan(
  res_all = res_all,
  ann_summary = NULL,
  chr_len = NULL,              # or readRDS(".../chr_len.rds")
  thr_main = 7,
  thr_sugg = 5,
  label_thr = 7,
  label_top_n = 20,
  label_mode = "marker",
  out_png = "GridLMM_introgression_bins_manhattan.png"
)
quartz()
print(out0$plot)

# Option B: with gene labels (use your annotation)
# ann <- annotate_introgression_gwas(...)  # from your working code
out1 <- plot_introgression_manhattan(
  res_all = res_all,
  ann_summary = ann$summary,    # <-- gene annotation table
  chr_len = readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/genotype/chr_len.rds"),
  thr_main = 5,
  thr_sugg = 3,
  label_thr = 5,
  label_top_n = 20,
  label_mode = "nearest_gene",
  out_png = "GridLMM_introgression_bins_manhattan_labeled.png"
)
quartz()
print(out1$plot)


# # 6. Combined analysis: Scan each chromosome using GridLMM
# library(lme4)   # ensure lme4 is loaded as GridLMM depends on it
# results_combined <- list()  # to collect results from each chromosome
# for(chr in unique(marker_chr)) {
#   # Subset the genotype matrix or window matrix for this chromosome
#   # We'll demonstrate both: first using marker probabilities, then using windows.
#   
#   ## (a) Using single-marker probabilities:
#   markers_on_chr <- which(marker_chr == chr)
#   X_chr <- geno_mat[, markers_on_chr, drop=FALSE]
#   
#   # Run GridLMM GWAS for this chromosome's markers
#   # We'll include Family as a covariate (fixed effect) and a random effect for line with LOCO kinship
#   res_chr <- GridLMM_GWAS(
#     formula = Flowering_TIME_BLUE ~ Family + (1|Line),  # trait ~ Family + random polygenic
#     test_formula = ~ 1,     # test an intercept for each marker (i.e., additive effect of the marker)
#     reduced_formula = ~ 0,  # reduced model has no marker effect
#     data = data,
#     weights = NULL,
#     X = X_chr,
#     X_ID = "Line",
#     relmat = list(Line = kinship_list[[as.character(chr)]]),
#     method = "ML",    # use ML for LRT
#     #test = "LRT",     # perform likelihood ratio test for marker effect
#     mc.cores = 1      # set number of cores; adjust if parallelizing
#   )
#   results_combined[[as.character(chr)]] <- res_chr$results  # store results data frame for this chr
#   # Note: res_chr$results should contain columns like Marker, Chromosome, Position, p-value, etc., for each marker.
#   
#   ## (b) Using window-based dosages (optional alternative):
#   windows_on_chr <- names(window_markers)[startsWith(names(window_markers), paste0("Chr", chr, "_"))]
#   Xw_chr <- window_mat[, windows_on_chr, drop=FALSE]
#   res_chr_w <- GridLMM_GWAS(
#     formula = Flowering_TIME_BLUE ~ Family + (1|Line),
#     test_formula = ~ 1,
#     reduced_formula = ~ 0,
#     data = data,
#     X = Xw_chr,
#     X_ID = "Line",
#     relmat = list(Line = kinship_list[[as.character(chr)]]),
#     method = "ML",
#     #test = "LRT",
#     mc.cores = 1
#   )
#   # Store window results similarly (if needed):
#   # results_combined_windows[[as.character(chr)]] <- res_chr_w$results
# }
# # After this loop, we can combine results from all chromosomes:
# all_marker_results <- bind_rows(results_combined)
# # all_marker_results now holds the p-values (and optionally effect estimates) for every marker tested across the genome.
# # If GridLMM appended marker map info, it may have Chromosome and Position columns as well.
# # We can inspect the most significant hits:
# head(arrange(all_marker_results, p_value), 10)  # top 10 marker signals

# # 7. Family-wise analysis: loop over each family and perform the scan within that subset.
# results_by_family <- list()
# for(fam in levels(data$Family)) {
#   # Subset phenotype data for this family
#   data_fam <- data %>% filter(Family == fam)
#   # Subset kinship for this family (take the submatrix of the LOCO kinship for lines in this family)
#   kinship_fam <- lapply(kinship_list, function(K) {
#     K_sub <- K[data_fam$Line, data_fam$Line]
#     return(K_sub)
#   })
#   # Subset genotype probabilities and window matrices for these lines:
#   X_fam <- geno_mat[data_fam$Line, , drop=FALSE]
#   Xw_fam <- window_mat[data_fam$Line, , drop=FALSE]
#   
#   # Perform genome scan for this family (markers):
#   fam_marker_results <- list()
#   for(chr in unique(marker_chr)) {
#     X_chr <- X_fam[, which(marker_chr == chr), drop=FALSE]
#     res_fam_chr <- GridLMM_GWAS(
#       formula = Flowering_TIME_BLUE ~ 1 + (1|Line),  # within one family, no need for Family covariate (all same family)
#       test_formula = ~ 1,
#       reduced_formula = ~ 0,
#       data = data_fam,
#       X = X_chr,
#       X_ID = "Line",
#       relmat = list(Line = kinship_fam[[as.character(chr)]]),
#       method = "ML",
#       test = "LRT",
#       mc.cores = 1
#     )
#     fam_marker_results[[as.character(chr)]] <- res_fam_chr$results
#   }
#   fam_results_df <- bind_rows(fam_marker_results)
#   results_by_family[[fam]] <- fam_results_df
#   # We can also output or examine the top hits per family:
#   top_hit <- fam_results_df %>% arrange(p_value) %>% slice(1)
#   cat("Top QTL in family", fam, ": ", top_hit$Marker, "at p =", top_hit$p_value, "\n")
# }
# Now results_by_family is a list of result tables for each family.

# # 8. (Optional) Extracting effect size estimates:
# # GridLMM results typically include effect estimates (e.g., beta) for the marker in the tested model.
# # You might find a column like beta or Effect in res_chr$results. If not, one approach is to refit a model with the marker.
# # For example, for a significant marker M on chromosome 3:
# sig_marker <- all_marker_results %>% filter(Chromosome == 3) %>% arrange(p_value) %>% slice(1) %>% pull(Marker)
# # Fit the mixed model with this marker to get effect size:
# marker_vec <- geno_mat[, sig_marker]
# data$marker_X <- marker_vec[ data$Line ]  # align marker genotype to data
# final_model <- lmer(Flowering_TIME_BLUE ~ Family + marker_X + (1|Line), data=data, REML=TRUE)
# summary(final_model)$coefficients  # the marker_X coefficient here is the estimated effect (in trait units per allele).
# # (We use REML for final estimation since we are no longer comparing models for p-value.)
