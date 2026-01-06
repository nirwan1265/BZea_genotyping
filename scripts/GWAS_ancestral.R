## ============================================================
## FAMILY-WISE INTROGRESSION GWAS (B73 vs Teosinte) using tracts
##   - Input: x = list of (chr,start,end,state) per sample
##   - Genotype predictor at each locus: ancestry dosage (0/1/2)
##   - Control: SNP-based kinship matrix K_snp.rel (from imputed VCF)
##   - Runs separately for Zd, Zl, Zv, Zx
##   - Meta combine across families with ACAT
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rrBLUP)
})

## ---------------------------
## 0) EDIT THESE PATHS
## ---------------------------
x_rds <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BZea_Introgression_Finder/results_list_new_name.rds"

geno_dir <- "/Users/nirwantandukar/Documents/Research/data/BZea/genotype"
kid_file <- file.path(geno_dir, "K_snp.rel.id")
K_file   <- file.path(geno_dir, "K_snp.rel")

pheno_file <- "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv"
trait_col  <- "trait"


out_dir <- file.path(geno_dir, "Introgression_GWAS_familywise")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

families_keep <- c("Zd","Zl","Zv","Zx")

## ---------------------------
## 1) Helper: family code
## ---------------------------
get_family <- function(taxa) sub("\\..*$", "", taxa)

## ---------------------------
## 2) Read phenotype (2 columns)
## ---------------------------
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE) %>%
  dplyr::select(1,3) %>%
  setNames(c("Taxa", trait_col)) %>%
  mutate(
    Taxa = as.character(Taxa),
    !!trait_col := suppressWarnings(as.numeric(.data[[trait_col]]))
  ) %>%
  filter(!is.na(.data[[trait_col]])) %>%
  mutate(family = get_family(Taxa))

## ---------------------------
## 3) Read kinship K + IDs
## ---------------------------
kid <- fread(kid_file, header = FALSE)
if (kid[1, V1] == "#FID" || kid[1, V2] == "IID") kid <- kid[-1]

K <- as.matrix(fread(K_file, header = FALSE))
stopifnot(nrow(kid) == nrow(K), ncol(K) == nrow(K))

taxa_all <- as.character(kid$V2)
rownames(K) <- taxa_all
colnames(K) <- taxa_all

## Optional: strip trailing ".B" everywhere if it exists
stripB <- function(z) sub("\\.B$", "", z)
pheno$Taxa <- stripB(pheno$Taxa)
taxa_all   <- stripB(taxa_all)
rownames(K) <- taxa_all
colnames(K) <- taxa_all

## ---------------------------
## 4) Read introgression tracts list x
##    Each element: tibble with V1=chr, V2=start, V3=end, V4=state
## ---------------------------
x <- readRDS(x_rds)

# standardize + strip .B on sample names
names(x) <- stripB(names(x))

std_tbl <- function(tb) {
  tb %>%
    transmute(
      chr   = as.character(V1),
      start = as.integer(V2),
      end   = as.integer(V3),
      state = as.character(V4)
    ) %>%
    mutate(
      chr = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr))
    ) %>%
    arrange(chr, start, end)
}
x <- lapply(x, std_tbl)

## Keep only samples that exist in phenotype + kinship + tract list
taxa_ok <- Reduce(intersect, list(pheno$Taxa, rownames(K), names(x)))
pheno   <- pheno %>% filter(Taxa %in% taxa_ok)

## ---------------------------
## 5) Build ancestry markers genome-wide
##    We place a marker every step_bp and call ancestry at that position.
## ---------------------------
step_bp <- 250000L   # <-- change (100k, 250k, 500k). 250k ~ ~12k markers total.

# get chr lengths from max end across samples
chr_len <- bind_rows(x, .id="Taxa") %>%
  group_by(chr) %>%
  summarise(chr_len = max(end, na.rm=TRUE), .groups="drop") %>%
  arrange(as.integer(str_remove(chr,"chr")))

make_markers <- function(chr, L, step_bp) {
  pos <- seq.int(from = step_bp, to = L, by = step_bp)
  tibble(chr = chr, pos = as.integer(pos),
         marker = paste0(chr, ":", pos))
}

markers <- pmap_dfr(list(chr_len$chr, chr_len$chr_len, list(step_bp)),
                    ~ make_markers(..1, ..2, ..3))

markers_by_chr <- split(markers, markers$chr)

## map tract state -> dosage
state_to_dosage <- function(state_chr) {
  s <- tolower(state_chr)
  
  # robust matching
  out <- rep(NA_integer_, length(s))
  out[grepl("^b73$", s)] <- 0L
  out[grepl("introgression", s)] <- 2L
  out[grepl("het|hetero|b73/intro", s)] <- 1L
  
  # if your file uses only "B73" and "Introgression", you'll get 0/2 only.
  out
}

## call ancestry dosage at fixed positions for ONE chromosome
call_chr <- function(tb_chr, pos_vec) {
  # tb_chr: chr/start/end/state sorted
  starts <- tb_chr$start
  ends   <- tb_chr$end
  states <- tb_chr$state
  
  idx <- findInterval(pos_vec, starts)
  idx[idx == 0] <- 1L
  
  # sanity: positions should fall within interval; if not, NA
  bad <- pos_vec > ends[idx]
  st  <- states[idx]
  st[bad] <- NA_character_
  
  state_to_dosage(st)
}

## Build full ancestry genotype matrix: samples x markers
build_sample_geno <- function(tb) {
  v <- integer(nrow(markers))
  v[] <- NA_integer_
  
  for (cc in names(markers_by_chr)) {
    pos_vec <- markers_by_chr[[cc]]$pos
    tb_cc   <- tb %>% filter(chr == cc)
    
    if (nrow(tb_cc) == 0) {
      # if missing chr, leave NA
      next
    }
    
    dos <- call_chr(tb_cc, pos_vec)
    
    idx_global <- match(markers_by_chr[[cc]]$marker, markers$marker)
    v[idx_global] <- dos
  }
  
  v
}

message("Building ancestry matrix (this can take a bit)...")
G_anc <- do.call(rbind, lapply(x[pheno$Taxa], build_sample_geno))

#saveRDS(G_anc,"G_anc_25k.RDS")

rownames(G_anc) <- pheno$Taxa
colnames(G_anc) <- markers$marker

## ---------------------------
## 6) Family-wise GWAS using rrBLUP LMM
##    phenotype + kinship + 3 PCs from kinship
## ---------------------------
calc_pcs_from_K <- function(Kmat, npc=3) {
  eig <- eigen(Kmat, symmetric=TRUE)
  pcs <- eig$vectors[, seq_len(min(npc, ncol(eig$vectors))), drop=FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  pcs
}

maf_from_dosage <- function(v) {
  # allele freq p = mean(dosage)/2
  v <- v[!is.na(v)]
  if (length(v) == 0) return(0)
  p <- mean(v) / 2
  min(p, 1 - p)
}

run_family <- function(fam_code, npc = 3, min_maf = 0.01, max_miss = 0.2) {
  
  # 1) pick lines in this family that exist everywhere
  taxa_f <- pheno %>%
    dplyr::filter(family == fam_code) %>%
    dplyr::pull(Taxa) %>%
    intersect(rownames(G_anc)) %>%
    intersect(rownames(K))
  
  if (length(taxa_f) < 30) {
    message("Skipping ", fam_code, " (n<30): ", length(taxa_f))
    return(NULL)
  }
  
  # 2) phenotype (rrBLUP pheno: first col = gid, other cols = trait and fixed effects)
  Y_f <- pheno %>%
    dplyr::filter(Taxa %in% taxa_f) %>%
    dplyr::arrange(match(Taxa, taxa_f))
  
  K_f <- K[taxa_f, taxa_f, drop = FALSE]
  pcs <- calc_pcs_from_K(K_f, npc = npc)
  pcs <- as.data.frame(pcs)
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  
  ph_rr <- data.frame(
    gid   = taxa_f,
    trait = Y_f[[trait_col]],
    stringsAsFactors = FALSE
  )
  
  # 3) ancestry "genotypes": individuals x markers (0/1/2, mostly 0/2 for BC lines)
  M <- G_anc[taxa_f, , drop = FALSE]
  
  # filter markers by missingness + MAF (computed on 0/1/2 scale)
  miss  <- colMeans(is.na(M))
  keep1 <- miss <= max_miss
  
  maf <- apply(M[, keep1, drop = FALSE], 2, maf_from_dosage)
  keep2 <- maf >= min_maf
  
  M2 <- M[, keep1, drop = FALSE]
  M2 <- M2[, keep2, drop = FALSE]
  
  if (ncol(M2) < 100) {
    message("Skipping ", fam_code, " (too few markers after filters): ", ncol(M2))
    return(NULL)
  }
  
  message("Running ", fam_code, ": n=", nrow(M2), " markers=", ncol(M2))
  
  # 4) rrBLUP::GWAS geno format:
  #    data.frame: [marker, chr, pos, <one column per line>]
  #    and marker scores coded as {-1,0,1}
  #
  # Convert ancestry dosage to rrBLUP allele coding:
  #   B73 = -1, het = 0, teosinte = +1
  #   (if you only have 0 and 2, this becomes -1 and +1)
  Gscore <- M2
  Gscore[Gscore == 0] <- -1
  Gscore[Gscore == 1] <-  0
  Gscore[Gscore == 2] <-  1
  
  # transpose to markers x individuals
  Gmat <- t(Gscore)  # rows=markers, cols=lines
  colnames(Gmat) <- taxa_f  # ensure exact IDs
  
  # mean-impute remaining NAs per marker (row)
  if (anyNA(Gmat)) {
    mu <- rowMeans(Gmat, na.rm = TRUE)
    ij <- which(is.na(Gmat), arr.ind = TRUE)
    Gmat[ij] <- mu[ij[, 1]]
  }
  
  # marker map for rrBLUP (needs chr + pos; your 'markers' has these)
  map <- markers %>%
    dplyr::filter(marker %in% rownames(Gmat)) %>%
    dplyr::select(marker, chr, pos)
  
  # enforce identical marker order between map and Gmat
  map <- map[match(rownames(Gmat), map$marker), ]
  stopifnot(identical(map$marker, rownames(Gmat)))
  
  geno_df <- data.frame(
    marker = map$marker,
    chr    = map$chr,
    pos    = map$pos,
    as.data.frame(Gmat, check.names = FALSE),
    check.names = FALSE
  )
  
  # 5) final alignment sanity (this is the critical part)
  #    rrBLUP requires geno column names (col 4+) to match pheno gid and K rownames
  stopifnot(identical(colnames(geno_df)[4:ncol(geno_df)], ph_rr$gid))
  stopifnot(identical(rownames(K_f), ph_rr$gid))
  
  # 6) GWAS with PCs as fixed effects
  fixed_cols <- paste0("PC", seq_len(npc))
  
  gwas_raw <- rrBLUP::GWAS(
    pheno   = ph_rr,
    geno    = geno_df,
    K       = K_f,
    n.PC    = npc,     # <-- THIS replaces your fixed PC columns
    min.MAF = 0,       # you already filtered, keep 0 here
    P3D     = TRUE,
    plot    = FALSE
  )
  
  # gwas_raw returns -log10(p) in column "trait"
  gwas_tbl <- dplyr::as_tibble(gwas_raw) %>%
    dplyr::rename(mlog10p = trait) %>%
    dplyr::mutate(
      p = 10^(-mlog10p),
      family = fam_code
    ) %>%
    dplyr::arrange(p)
  
}


res_list <- lapply(families_keep, run_family)


names(res_list) <- families_keep
res_list <- res_list[!vapply(res_list, is.null, logical(1))]

res_all <- bind_rows(res_list)
str(res_all)
write.csv(res_all, file.path(out_dir, "introgression_GWAS_all_families_long_DTA_250k.csv"),
          row.names = FALSE)


## ============================================================
## Annotate your family-wise introgression GWAS hits to B73 genes
##   - Input: res_list (list of rrBLUP GWAS outputs per family)
##            OR the CSVs you wrote out per family
##   - Output: per-marker annotation + per-gene summary
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

## ---------------------------
## 0) Read / combine GWAS results
## ---------------------------
## If you already have res_list in memory:
gwas_all <- bind_rows(res_list, .id = "family") %>%
  #dplyr::rename(mlog10p = trait) %>%              # rrBLUP names the score column by the trait name
  mutate(p = 10^(-mlog10p)) %>%            # convert -log10(p) to p
  select(family, marker, chr, pos, p, mlog10p)

## If instead you saved CSVs:
# files <- list.files(out_dir, pattern="^introgression_GWAS_.*\\.csv$", full.names=TRUE)
# gwas_all <- bind_rows(lapply(files, read.csv), .id="file") %>%
#   mutate(family = sub("^introgression_GWAS_(Zd|Zl|Zv|Zx)\\.csv$", "\\1", basename(file))) %>%
#   rename(mlog10p = trait) %>%
#   mutate(p = 10^(-mlog10p)) %>%
#   select(family, marker, chr, pos, p, mlog10p)

## ---------------------------
## 1) Load B73 gene annotation (GFF3) -> genes GRanges
## ---------------------------
# genes only
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# helper: get gene ID robustly (GFF3 attribute names vary)
get_gene_id <- function(gr) {
  m <- mcols(gr)
  id <- NULL
  if ("ID" %in% names(m)) id <- as.character(m$ID)
  if (is.null(id) && "gene_id" %in% names(m)) id <- as.character(m$gene_id)
  if (is.null(id) && "Name" %in% names(m)) id <- as.character(m$Name)
  if (is.null(id)) id <- rep(NA_character_, length(gr))
  id
}

# helper: get gene name/symbol if present
get_gene_name <- function(gr) {
  m <- mcols(gr)
  if ("Name" %in% names(m)) return(as.character(m$Name))
  if ("gene" %in% names(m)) return(as.character(m$gene))
  rep(NA_character_, length(gr))
}

gene_id <- get_gene_id(genes)
gene_id <- sub("^gene:", "", gene_id)   # common in RefGen_v5

gene_name <- get_gene_name(genes)

mcols(genes)$gene_id   <- gene_id
mcols(genes)$gene_name <- gene_name

# harmonize chromosome naming to match your markers (chr1..chr10)
seqlevels(genes) <- ifelse(grepl("^chr", seqlevels(genes)),
                           seqlevels(genes),
                           paste0("chr", seqlevels(genes)))

# drop non-primary seqlevels if needed
genes <- keepSeqlevels(genes, paste0("chr", 1:10), pruning.mode = "coarse")

## ---------------------------
## 2) Annotation function (nearest gene + distance, plus "within window" flag)
## ---------------------------
annotate_markers_to_genes <- function(gwas_df, genes_gr, window_bp = 25000L) {
  
  # marker GRanges
  mk <- GRanges(
    seqnames = gwas_df$chr,
    ranges   = IRanges(start = gwas_df$pos, end = gwas_df$pos)
  )
  mcols(mk)$marker <- gwas_df$marker
  
  # nearest gene (distance 0 means inside gene or at boundary)
  dn <- distanceToNearest(mk, genes_gr, ignore.strand = TRUE)
  q  <- queryHits(dn)
  s  <- subjectHits(dn)
  dist_bp <- mcols(dn)$distance
  
  ann <- gwas_df
  ann$nearest_gene_id   <- NA_character_
  ann$nearest_gene_name <- NA_character_
  ann$gene_chr          <- NA_character_
  ann$gene_start        <- NA_integer_
  ann$gene_end          <- NA_integer_
  ann$gene_strand       <- NA_character_
  ann$dist_to_gene_bp   <- NA_integer_
  ann$in_window         <- FALSE
  ann$in_gene           <- FALSE
  
  ann$nearest_gene_id[q]   <- mcols(genes_gr)$gene_id[s]
  ann$nearest_gene_name[q] <- mcols(genes_gr)$gene_name[s]
  ann$gene_chr[q]          <- as.character(seqnames(genes_gr)[s])
  ann$gene_start[q]        <- start(genes_gr)[s]
  ann$gene_end[q]          <- end(genes_gr)[s]
  ann$gene_strand[q]       <- as.character(strand(genes_gr)[s])
  ann$dist_to_gene_bp[q]   <- as.integer(dist_bp)
  
  ann$in_gene   <- ann$dist_to_gene_bp == 0L
  ann$in_window <- ann$dist_to_gene_bp <= window_bp
  
  ann
}

## ---------------------------
## 3) Choose your significance threshold + annotate
## ---------------------------
sig_thr_mlog10 <- 3   # e.g. -log10(p) >= 7  (p <= 1e-7)
window_bp <- 25000L   # +/- 25 kb around nearest gene

gwas_sig <- gwas_all %>%
  filter(mlog10p >= sig_thr_mlog10)

ann_sig <- annotate_markers_to_genes(gwas_sig, genes, window_bp = window_bp)

# write per-marker annotated hits
write.csv(ann_sig,
          file.path(out_dir, paste0("introgression_GWAS_hits_annotated_mlog10ge", sig_thr_mlog10, "_win", window_bp, "_250k_DTA.csv")),
          row.names = FALSE)


## ---------------------------
## 4) Per-gene summary (best marker per gene per family)
## ---------------------------
gene_summary <- ann_sig %>%
  filter(in_window) %>%   # keep only if within window of nearest gene
  group_by(family, nearest_gene_id, nearest_gene_name) %>%
  summarise(
    best_p       = min(p, na.rm = TRUE),
    best_mlog10p = max(mlog10p, na.rm = TRUE),
    best_marker  = marker[which.min(p)][1],
    best_chr     = chr[which.min(p)][1],
    best_pos     = pos[which.min(p)][1],
    dist_bp      = dist_to_gene_bp[which.min(p)][1],
    .groups = "drop"
  ) %>%
  arrange(best_p)

write.csv(gene_summary,
          file.path(out_dir, paste0("introgression_GWAS_gene_summary_mlog10ge", sig_thr_mlog10, "_win", window_bp, ".csv")),
          row.names = FALSE)

## Done: ann_sig = per-marker annotation
##       gene_summary = per-gene best-hit summary






## ============================================================
## HYBRID INTROGRESSION GWAS (most robust)
##   Predictor: 25kb BIN FRACTION introgressed (0..1) -> scaled to [-1,+1]
##   Kinship:    SNP-based K_snp.rel (from your imputed VCF)
##   Runs:       separately within each family (Zd, Zl, Zv, Zx)
##   Output:     per-family GWAS tables + long combined table
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(IRanges)
  library(rrBLUP)
})

## ---------------------------
## 0) EDIT PATHS
## ---------------------------
x_rds <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BZea_Introgression_Finder/results_list_new_name.rds"

geno_dir <- "/Users/nirwantandukar/Documents/Research/data/BZea/genotype"
kid_file <- file.path(geno_dir, "K_snp.rel.id")
K_file   <- file.path(geno_dir, "K_snp.rel")

pheno_file <- "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv"
trait_col  <- "trait"

out_dir <- file.path(geno_dir, "Introgression_GWAS_familywise_HYBRID_binFrac_Ksnp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

families_keep <- c("Zd","Zl","Zv","Zx")

## ---------------------------
## 1) Helpers
## ---------------------------
stripB <- function(z) sub("\\.B$", "", z)
get_family <- function(taxa) sub("\\..*$", "", taxa)

std_tbl <- function(tb) {
  tb %>%
    transmute(
      chr   = as.character(V1),
      start = as.integer(V2),
      end   = as.integer(V3),
      state = as.character(V4)
    ) %>%
    mutate(chr = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr))) %>%
    arrange(chr, start, end)
}

acat <- function(p, w = NULL, min_p = 1e-15) {
  p <- suppressWarnings(as.numeric(p))
  p <- p[!is.na(p)]
  if (length(p) == 0) return(NA_real_)
  p <- pmin(pmax(p, min_p), 1 - min_p)
  if (is.null(w)) w <- rep(1 / length(p), length(p))
  w <- w / sum(w)
  tstat <- sum(w * tan((0.5 - p) * pi))
  0.5 - atan(tstat) / pi
}

## ---------------------------
## 2) Load phenotype
##    NOTE: you had select(1,3) in your pasted script — keep that if trait is col 3.
## ---------------------------
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE) %>%
  dplyr::select(1, 3) %>%                  # <-- change to (1,2) if your trait is column 2
  setNames(c("Taxa", trait_col)) %>%
  mutate(
    Taxa = stripB(as.character(Taxa)),
    !!trait_col := suppressWarnings(as.numeric(.data[[trait_col]]))
  ) %>%
  filter(!is.na(.data[[trait_col]])) %>%
  mutate(family = get_family(Taxa))

## ---------------------------
## 3) Load SNP kinship (K_snp.rel) + IDs
## ---------------------------
kid <- fread(kid_file, header = FALSE)
if (kid[1, V1] == "#FID" || kid[1, V2] == "IID") kid <- kid[-1]

K <- as.matrix(fread(K_file, header = FALSE))
stopifnot(nrow(kid) == nrow(K), ncol(K) == nrow(K))

taxa_all <- stripB(as.character(kid$V2))
rownames(K) <- taxa_all
colnames(K) <- taxa_all

## ---------------------------
## 4) Load introgression tracts list x
## ---------------------------
x <- readRDS(x_rds)
names(x) <- stripB(names(x))
x <- lapply(x, std_tbl)

## Keep only taxa that exist in all pieces
taxa_ok <- Reduce(intersect, list(pheno$Taxa, rownames(K), names(x)))
pheno   <- pheno %>% filter(Taxa %in% taxa_ok)
x       <- x[taxa_ok]
K       <- K[taxa_ok, taxa_ok, drop = FALSE]

## ---------------------------
## 5) Build bins genome-wide + BIN-FRACTION matrix A (samples x bins)
## ---------------------------
make_bins <- function(x_list, window_bp = 25000L) {
  chr_len <- bind_rows(x_list, .id = "Taxa") %>%
    group_by(chr) %>%
    summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") %>%
    arrange(as.integer(str_remove(chr, "chr")))
  
  make_bins_one_chr <- function(chr, L, w) {
    starts <- seq.int(1L, L, by = w)
    tibble(
      chr = chr,
      bin_start = starts,
      bin_end   = pmin(starts + w - 1L, L)
    )
  }
  
  bins <- pmap_dfr(list(chr_len$chr, chr_len$chr_len, list(window_bp)),
                   ~ make_bins_one_chr(..1, ..2, ..3)) %>%
    mutate(
      bin_width = bin_end - bin_start + 1L,
      bin_id    = paste0(chr, ":", bin_start, "-", bin_end)
    )
  
  bins_by_chr <- split(bins, bins$chr)
  bins_ir_by_chr <- lapply(bins_by_chr, \(df) IRanges(df$bin_start, df$bin_end))
  
  list(bins = bins, bins_by_chr = bins_by_chr, bins_ir_by_chr = bins_ir_by_chr)
}

build_bin_fraction_matrix <- function(x_list, bins_obj, intro_label = "introgression") {
  bins <- bins_obj$bins
  bins_by_chr <- bins_obj$bins_by_chr
  bins_ir_by_chr <- bins_obj$bins_ir_by_chr
  
  sample_to_vec <- function(tb) {
    v <- numeric(nrow(bins))
    
    tb_intro <- tb %>% filter(tolower(state) == intro_label)
    if (nrow(tb_intro) == 0) return(v)
    
    for (cc in unique(tb_intro$chr)) {
      seg <- tb_intro %>% filter(chr == cc)
      if (!cc %in% names(bins_by_chr)) next
      
      seg_ir <- IRanges(seg$start, seg$end)
      bin_df <- bins_by_chr[[cc]]
      bin_ir <- bins_ir_by_chr[[cc]]
      
      hits <- findOverlaps(bin_ir, seg_ir)
      if (length(hits) == 0) next
      
      inter <- pintersect(bin_ir[queryHits(hits)], seg_ir[subjectHits(hits)])
      w_int <- width(inter)
      
      local_bin_ids <- bin_df$bin_id[queryHits(hits)]
      global_idx <- match(local_bin_ids, bins$bin_id)
      
      tmp <- tapply(w_int, global_idx, sum)
      idx <- as.integer(names(tmp))
      v[idx] <- as.numeric(tmp) / bins$bin_width[idx]
    }
    v
  }
  
  A <- do.call(rbind, lapply(x_list, sample_to_vec))
  rownames(A) <- names(x_list)
  colnames(A) <- bins$bin_id
  
  # drop bins with zero variance across all samples (saves time)
  keep_var <- apply(A, 2, sd, na.rm = TRUE) > 0
  A <- A[, keep_var, drop = FALSE]
  bins_use <- bins %>% filter(bin_id %in% colnames(A))
  
  list(A = A, bins_use = bins_use)
}


## ---------------------------
## 6) Run family-wise rrBLUP GWAS
##    Predictor: bin fraction -> scaled [-1,+1]
##    Kinship:   subset of SNP K_snp.rel
## ---------------------------
## ---- drop-in helper: standardize rrBLUP GWAS output ----
standardize_rrblup_gwas <- function(g_raw, fam_code, bins_use, trait_col = "trait") {
  
  # rrBLUP sometimes returns a list (named by trait); grab the first/target element
  if (is.list(g_raw) && !is.data.frame(g_raw)) {
    if (trait_col %in% names(g_raw)) {
      g_raw <- g_raw[[trait_col]]
    } else if ("trait" %in% names(g_raw)) {
      g_raw <- g_raw[["trait"]]
    } else {
      g_raw <- g_raw[[1]]
    }
  }
  
  g_tbl <- tibble::as_tibble(g_raw)
  
  # marker column name can vary
  mcol <- intersect(c("Marker","marker","SNP","snp","ID","id"), names(g_tbl))[1]
  if (is.na(mcol)) {
    stop("No marker column found in rrBLUP output. Columns: ",
         paste(names(g_tbl), collapse = ", "))
  }
  
  # score column (rrBLUP stores -log10(p) under the trait name)
  score_col <- intersect(names(g_tbl), c(trait_col, "trait", "y"))[1]
  
  if (is.na(score_col)) {
    # fallback: pick the first numeric column that isn't known metadata
    known <- c(mcol, "Chrom","chrom","Chr","chr","Pos","pos","Position","position",
               "Effect","effect","P.value","p.value","R2","r2","MAF","maf")
    cand <- setdiff(names(g_tbl), known)
    cand <- cand[sapply(g_tbl[cand], is.numeric)]
    score_col <- cand[1]
  }
  
  if (is.na(score_col)) {
    stop("Could not identify rrBLUP score column (-log10p). Columns: ",
         paste(names(g_tbl), collapse = ", "))
  }
  
  # Build clean output
  out <- g_tbl %>%
    dplyr::rename(marker = !!rlang::sym(mcol)) %>%
    dplyr::mutate(
      family  = fam_code,
      mlog10p = as.numeric(.data[[score_col]]),
      p       = 10^(-mlog10p)
    ) %>%
    dplyr::left_join(
      bins_use %>% dplyr::select(bin_id, chr, bin_start, bin_end),
      by = c("marker" = "bin_id")
    ) %>%
    dplyr::arrange(p)
  
  out
}
## ---- family-wise hybrid GWAS: bin-fraction predictor + SNP kinship (K_snp.rel) ----
run_family_binFrac_Ksnp <- function(fam_code,
                                    pheno, A, K_snp, bins_use,
                                    trait_col = "trait",
                                    npc = 3,
                                    min_maf = 0.01,
                                    max_miss = 0.2) {
  
  taxa_f <- pheno %>%
    dplyr::filter(family == fam_code) %>%
    dplyr::pull(Taxa) %>%
    intersect(rownames(A)) %>%
    intersect(rownames(K_snp))
  
  if (length(taxa_f) < 30) {
    message("Skipping ", fam_code, " (n<30): ", length(taxa_f))
    return(NULL)
  }
  
  # phenotype
  Y_f <- pheno %>%
    dplyr::filter(Taxa %in% taxa_f) %>%
    dplyr::arrange(match(Taxa, taxa_f))
  
  ph_rr <- data.frame(
    gid   = taxa_f,
    trait = suppressWarnings(as.numeric(Y_f[[trait_col]])),
    stringsAsFactors = FALSE
  )
  
  # kinship subset (SNP-based)
  K_f <- K_snp[taxa_f, taxa_f, drop = FALSE]
  
  # ancestry matrix subset: samples x bins (0..1)
  M <- A[taxa_f, , drop = FALSE]
  
  # missingness filter (rare unless NA exists)
  miss  <- colMeans(is.na(M))
  keep1 <- miss <= max_miss
  
  # scale to [-1,+1] where 0 -> -1 (all B73), 1 -> +1 (fully introgressed)
  M_sc <- 2 * M - 1
  
  maf_bin <- function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) return(0)
    d <- v + 1              # [-1,+1] -> [0,2]
    p <- mean(d) / 2
    min(p, 1 - p)
  }
  
  maf   <- apply(M_sc[, keep1, drop = FALSE], 2, maf_bin)
  keep2 <- maf >= min_maf
  
  M2 <- M_sc[, keep1, drop = FALSE]
  M2 <- M2[, keep2, drop = FALSE]
  
  if (ncol(M2) < 100) {
    message("Skipping ", fam_code, " (too few bins after filters): ", ncol(M2))
    return(NULL)
  }
  
  message("Running ", fam_code, ": n=", length(taxa_f), " bins=", ncol(M2))
  
  # rrBLUP geno format: rows=markers, cols=lines
  Gmat <- t(M2)     # rows=bins, cols=lines
  colnames(Gmat) <- taxa_f
  
  # impute NAs per bin (row)
  if (anyNA(Gmat)) {
    mu <- rowMeans(Gmat, na.rm = TRUE)
    ij <- which(is.na(Gmat), arr.ind = TRUE)
    Gmat[ij] <- mu[ij[, 1]]
  }
  
  # map for rrBLUP (Chrom integer, Pos midpoint)
  map <- bins_use %>%
    dplyr::filter(bin_id %in% rownames(Gmat)) %>%
    dplyr::mutate(
      Chrom = as.integer(stringr::str_remove(chr, "^chr")),
      Pos   = as.integer((bin_start + bin_end) / 2)
    ) %>%
    dplyr::select(bin_id, Chrom, Pos)
  
  map <- map[match(rownames(Gmat), map$bin_id), ]
  stopifnot(identical(map$bin_id, rownames(Gmat)))
  
  geno_df <- data.frame(
    marker = map$bin_id,
    chr    = map$Chrom,
    pos    = map$Pos,
    as.data.frame(Gmat, check.names = FALSE),
    check.names = FALSE
  )
  
  # sanity alignment
  stopifnot(identical(colnames(geno_df)[4:ncol(geno_df)], ph_rr$gid))
  stopifnot(identical(rownames(K_f), ph_rr$gid))
  
  # GWAS
  g_raw <- rrBLUP::GWAS(
    pheno   = ph_rr,
    geno    = geno_df,
    K       = K_f,
    n.PC    = npc,
    min.MAF = 0,
    P3D     = TRUE,
    plot    = FALSE
  )
  
  # standardize to (marker, family, mlog10p, p, chr/bin coords)
  out <- standardize_rrblup_gwas(g_raw, fam_code = fam_code, bins_use = bins_use, trait_col = "trait")
  out
}



## ---------------------------
## 7) MASTER: run everything in one go
## ---------------------------
run_hybrid_introgression_gwas <- function(window_bp = 25000L,
                                          npc = 3,
                                          min_maf = 0.01,
                                          max_miss = 0.2) {
  
  # Build bins + A
  bins_obj <- make_bins(x, window_bp = window_bp)
  built    <- build_bin_fraction_matrix(x, bins_obj, intro_label = "introgression")
  A        <- built$A
  bins_use <- built$bins_use
  
  # run each family
  res_list <- lapply(families_keep, run_family_binFrac_Ksnp,
                     pheno = pheno,
                     A = A,
                     K_snp = K,
                     bins_use = bins_use,
                     trait_col = trait_col,
                     npc = npc,
                     min_maf = min_maf,
                     max_miss = max_miss)
  
  names(res_list) <- families_keep
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  
  res_all <- bind_rows(res_list)
  
  # save
  tag <- paste0("bin", window_bp, "_npc", npc, "_maf", min_maf, "_miss", max_miss)
  saveRDS(res_list, file.path(out_dir, paste0("res_list_", tag, ".rds")))
  write.csv(res_all, file.path(out_dir, paste0("introgression_GWAS_hybrid_all_", tag, ".csv")),
            row.names = FALSE)
  
  list(res_list = res_list, res_all = res_all, A = A, bins_use = bins_use)
}

## ---------------------------
## 8) RUN IT
## ---------------------------
ans <- run_hybrid_introgression_gwas(
  window_bp = 25000L,  # try 25k, 50k, 100k
  npc       = 3,
  min_maf   = 0.01,
  max_miss  = 0.2
)

res_list <- ans$res_list
res_all  <- ans$res_all


## ============================================================
## Annotate BIN-based introgression GWAS (res_all) to B73 genes
##   - Works with your res_all columns: marker, chr.x, pos, family, mlog10p, p, chr.y, bin_start, bin_end
##   - Outputs:
##       1) ann_bins: per-bin annotations (all overlapping genes within ±window)
##       2) gene_summary: per-gene best hit per family
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

## ---------------------------
## 0) Inputs
## ---------------------------
# Your GWAS results tibble:
# res_all <- ...

gff_file <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

window_bp <- 25000L         # +/- around bin when annotating genes
sig_thr_mlog10 <- 3         # optional filter threshold (set NULL to keep all)

out_dir <- "/Users/nirwantandukar/Documents/Research/data/BZea/genotype/Introgression_GWAS_familywise"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------
## 1) Harmonize res_all columns
## ---------------------------
gwas_bins <- res_all %>%
  mutate(
    chr = ifelse(!is.na(chr.y), as.character(chr.y),
                 ifelse(!is.na(chr.x), paste0("chr", chr.x), NA_character_)),
    bin_start = as.integer(bin_start),
    bin_end   = as.integer(bin_end),
    pos = as.integer(pos),
    mlog10p = as.numeric(mlog10p),
    p = as.numeric(p),
    family = as.character(family),
    marker = as.character(marker)
  ) %>%
  select(marker, family, chr, pos, bin_start, bin_end, mlog10p, p) %>%
  filter(!is.na(chr), !is.na(bin_start), !is.na(bin_end))

if (!is.null(sig_thr_mlog10)) {
  gwas_bins <- gwas_bins %>% filter(mlog10p >= sig_thr_mlog10)
}

## ---------------------------
## 2) Load genes from GFF3
## ---------------------------
ref_gr <- rtracklayer::import(gff_file)
genes  <- ref_gr[mcols(ref_gr)$type == "gene"]

# robust gene id/name getters
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

# force chr naming like chr1..chr10
seqlevels(genes) <- ifelse(grepl("^chr", seqlevels(genes)),
                           seqlevels(genes),
                           paste0("chr", seqlevels(genes)))
genes <- keepSeqlevels(genes, paste0("chr", 1:10), pruning.mode = "coarse")

## ---------------------------
## 3) Annotation function: overlap within expanded window
## ---------------------------
annotate_bins_to_genes <- function(gwas_df, genes_gr, window_bp = 25000L) {
  
  bins_gr <- GRanges(
    seqnames = gwas_df$chr,
    ranges   = IRanges(start = gwas_df$bin_start, end = gwas_df$bin_end)
  )
  mcols(bins_gr)$marker  <- gwas_df$marker
  mcols(bins_gr)$family  <- gwas_df$family
  mcols(bins_gr)$mlog10p <- gwas_df$mlog10p
  mcols(bins_gr)$p       <- gwas_df$p
  mcols(bins_gr)$pos     <- gwas_df$pos
  
  # expand bins by +/- window
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
  
  # distance from ORIGINAL bin to gene (0 if overlap)
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
## 4) Run annotation
## ---------------------------
ann_bins <- annotate_bins_to_genes(gwas_bins, genes, window_bp = window_bp)

write.csv(
  ann_bins,
  file.path(out_dir, paste0("mixed_GWAS_family_bins", window_bp,
                            if (!is.null(sig_thr_mlog10)) paste0("_mlog10ge", sig_thr_mlog10) else "",
                            ".csv")),
  row.names = FALSE
)

## ---------------------------
## 5) Per-gene best-hit summary (by family)
## ---------------------------
gene_summary <- ann_bins %>%
  filter(in_window) %>%
  group_by(family, gene_id, gene_name) %>%
  summarise(
    best_p       = min(p, na.rm = TRUE),
    best_mlog10p = max(mlog10p, na.rm = TRUE),
    best_marker  = marker[which.min(p)][1],
    best_chr     = chr[which.min(p)][1],
    best_bin_start = bin_start[which.min(p)][1],
    best_bin_end   = bin_end[which.min(p)][1],
    best_pos       = pos[which.min(p)][1],
    best_dist_bp   = dist_bp[which.min(p)][1],
    .groups = "drop"
  ) %>%
  arrange(best_p)

write.csv(
  gene_summary,
  file.path(out_dir, paste0("introgression_bins_GWAS_gene_summary_win", window_bp,
                            if (!is.null(sig_thr_mlog10)) paste0("_mlog10ge", sig_thr_mlog10) else "",
                            ".csv")),
  row.names = FALSE
)

ann_bins
gene_summary








## ============================================================
## DTS vs DTA: correlations + plots (overall + by family)
##  - Viridis colors, Nature-friendly styling
##  - Inputs: pheno_file OR existing pheno data.frame
##  - Expects columns: new_genotype, DTS, DTA
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

## ----------------------------
## 0) Load phenotype (edit pheno_file if needed)
## ----------------------------
# pheno_file <- "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv"
# pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)

stopifnot(all(c("new_genotype","DTS","DTA") %in% names(pheno)))

stripB <- function(z) sub("\\.B$", "", z)
get_family <- function(taxa) sub("\\..*$", "", taxa)

ph <- pheno %>%
  mutate(
    new_genotype = stripB(as.character(new_genotype)),
    family = get_family(new_genotype),
    DTS = suppressWarnings(as.numeric(DTS)),
    DTA = suppressWarnings(as.numeric(DTA))
  ) %>%
  filter(is.finite(DTS), is.finite(DTA))

## ----------------------------
## 1) Correlation (overall + by family)
## ----------------------------
safe_cor_test <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3) return(tibble(n = n, r = NA_real_, p = NA_real_))
  if (sd(x) == 0 || sd(y) == 0) return(tibble(n = n, r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x, y, method = method))
  tibble(n = n, r = unname(ct$estimate), p = ct$p.value)
}

cor_all <- safe_cor_test(ph$DTS, ph$DTA) %>%
  mutate(group = "All")

cor_by_family <- ph %>%
  group_by(family) %>%
  summarise(safe_cor_test(DTS, DTA), .groups = "drop") %>%
  mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
  arrange(p)

print(cor_all)
print(cor_by_family)

## ----------------------------
## 2) Scatter plots: overall + by family
## ----------------------------
label_all <- cor_all %>%
  mutate(
    label = sprintf("Pearson r = %.2f\np = %.2g\nn = %d", r, p, n),
    x = Inf, y = -Inf
  )

label_fam <- cor_by_family %>%
  mutate(
    label = ifelse(
      is.na(r),
      sprintf("n = %d\nr = NA", n),
      sprintf("Pearson r = %.2f\np = %.2g\nn = %d", r, p, n)
    ),
    x = Inf, y = -Inf
  )

p_scatter_all <- ggplot(ph, aes(x = DTS, y = DTA)) +
  geom_point(aes(color = family), alpha = 0.75, size = 1.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  geom_text(
    data = label_all,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = -0.15, size = 3.3
  ) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  labs(
    x = "Days to Silking (DTS)",
    y = "Days to Anthesis (DTA)",
    color = "Family"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold")
  )

p_scatter_byfam <- ggplot(ph, aes(x = DTS, y = DTA)) +
  geom_point(aes(color = family), alpha = 0.75, size = 1.5, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "black") +
  facet_wrap(~family, scales = "free", ncol = 2) +
  geom_text(
    data = label_fam,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = -0.15, size = 3.0
  ) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  labs(
    x = "Days to Silking (DTS)",
    y = "Days to Anthesis (DTA)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

## ----------------------------
## 3) Box plots: overall + by family
## ----------------------------
ph_long <- ph %>%
  select(new_genotype, family, DTS, DTA) %>%
  pivot_longer(cols = c(DTS, DTA), names_to = "Trait", values_to = "Value") %>%
  mutate(Trait = factor(Trait, levels = c("DTS","DTA")))

p_box_all <- ggplot(ph_long, aes(x = Trait, y = Value, fill = Trait)) +
  geom_boxplot(width = 0.62, outlier.alpha = 0.25, linewidth = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.28, size = 0.9, show.legend = FALSE) +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  labs(x = NULL, y = "Days") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

p_box_byfam <- ggplot(ph_long, aes(x = family, y = Value, fill = Trait)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.alpha = 0.25,
    linewidth = 0.6
  ) +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  labs(x = "Family", y = "Days", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )


quartz()
p_scatter_all

quartz()
p_scatter_byfam

quartz()
p_box_all

quartz()
p_box_byfam
## ----------------------------
## 4) Save (PDF for vector; TIFF for Nature-style raster)
## ----------------------------
fig_dir <- file.path(getwd(), "Figures_DTS_DTA")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(fig_dir, "DTS_DTA_scatter_all.pdf"),   p_scatter_all,  width = 6.6, height = 4.6)
ggsave(file.path(fig_dir, "DTS_DTA_scatter_byFamily.pdf"), p_scatter_byfam, width = 7.2, height = 6.6)
ggsave(file.path(fig_dir, "DTS_DTA_box_a1`ll.pdf"),       p_box_all,      width = 4.8, height = 4.6)
ggsave(file.path(fig_dir, "DTS_DTA_box_byFamily.pdf"),  p_box_byfam,    width = 6.8, height = 4.8)

ggsave(file.path(fig_dir, "DTS_DTA_scatter_all.tiff"),     p_scatter_all,  width = 6.6, height = 4.6, dpi = 600, compression = "lzw")
ggsave(file.path(fig_dir, "DTS_DTA_scatter_byFamily.tiff"),p_scatter_byfam,width = 7.2, height = 6.6, dpi = 600, compression = "lzw")
ggsave(file.path(fig_dir, "DTS_DTA_box_all.tiff"),         p_box_all,      width = 4.8, height = 4.6, dpi = 600, compression = "lzw")
ggsave(file.path(fig_dir, "DTS_DTA_box_byFamily.tiff"),    p_box_byfam,    width = 6.8, height = 4.8, dpi = 600, compression = "lzw")

## View in RStudio:
p_scatter_all
p_scatter_byfam
p_box_all
p_box_byfam
