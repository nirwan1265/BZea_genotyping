## ============================================================
## FLEXIBLE INTROGRESSION GWAS (step_bp grid) with options:
##   - ancestry_mode:
##       "hard"         : state -> dosage {0,1,2}
##       "soft_state"   : state -> P(teo) -> dosage = 2*P(teo)
##       "posterior_mat": use posterior matrix (P(teo) or dosage) aligned to markers
##   - add_seg_covariates: add global introgression burden covariates to pheno
##   - weight_by_seglen  : weight marker predictor by segment length (local)
##   - family-wise GWAS with rrBLUP::GWAS and SNP-kinship K_snp.rel
##   - optional ACAT meta across families per marker
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rrBLUP)
})

## ---------------------------
## Small helpers
## ---------------------------
stripB <- function(z) sub("\\.B$", "", z)
get_family <- function(taxa) sub("\\..*$", "", taxa)

calc_pcs_from_K <- function(Kmat, npc = 3) {
  eig <- eigen(Kmat, symmetric = TRUE)
  pcs <- eig$vectors[, seq_len(min(npc, ncol(eig$vectors))), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  pcs
}

maf_from_dosage <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(0)
  p <- mean(v) / 2
  min(p, 1 - p)
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
## Standardize tract tibble
## ---------------------------
std_tbl <- function(tb) {
  tb %>%
    transmute(
      chr   = as.character(V1),
      start = as.integer(V2),
      end   = as.integer(V3),
      state = as.character(V4),
      seglen_bp = as.integer(end - start + 1L)
    ) %>%
    mutate(chr = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr))) %>%
    arrange(chr, start, end)
}

## ---------------------------
## State -> dosage factories
## ---------------------------
state_to_dosage_hard <- function(state_chr) {
  s <- tolower(state_chr)
  out <- rep(NA_real_, length(s))
  out[grepl("^b73$", s)] <- 0
  out[grepl("het|hetero|b73/intro", s)] <- 1
  out[grepl("introgression", s)] <- 2
  out
}

# soft mapping: state -> P(teo) -> dosage = 2*P(teo)
# (edit these probabilities if your caller is biased)
state_to_dosage_soft <- function(state_chr,
                                 p_teo_b73 = 0.0,
                                 p_teo_het = 0.5,
                                 p_teo_intro = 1.0) {
  s <- tolower(state_chr)
  p <- rep(NA_real_, length(s))
  p[grepl("^b73$", s)] <- p_teo_b73
  p[grepl("het|hetero|b73/intro", s)] <- p_teo_het
  p[grepl("introgression", s)] <- p_teo_intro
  2 * p
}


## ---------------------------
## Call ancestry dosage at positions for one chr (fast)
## Optionally return local seglen weights
## ---------------------------
call_chr <- function(tb_chr, pos_vec,
                     dosage_mode = c("hard","soft_state"),
                     soft_probs = list(p_teo_b73 = 0, p_teo_het = 0.5, p_teo_intro = 1),
                     weight_by_seglen = FALSE,
                     weight_transform = c("none","log1p","sqrt"),
                     weight_center = c("median","mean","none")) {
  
  dosage_mode <- match.arg(dosage_mode)
  weight_transform <- match.arg(weight_transform)
  weight_center <- match.arg(weight_center)
  
  starts <- tb_chr$start
  ends   <- tb_chr$end
  states <- tb_chr$state
  seglen <- tb_chr$seglen_bp
  
  idx <- findInterval(pos_vec, starts)
  idx[idx == 0] <- 1L
  
  bad <- pos_vec > ends[idx]
  st  <- states[idx]
  st[bad] <- NA_character_
  
  if (dosage_mode == "hard") {
    dos <- state_to_dosage_hard(st)
  } else {
    dos <- state_to_dosage_soft(
      st,
      p_teo_b73   = soft_probs$p_teo_b73,
      p_teo_het   = soft_probs$p_teo_het,
      p_teo_intro = soft_probs$p_teo_intro
    )
  }
  
  if (!weight_by_seglen) return(dos)
  
  w <- as.numeric(seglen[idx])
  w[bad] <- NA_real_
  
  if (weight_transform == "log1p") w <- log1p(w)
  if (weight_transform == "sqrt")  w <- sqrt(w)
  
  if (weight_center == "median") {
    m <- median(w, na.rm = TRUE); if (is.finite(m) && m > 0) w <- w / m
  } else if (weight_center == "mean") {
    m <- mean(w, na.rm = TRUE); if (is.finite(m) && m > 0) w <- w / m
  }
  
  list(dosage = dos, weight = w)
}

## ---------------------------
## Build marker grid and map
## ---------------------------
make_markers <- function(x_list, step_bp = 100000L) {
  
  chr_len <- bind_rows(x_list, .id = "Taxa") %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(as.integer(stringr::str_remove(chr, "^chr")))
  
  purrr::pmap_dfr(
    list(chr_len$chr, chr_len$chr_len, list(step_bp)),
    function(chr, L, step_bp) {
      pos <- seq.int(from = step_bp, to = L, by = step_bp)
      tibble::tibble(
        chr    = as.character(chr),                         # "chr3"
        Chrom  = as.integer(stringr::str_remove(chr, "^chr")),# 3
        pos    = as.integer(pos),
        Marker = paste0(chr, ":", pos)                      # "chr3:100000"
      )
    }
  )
}
## ---------------------------
## Helper: parse marker if needed
## Works for:
##   "chr3:100000"
##   "chr3:118700001-118725000"  (midpoint)
## ---------------------------
parse_marker <- function(m) {
  m <- as.character(m)
  chr <- sub(":.*$", "", m)
  
  rest <- sub("^.*?:", "", m)
  if (grepl("-", rest)) {
    a <- as.integer(sub("-.*$", "", rest))
    b <- as.integer(sub("^.*?-", "", rest))
    pos <- as.integer(floor((a + b) / 2))
    start <- a; end <- b
  } else {
    pos <- as.integer(rest)
    start <- pos; end <- pos
  }
  
  Chrom <- suppressWarnings(as.integer(sub("^chr", "", chr)))
  list(chr = chr, Chrom = Chrom, pos = pos, start = start, end = end)
}


## ---------------------------
## Global introgression covariates (per individual)
## ---------------------------
calc_intro_covariates <- function(x_std, intro_regex = "introgression") {
  ids <- names(x_std)
  out <- lapply(ids, function(id) {
    tb <- x_std[[id]]
    s  <- tolower(tb$state)
    is_intro <- grepl(intro_regex, s)
    intro_bp <- sum(tb$seglen_bp[is_intro], na.rm = TRUE)
    n_intro  <- sum(is_intro, na.rm = TRUE)
    mean_intro_kb <- if (n_intro > 0) mean(tb$seglen_bp[is_intro], na.rm = TRUE) / 1000 else 0
    tibble(
      Taxa = id,
      intro_bp = intro_bp,
      n_intro_segments = n_intro,
      mean_intro_seg_kb = mean_intro_kb
    )
  }) %>% bind_rows()
  
  # (optional) normalize by an approximate genome length = max end across all tracts
  genome_bp <- bind_rows(x_std) %>% summarise(L = max(end, na.rm = TRUE)) %>% pull(L)
  if (is.finite(genome_bp) && genome_bp > 0) {
    out <- out %>% mutate(intro_frac = intro_bp / genome_bp)
  } else {
    out <- out %>% mutate(intro_frac = NA_real_)
  }
  out
}

## ---------------------------
## Standardize rrBLUP GWAS output (robust)
## rrBLUP puts -log10(p) under the TRAIT column name
## ---------------------------
standardize_rrblup_gwas <- function(g_raw, fam_code, trait_col) {
  
  # rrBLUP sometimes returns list-of-dataframes
  if (is.list(g_raw) && !is.data.frame(g_raw)) {
    if (trait_col %in% names(g_raw)) g_raw <- g_raw[[trait_col]] else g_raw <- g_raw[[1]]
  }
  
  g_tbl <- tibble::as_tibble(g_raw)
  
  # marker col
  mcol <- intersect(c("Marker","marker","SNP","snp","ID","id"), names(g_tbl))[1]
  if (is.na(mcol)) stop("No marker column found in rrBLUP output. Columns: ",
                        paste(names(g_tbl), collapse = ", "))
  
  # score col (-log10p)
  score_col <- intersect(names(g_tbl), c(trait_col, "trait", "y"))[1]
  if (is.na(score_col)) {
    known <- c(mcol, "Chrom","chrom","Chr","chr","Pos","pos","Position","position",
               "Effect","effect","P.value","p.value","R2","r2","MAF","maf")
    cand <- setdiff(names(g_tbl), known)
    cand <- cand[sapply(g_tbl[cand], is.numeric)]
    score_col <- cand[1]
  }
  if (is.na(score_col)) stop("Could not identify -log10(p) column in rrBLUP output.")
  
  out <- g_tbl %>%
    dplyr::rename(Marker = !!rlang::sym(mcol)) %>%
    dplyr::mutate(
      family  = fam_code,
      mlog10p = suppressWarnings(as.numeric(.data[[score_col]])),
      p       = 10^(-mlog10p)
    )
  
  # If rrBLUP gave Chrom/Pos, use them; otherwise parse from Marker
  chrom_col <- intersect(names(out), c("Chrom","chrom","Chr","chr"))[1]
  pos_col   <- intersect(names(out), c("Pos","pos","Position","position"))[1]
  
  if (!is.na(chrom_col) && !is.na(pos_col)) {
    chrom_raw <- as.character(out[[chrom_col]])
    # handle both "3" and "chr3"
    Chrom <- suppressWarnings(as.integer(sub("^chr", "", chrom_raw)))
    pos   <- suppressWarnings(as.integer(out[[pos_col]]))
    chr   <- ifelse(grepl("^chr", chrom_raw), chrom_raw, paste0("chr", Chrom))
    
    out$Chrom <- Chrom
    out$pos   <- pos
    out$chr   <- chr
  } else {
    pm <- lapply(out$Marker, parse_marker)
    out$chr   <- vapply(pm, `[[`, character(1), "chr")
    out$Chrom <- vapply(pm, `[[`, integer(1), "Chrom")
    out$pos   <- vapply(pm, `[[`, integer(1), "pos")
  }
  
  out %>%
    dplyr::select(Marker, chr, Chrom, pos, family, mlog10p, p, dplyr::everything()) %>%
    dplyr::arrange(p)
}


## ---------------------------
## Main runner: family-wise introgression GWAS with options
## ---------------------------
run_introgression_gwas_flexible <- function(
    x_rds,
    pheno_file,
    geno_dir,
    trait_col,
    pheno_id_col = 1,
    pheno_trait_col = 3,
    families_keep = c("Zd","Zl","Zv","Zx"),
    step_bp = 100000L,
    
    ancestry_mode = c("hard","soft_state","posterior_mat"),
    # for soft_state:
    soft_probs = list(p_teo_b73 = 0.0, p_teo_het = 0.5, p_teo_intro = 1.0),
    
    # for posterior_mat:
    posterior_mat = NULL,          # matrix of P(teo) OR dosage, rows=Taxa, cols=Marker (chr:pos)
    posterior_is_prob = TRUE,      # TRUE => dosage = 2*P(teo); FALSE => values already dosage 0..2
    
    add_seg_covariates = TRUE,     # add intro_frac, intro_bp, n_intro_segments, mean_intro_seg_kb
    weight_by_seglen = FALSE,      # local weight by segment length at each marker
    weight_transform = c("none","log1p","sqrt"),
    weight_center = c("median","mean","none"),
    
    npc = 3,                       # PCs from K_snp.rel
    min_maf = 0.01,
    max_miss = 0.2,
    
    do_meta_acat = TRUE,
    out_dir = file.path(geno_dir, "Introgression_GWAS_familywise_flexible")
) {
  ancestry_mode <- match.arg(ancestry_mode)
  weight_transform <- match.arg(weight_transform)
  weight_center <- match.arg(weight_center)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ---- read phenotype ----
  pheno_raw <- read.csv(pheno_file, stringsAsFactors = FALSE)
  
  pheno <- pheno_raw %>%
    dplyr::select(!!pheno_id_col, !!pheno_trait_col) %>%
    setNames(c("Taxa", trait_col)) %>%
    mutate(
      Taxa = stripB(as.character(Taxa)),
      !!trait_col := suppressWarnings(as.numeric(.data[[trait_col]]))
    ) %>%
    filter(!is.na(.data[[trait_col]])) %>%
    mutate(family = get_family(Taxa))
  
  ## ---- read kinship ----
  kid_file <- file.path(geno_dir, "K_snp.rel.id")
  K_file   <- file.path(geno_dir, "K_snp.rel")
  
  kid <- fread(kid_file, header = FALSE)
  if (kid[1, V1] == "#FID" || kid[1, V2] == "IID") kid <- kid[-1]
  taxa_all <- stripB(as.character(kid$V2))
  
  K <- as.matrix(fread(K_file, header = FALSE))
  stopifnot(nrow(K) == length(taxa_all), ncol(K) == length(taxa_all))
  
  rownames(K) <- taxa_all
  colnames(K) <- taxa_all
  
  ## ---- read + standardize tract list ----
  x <- readRDS(x_rds)
  names(x) <- stripB(names(x))
  x_std <- lapply(x, std_tbl)
  
  ## ---- restrict to common taxa ----
  taxa_ok <- Reduce(intersect, list(pheno$Taxa, rownames(K), names(x_std)))
  pheno   <- pheno %>% filter(Taxa %in% taxa_ok)
  x_std   <- x_std[pheno$Taxa]
  K       <- K[taxa_ok, taxa_ok, drop = FALSE]
  
  ## ---- marker grid ----
  markers_map <- make_markers(x_std, step_bp = step_bp)
  markers_by_chr <- split(markers_map, markers_map$chr)
  
  ## ---- build ancestry dosage matrix (unless posterior_mat supplied) ----
  if (ancestry_mode != "posterior_mat") {
    
    build_sample <- function(tb) {
      v <- rep(NA_real_, nrow(markers_map))
      if (weight_by_seglen) wv <- rep(NA_real_, nrow(markers_map))
      
      for (cc in names(markers_by_chr)) {
        pos_vec <- markers_by_chr[[cc]]$pos
        tb_cc   <- tb %>% filter(chr == cc)
        if (nrow(tb_cc) == 0) next
        
        ans <- call_chr(
          tb_cc, pos_vec,
          dosage_mode = ancestry_mode,
          soft_probs = soft_probs,
          weight_by_seglen = weight_by_seglen,
          weight_transform = weight_transform,
          weight_center = weight_center
        )
        
        idx_global <- match(markers_by_chr[[cc]]$Marker, markers_map$Marker)
        
        if (!weight_by_seglen) {
          v[idx_global] <- ans
        } else {
          v[idx_global]  <- ans$dosage
          wv[idx_global] <- ans$weight
        }
      }
      
      if (!weight_by_seglen) return(v)
      
      # length-weighted dosage (optional): dosage * weight
      v_w <- v
      ok <- !is.na(v) & !is.na(wv)
      v_w[ok] <- v[ok] * wv[ok]
      v_w
    }
    
    message("Building ancestry matrix from tract states...")
    G_dos <- do.call(rbind, lapply(x_std, build_sample))
    rownames(G_dos) <- names(x_std)
    colnames(G_dos) <- markers_map$Marker
    
  } else {
    if (is.null(posterior_mat)) {
      stop("ancestry_mode='posterior_mat' requires posterior_mat (rows=Taxa, cols=Marker).")
    }
    # normalize IDs
    rownames(posterior_mat) <- stripB(rownames(posterior_mat))
    colnames(posterior_mat) <- stripB(colnames(posterior_mat))
    
    common_ids <- intersect(pheno$Taxa, rownames(posterior_mat))
    common_mk  <- intersect(markers_map$Marker, colnames(posterior_mat))
    if (length(common_ids) < 30) stop("Too few samples in posterior_mat after matching.")
    if (length(common_mk) < 100) stop("Too few markers in posterior_mat after matching.")
    
    G_dos <- posterior_mat[common_ids, common_mk, drop = FALSE]
    if (posterior_is_prob) G_dos <- 2 * G_dos  # P(teo)->dosage 0..2
    
    # expand to full marker set (optional) by keeping only common_mk
    pheno <- pheno %>% filter(Taxa %in% common_ids)
    K <- K[common_ids, common_ids, drop = FALSE]
    markers_map <- markers_map %>% filter(Marker %in% common_mk)
  }
  
  ## ---- optional global covariates from tracts ----
  intro_cov <- NULL
  if (add_seg_covariates) {
    intro_cov <- calc_intro_covariates(x_std)
    pheno <- pheno %>%
      left_join(intro_cov, by = c("Taxa" = "Taxa"))
  }
  
  ## ---------------------------
  ## Family-wise GWAS
  ## ---------------------------
  run_one_family <- function(fam_code) {
    
    taxa_f <- pheno %>%
      filter(family == fam_code) %>%
      pull(Taxa) %>%
      intersect(rownames(G_dos)) %>%
      intersect(rownames(K))
    
    if (length(taxa_f) < 30) {
      message("Skipping ", fam_code, " (n<30): ", length(taxa_f))
      return(NULL)
    }
    
    # subset phenotype
    Y_f <- pheno %>%
      filter(Taxa %in% taxa_f) %>%
      arrange(match(Taxa, taxa_f))
    
    # kinship subset
    K_f <- K[taxa_f, taxa_f, drop = FALSE]
    
    # PCs from K
    pcs <- as.data.frame(calc_pcs_from_K(K_f, npc = npc))
    
    # pheno df for rrBLUP
    ph_rr <- data.frame(
      gid = taxa_f,
      y   = suppressWarnings(as.numeric(Y_f[[trait_col]])),
      pcs,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    fixed_cols <- colnames(pcs)
    
    if (add_seg_covariates) {
      # choose which covariates you actually want
      ph_rr$intro_frac <- suppressWarnings(as.numeric(Y_f$intro_frac))
      ph_rr$intro_bp   <- suppressWarnings(as.numeric(Y_f$intro_bp))
      ph_rr$n_intro_segments <- suppressWarnings(as.numeric(Y_f$n_intro_segments))
      ph_rr$mean_intro_seg_kb <- suppressWarnings(as.numeric(Y_f$mean_intro_seg_kb))
      
      # include a conservative subset as fixed effects (avoid overfitting)
      fixed_cols <- c(fixed_cols, "intro_frac")
    }
    
    # ancestry dosage matrix for this family
    M <- G_dos[taxa_f, , drop = FALSE]
    
    # missingness filter
    miss  <- colMeans(is.na(M))
    keep1 <- miss <= max_miss
    
    # MAF filter on dosage 0..2
    maf <- apply(M[, keep1, drop = FALSE], 2, maf_from_dosage)
    keep2 <- maf >= min_maf
    
    M2 <- M[, keep1, drop = FALSE]
    M2 <- M2[, keep2, drop = FALSE]
    
    if (ncol(M2) < 100) {
      message("Skipping ", fam_code, " (too few markers after filters): ", ncol(M2))
      return(NULL)
    }
    
    message("Running ", fam_code, ": n=", length(taxa_f), " markers=", ncol(M2))
    
    # rrBLUP marker score in [-1,+1]:
    # dosage in [0,2] -> score = dosage - 1
    Score <- M2 - 1
    
    # markers x lines
    Gmat <- t(Score)
    colnames(Gmat) <- taxa_f
    
    # impute NA by marker mean
    if (anyNA(Gmat)) {
      mu <- rowMeans(Gmat, na.rm = TRUE)
      ij <- which(is.na(Gmat), arr.ind = TRUE)
      Gmat[ij] <- mu[ij[, 1]]
    }
    
    # marker map subset + order
    map <- markers_map %>%
      dplyr::filter(Marker %in% rownames(Gmat)) %>%
      dplyr::arrange(match(Marker, rownames(Gmat)))
    
    stopifnot(identical(map$Marker, rownames(Gmat)))
    
    geno_df <- data.frame(
      Marker = map$Marker,
      Chrom  = map$Chrom,
      Pos    = map$pos,
      as.data.frame(Gmat, check.names = FALSE),
      check.names = FALSE
    )
    
    # alignment checks
    stopifnot(identical(colnames(geno_df)[4:ncol(geno_df)], ph_rr$gid))
    stopifnot(identical(rownames(K_f), ph_rr$gid))
    
    # GWAS
    g_raw <- rrBLUP::GWAS(
      pheno   = ph_rr,
      geno    = geno_df,
      fixed   = fixed_cols,
      K       = K_f,
      n.PC    = 0,
      min.MAF = 0,
      P3D     = TRUE,
      plot    = FALSE
    )
    
    out <- standardize_rrblup_gwas(g_raw, fam_code = fam_code, trait_col = trait_col)
  
    
    # write per-family
    write.csv(out,
              file.path(out_dir, paste0("introgression_GWAS_", fam_code, "_step", step_bp, "_", ancestry_mode, ".csv")),
              row.names = FALSE)
    
    out
  }
  
  res_list <- lapply(families_keep, run_one_family)
  names(res_list) <- families_keep
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  
  res_all <- bind_rows(res_list)
  
  write.csv(res_all,
            file.path(out_dir, paste0("introgression_GWAS_allFamilies_step", step_bp, "_", ancestry_mode, ".csv")),
            row.names = FALSE)
  
  meta <- NULL
  if (do_meta_acat && nrow(res_all) > 0) {
    meta <- res_all %>%
      dplyr::group_by(Marker, chr, Chrom, pos) %>%
      dplyr::summarise(
        n_fam  = dplyr::n_distinct(family),
        min_p  = min(p, na.rm = TRUE),
        acat_p = acat(p),
        .groups = "drop"
      ) %>%
      dplyr::arrange(acat_p)
    
    write.csv(meta,
              file.path(out_dir, paste0("introgression_GWAS_meta_ACAT_step", step_bp, "_", ancestry_mode, ".csv")),
              row.names = FALSE)
  }
  
  list(
    res_list = res_list,
    res_all  = res_all,
    meta     = meta,
    G_dosage = G_dos,
    markers  = markers_map,
    pheno    = pheno,
    K_snp    = K
  )
}

## ============================================================
## EXAMPLES
## ============================================================

## 1) Your current workflow (hard 0/1/2 -> score -1/0/+1)
ans <- run_introgression_gwas_flexible(
  x_rds      = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BZea_Introgression_Finder/results_list_new_name.rds",
  pheno_file = "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv",
  geno_dir   = "/Users/nirwantandukar/Documents/Research/data/BZea/genotype",
  trait_col  = "DTS",               # or "DTS" or "trait"
  pheno_id_col = 1,
  pheno_trait_col = 2,
  
  # Windows.
  step_bp    = 5000000L,
  
  #ancestry_mode = "hard", # hard or soft_state or posterior_mat
  #ancestry_mode = "soft_state", # soft requireds soft probabilities as well
  
  #soft_state
  ancestry_mode = "soft_state",
  soft_probs = list(p_teo_b73 = 0.02, p_teo_het = 0.5, p_teo_intro = 0.98),
  
  
  # posterior_mat
  #ancestry_mode = "posterior_mat",
  #posterior_mat = post,
  #posterior_is_prob = TRUE,
  
  add_seg_covariates = FALSE,
  weight_by_seglen = FALSE,
  npc = 3,
  min_maf = 0.005,
  max_miss = 0.5,
  do_meta_acat = TRUE
)


str(ans)
## 2) Soft mapping from states (probabilities -> dosage = 2*P(teo))
# ans <- run_introgression_gwas_flexible(
#   ...,
#   ancestry_mode = "soft_state",
#   soft_probs = list(p_teo_b73 = 0.02, p_teo_het = 0.5, p_teo_intro = 0.98)
# )

## 3) HMM posterior matrix (P(teo) per marker) -> dosage=2*P(teo)
##    You build this elsewhere (could be from your HMM/VCF pipeline),
##    then pass it in as posterior_mat.
# post <- readRDS("/path/to/posterior_Pteo_matrix_step10k.rds")  # rows=Taxa, cols=chr:pos
# ans <- run_introgression_gwas_flexible(
#   ...,
#   ancestry_mode = "posterior_mat",
#   posterior_mat = post,
#   posterior_is_prob = TRUE
# )

## 4) Weight marker predictor by segment length (local weighting)
# ans <- run_introgression_gwas_flexible(
#   ...,
#   ancestry_mode = "hard",
#   weight_by_seglen = TRUE,
#   weight_transform = "log1p",
#   weight_center = "median"
# )




## ============================================================
## Example: annotate your already-generated res_all
## ============================================================
gff3 <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

ann <- annotate_introgression_bins(
  res_all = ans$res_all,      # OR your res_all tibble in memory
  gff3_file = gff3,
  window_bp = 25000L,
  sig_mlog10p = 3,
  out_prefix = "/Users/nirwantandukar/Documents/Research/data/BZea/genotype/Introgression_GWAS_familywise_flexible/DTA_soft_100k"
)

ann$ann_hits %>% head()
ann$gene_summary %>% head()


