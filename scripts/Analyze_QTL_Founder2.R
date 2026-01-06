## ============================================================
## BZea GridLMM QTL scan (PACKAGE-READY CORE)
##
## One main function:
##   run_bzea_gridlmm_scan(scan_mode = "single" or "window", step_bp = 100000L, ...)
##
## scan_mode:
##   - "single": test EVERY marker column in geno_mat (chr:pos)
##   - "window": first aggregate markers into fixed genomic windows (step_bp),
##               then test ONE value per window (mean/median across markers in window)
## analysis_mode:
##   - "combined" | "familywise" | "both"
##   - include_family_fixed = TRUE/FALSE   (ONLY affects COMBINED model)
## Output:
##   $res_all       : combined scan results (all families pooled)
##   $res_by_family : list of family-wise scan results
##   $marker_map    : marker/bin map used in the scan (chr, pos, bin_start/bin_end if windowed)
##
## ============================================================


suppressPackageStartupMessages({
  library(GridLMM)
  library(data.table)
  library(dplyr)
  library(stringr)
})

stripB <- function(x) sub("\\.B$", "", as.character(x))

fix_id <- function(x) {
  x <- as.character(x)
  x <- gsub("-", ".", x)
  stripB(x)
}

stopifnot1 <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

impute_colmean <- function(X) {
  X <- as.matrix(X)
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (anyNA(v)) {
      mu <- mean(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      X[, j] <- v
    }
  }
  X
}

safe_p <- function(p, min_p = 1e-300) {
  p <- suppressWarnings(as.numeric(p))
  p[is.na(p)] <- 1
  p <- pmin(pmax(p, min_p), 1)
  p
}

parse_marker_map <- function(markers) {
  markers <- as.character(markers)
  chr <- sub(":.*$", "", markers)
  pos <- suppressWarnings(as.integer(sub("^.*?:", "", markers)))
  tibble(Marker = markers, chr = chr, pos = pos)
}

make_windowed_X <- function(geno_mat,
                            marker_map,
                            step_bp = 100000L,
                            agg = c("mean", "median"),
                            min_markers = 1L) {
  
  agg <- match.arg(agg)
  step_bp <- as.integer(step_bp)
  stopifnot1(step_bp >= 1L, "step_bp must be >= 1")
  
  geno_mat <- as.matrix(geno_mat)
  
  # keep as tibble so dplyr sees columns cleanly
  marker_map <- dplyr::as_tibble(marker_map)
  
  stopifnot1(ncol(geno_mat) == nrow(marker_map),
             "make_windowed_X: ncol(geno_mat) must equal nrow(marker_map)")
  
  marker_map <- marker_map %>%
    dplyr::mutate(
      chr = as.character(.data$chr),
      chr = ifelse(stringr::str_detect(.data$chr, "^chr"), .data$chr, paste0("chr", .data$chr)),
      pos = as.integer(.data$pos)
    )
  
  bin_id <- (pmax(marker_map$pos, 1L) - 1L) %/% step_bp
  marker_map <- marker_map %>%
    dplyr::mutate(
      bin_start = bin_id * step_bp + 1L,
      bin_end   = .data$bin_start + step_bp - 1L,
      Bin       = paste0(.data$chr, ":", as.integer((.data$bin_start + .data$bin_end) / 2))
    )
  
  idx_by_bin <- split(seq_len(nrow(marker_map)), marker_map$Bin)
  
  if (!is.null(min_markers) && as.integer(min_markers) > 1L) {
    keep_bins <- names(idx_by_bin)[vapply(idx_by_bin, length, integer(1)) >= as.integer(min_markers)]
    idx_by_bin <- idx_by_bin[keep_bins]
  }
  
  # columns = bins, rows = lines
  Xbin <- vapply(idx_by_bin, function(idx) {
    Xsub <- geno_mat[, idx, drop = FALSE]
    if (ncol(Xsub) == 1) return(as.numeric(Xsub[, 1]))
    if (agg == "mean") return(rowMeans(Xsub, na.rm = TRUE))
    apply(Xsub, 1, median, na.rm = TRUE)
  }, FUN.VALUE = numeric(nrow(geno_mat)))
  
  Xbin <- as.matrix(Xbin)
  rownames(Xbin) <- rownames(geno_mat)
  
  mm2 <- marker_map %>%
    dplyr::group_by(Bin) %>%                      # <- no .data
    dplyr::summarise(
      chr = dplyr::first(chr),
      bin_start = min(bin_start, na.rm = TRUE),
      bin_end   = max(bin_end,   na.rm = TRUE),
      pos       = as.integer((dplyr::first(bin_start) + dplyr::first(bin_end)) / 2),
      n_markers = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      chr_num = suppressWarnings(as.integer(stringr::str_extract(chr, "\\d+")))
    ) %>%
    dplyr::arrange(chr_num, pos) %>%
    dplyr::select(-chr_num) %>%
    dplyr::mutate(Marker = Bin) %>%
    dplyr::select(Marker, chr, pos, bin_start, bin_end, n_markers)
  
  colnames(Xbin) <- mm2$Marker
  
  list(X = Xbin, map = mm2)
}

read_plink2_rel <- function(prefix) {
  f_rel <- paste0(prefix, ".rel")
  f_id  <- paste0(prefix, ".rel.id")
  stopifnot1(file.exists(f_rel), paste0("Missing: ", f_rel))
  stopifnot1(file.exists(f_id),  paste0("Missing: ", f_id))
  
  ids <- data.table::fread(f_id, data.table = FALSE)
  if ("IID" %in% names(ids)) {
    id <- ids[["IID"]]
  } else {
    id <- ids[[ncol(ids)]]
  }
  id <- fix_id(id)
  
  Kdt <- data.table::fread(f_rel, header = FALSE, data.table = FALSE)
  K <- as.matrix(Kdt)
  stopifnot1(nrow(K) == ncol(K), "K is not square — .rel read failed?")
  stopifnot1(nrow(K) == length(id), "ID count != K dimension — wrong .rel/.rel.id pairing?")
  
  rownames(K) <- id
  colnames(K) <- id
  K
}

read_K_matrix_txt <- function(f) {
  stopifnot1(file.exists(f), paste0("Missing: ", f))
  Kdf <- data.table::fread(f, data.table = FALSE)
  ids <- fix_id(Kdf[[1]])
  K <- as.matrix(Kdf[, -1, drop = FALSE])
  storage.mode(K) <- "double"
  rownames(K) <- ids
  if (nrow(K) == ncol(K)) colnames(K) <- ids
  K
}

read_loco_K <- function(K_dir,
                        chr,
                        K_format = c("plink2_rel", "text"),
                        K_prefix_pattern = "LOCO_chr{chr}",
                        K_text_pattern   = "K_matrix_chr{chr}.txt") {
  K_format <- match.arg(K_format)
  chr <- as.character(chr)
  
  if (K_format == "plink2_rel") {
    prefix <- file.path(K_dir, gsub("\\{chr\\}", chr, K_prefix_pattern))
    return(read_plink2_rel(prefix))
  } else {
    f <- file.path(K_dir, gsub("\\{chr\\}", chr, K_text_pattern))
    return(read_K_matrix_txt(f))
  }
}

standardize_gridlmm_out <- function(out_df) {
  out_df <- as.data.frame(out_df)
  stopifnot1("p_value_ML" %in% names(out_df), "No p_value_ML column found in GridLMM results.")
  out_df$p <- safe_p(out_df$p_value_ML)
  out_df$mlog10p <- -log10(out_df$p)
  out_df
}

## ============================================================
## MAIN (PACKAGE-READY)
## ============================================================
run_bzea_gridlmm_scan <- function(
    pheno_file,
    trait_col,
    id_col = 1,
    
    geno_file,                     # first col IDs, remaining cols markers "chr:pos"
    scan_mode = c("single", "window"),
    step_bp = 100000L,             # only if window
    window_agg = c("mean", "median"),
    min_markers_per_window = 1L,
    
    K_dir,
    K_format = c("plink2_rel", "text"),
    K_prefix_pattern = "LOCO_chr{chr}",
    K_text_pattern   = "K_matrix_chr{chr}.txt",
    
    analysis_mode = c("both", "combined", "familywise"),
    include_family_fixed = TRUE,   # only for COMBINED
    family_col = NULL,             # optional explicit family column in pheno
    covar_cols = character(),
    
    method = "ML",
    cores = 1,
    
    chr_set = paste0("chr", 1:10),
    min_n_family = 20,
    verbose = TRUE
) {
  scan_mode <- match.arg(scan_mode)
  window_agg <- match.arg(window_agg)
  K_format <- match.arg(K_format)
  analysis_mode <- match.arg(analysis_mode)
  
  do_combined  <- analysis_mode %in% c("both", "combined")
  do_familywise <- analysis_mode %in% c("both", "familywise")
  
  ## ---- phenotype
  ph <- data.table::fread(pheno_file, data.table = FALSE)
  stopifnot1(id_col >= 1 && id_col <= ncol(ph), "id_col index out of range for phenotype table.")
  idname <- names(ph)[id_col]
  stopifnot1(trait_col %in% names(ph), paste0("trait_col not found: ", trait_col))
  
  keep_cols <- unique(c(idname, trait_col, covar_cols, family_col))
  keep_cols <- keep_cols[!is.na(keep_cols)]
  ph <- ph[, intersect(keep_cols, names(ph)), drop = FALSE]
  
  dat <- data.frame(
    Line = fix_id(ph[[idname]]),
    y    = suppressWarnings(as.numeric(ph[[trait_col]])),
    stringsAsFactors = FALSE
  )
  
  # always compute Family labels if familywise OR combined-with-family-fixed
  need_family_labels <- do_familywise || (do_combined && include_family_fixed)
  if (need_family_labels) {
    if (!is.null(family_col) && family_col %in% names(ph)) {
      dat$Family <- as.factor(as.character(ph[[family_col]]))
    } else {
      dat$Family <- as.factor(sub("\\..*$", "", dat$Line))
    }
  }
  
  # add covariates
  if (length(covar_cols) > 0) {
    for (cc in covar_cols) {
      stopifnot1(cc %in% names(ph), paste0("covar_cols includes missing column: ", cc))
      dat[[cc]] <- suppressWarnings(as.numeric(ph[[cc]]))
    }
  }
  
  dat <- dat[!is.na(dat$Line) & !is.na(dat$y), , drop = FALSE]
  dat$y <- dat$y - mean(dat$y, na.rm = TRUE)
  
  ## ---- genotype
  gdf <- data.table::fread(geno_file, data.table = FALSE)
  stopifnot1(ncol(gdf) >= 2, "geno_file must be: ID column + >=1 marker columns")
  
  ids_g <- fix_id(gdf[[1]])
  geno_mat <- as.matrix(gdf[, -1, drop = FALSE])
  storage.mode(geno_mat) <- "double"
  rownames(geno_mat) <- ids_g
  
  marker_map <- parse_marker_map(colnames(geno_mat))
  
  ## ---- choose X (single vs window)  <<<< THIS is the selector
  if (scan_mode == "single") {
    X_all <- geno_mat
    map_use <- marker_map
    map_use$bin_start <- map_use$pos
    map_use$bin_end   <- map_use$pos
    map_use$n_markers <- 1L
  } else {
    win <- make_windowed_X(
      geno_mat = geno_mat,
      marker_map = marker_map,
      step_bp = as.integer(step_bp),
      agg = window_agg,
      min_markers = as.integer(min_markers_per_window)
    )
    X_all <- win$X
    map_use <- win$map
  }
  
  ## ---- keep chr_set
  map_use$chr <- as.character(map_use$chr)
  map_use <- map_use[map_use$chr %in% chr_set, , drop = FALSE]
  X_all <- X_all[, map_use$Marker, drop = FALSE]
  
  ## ---- intersect IDs (pheno/genotype)
  keep_ids <- intersect(dat$Line, rownames(X_all))
  stopifnot1(length(keep_ids) >= 10, "Too few overlapping IDs between phenotype and genotype.")
  dat <- dat[dat$Line %in% keep_ids, , drop = FALSE]
  X_all <- X_all[keep_ids, , drop = FALSE]
  
  ## ---- chromosomes present
  chr_levels <- unique(map_use$chr)
  chr_levels <- chr_levels[chr_levels %in% chr_set]
  
  if (verbose) {
    message(
      "analysis_mode=", analysis_mode,
      " | scan_mode=", scan_mode,
      if (scan_mode == "window") paste0(" (step_bp=", step_bp, ", agg=", window_agg, ")") else "",
      "\nN lines=", nrow(dat),
      " | N tests=", ncol(X_all),
      "\nchr: ", paste(chr_levels, collapse = ", ")
    )
  }
  
  ## ---- helper: run one scan set (combined or one family)
  run_one_scanset <- function(dat_sub, include_family_fixed_here, fam_label = NA_character_) {
    out_by_chr <- list()
    
    for (chr_key in chr_levels) {
      chr_num <- suppressWarnings(as.integer(str_extract(chr_key, "\\d+")))
      if (!is.finite(chr_num)) next
      
      mk <- map_use$Marker[map_use$chr == chr_key]
      if (length(mk) < 2) next
      
      # LOCO K for this chr
      K <- read_loco_K(
        K_dir = K_dir,
        chr = chr_num,
        K_format = K_format,
        K_prefix_pattern = K_prefix_pattern,
        K_text_pattern = K_text_pattern
      )
      
      rownames(K) <- fix_id(rownames(K))
      colnames(K) <- fix_id(colnames(K))
      
      keep <- intersect(dat_sub$Line, rownames(K))
      if (length(keep) < min_n_family) next
      
      # order everything
      dat_chr <- dat_sub[match(keep, dat_sub$Line), , drop = FALSE]
      keep <- dat_chr$Line
      
      K_use <- K[keep, keep, drop = FALSE]
      X_use <- X_all[keep, mk, drop = FALSE]
      X_use <- impute_colmean(X_use)
      
      # build ydat
      ydat <- data.frame(Line = keep, y = dat_chr$y, stringsAsFactors = FALSE)
      if (include_family_fixed_here) ydat$Family <- dat_chr$Family
      if (length(covar_cols) > 0) {
        for (cc in covar_cols) ydat[[cc]] <- dat_chr[[cc]]
      }
      
      fixed_terms <- c("1")
      if (include_family_fixed_here) fixed_terms <- c(fixed_terms, "Family")
      if (length(covar_cols) > 0) fixed_terms <- c(fixed_terms, covar_cols)
      fml <- as.formula(paste0("y ~ ", paste(fixed_terms, collapse = " + "), " + (1|Line)"))
      
      res <- GridLMM_GWAS(
        formula         = fml,
        test_formula    = ~ 1,
        reduced_formula = ~ 0,
        data            = ydat,
        X               = X_use,
        X_ID            = "Line",
        relmat          = list(Line = K_use),
        method          = method,
        centerX         = FALSE,
        scaleX          = FALSE,
        fillNAX         = FALSE,
        verbose         = FALSE,
        mc.cores        = cores
      )
      
      out <- as.data.frame(res$results)
      out$Marker <- colnames(X_use)
      out$chr <- chr_key
      if (!is.na(fam_label)) out$Family <- fam_label
      
      # attach positions / bin bounds
      mm <- map_use %>%
        filter(Marker %in% out$Marker) %>%
        distinct(Marker, chr, pos, bin_start, bin_end, n_markers)
      
      out <- out %>%
        left_join(mm, by = "Marker") %>%
        standardize_gridlmm_out()
      
      out_by_chr[[chr_key]] <- out
    }
    
    dplyr::bind_rows(out_by_chr) %>% arrange(p)
  }
  
  ## ---- COMBINED
  res_all <- NULL
  if (do_combined) {
    res_all <- run_one_scanset(dat_sub = dat, include_family_fixed_here = need_family_labels && include_family_fixed)
  }
  
  ## ---- FAMILY-WISE (separate scans per family)
  res_by_family <- list()
  if (do_familywise) {
    stopifnot1("Family" %in% names(dat), "Family-wise requested but Family labels are missing.")
    fams <- levels(dat$Family)
    
    for (fam in fams) {
      dat_f <- dat[dat$Family == fam, , drop = FALSE]
      if (nrow(dat_f) < min_n_family) next
      
      # within-family: DO NOT include Family as fixed effect
      fam_res <- run_one_scanset(dat_sub = dat_f, include_family_fixed_here = FALSE, fam_label = fam)
      if (nrow(fam_res) > 0) res_by_family[[fam]] <- fam_res
    }
  }
  
  list(
    analysis_mode = analysis_mode,
    scan_mode = scan_mode,
    step_bp = if (scan_mode == "window") as.integer(step_bp) else NA_integer_,
    window_agg = if (scan_mode == "window") window_agg else NA_character_,
    marker_map = map_use,
    res_all = res_all,
    res_by_family = res_by_family
  )
}

## ============================================================
## HOW YOU CALL IT (THIS is the “family-wise instead of combined”)
## ============================================================

## A) BOTH (combined + familywise)
fit <- run_bzea_gridlmm_scan(
  pheno_file = "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv",
  trait_col  = "DTS",
  id_col     = 1,
  geno_file  = "RTiger_binned_dosage_genome_step1e+05.csv",
  scan_mode  = "window",
  step_bp    = 100000L,
  window_agg = "mean",
  K_dir      = "/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K/tmp",
  K_format   = "plink2_rel",
  analysis_mode = "both",
  include_family_fixed = TRUE,
  cores = 8
)
str(fit)
namefitnames(fit$res_by_family)
head(fit$res_by_family$Zd[, c("Marker","chr","pos","p","mlog10p")])

## B) FAMILYWISE ONLY
fit_fam <- run_bzea_gridlmm_scan(
  pheno_file = "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv",
  trait_col  = "DTS",
  id_col     = 1,
  geno_file  = "RTiger_markerlevel_dosage.csv",
  scan_mode  = "single",
  K_dir      = "/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K/tmp",
  K_format   = "plink2_rel",
  analysis_mode = "familywise",
  cores = 8
)
names(fit_fam$res_by_family)
lapply(fit_fam$res_by_family, function(df) head(df[, c("Family","Marker","p","mlog10p")], 3))



