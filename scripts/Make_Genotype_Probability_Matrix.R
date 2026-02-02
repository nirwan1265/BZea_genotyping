library(dplyr)
library(data.table)
library(stringr)

obj <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/RTIGER_RDS/Samples_chr1_2_1_100_RTIGER_results.rds")
obj2 <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/RTIGER_RDS/Samples_chr1_2_101_200_RTIGER_results.rds")
#obj3 <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/RTIGER_RDS/Samples_chr1_2_201_300_RTIGER_results.rds")
obj3 <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/RTIGER_RDS/Samples_chr1_2_301_400_RTIGER_results.rds")
obj4 <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/RTIGER_RDS/Samples_chr1_2_401_500_RTIGER_results.rds")
obj5 <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/RTIGER_RDS/Samples_chr1_2_501_611_RTIGER_results.rds")


gamma <- obj@Probabilities[["gamma"]]
gamma2 <- obj2@Probabilities[["gamma"]]
gamma3 <- obj3@Probabilities[["gamma"]]
gamma4 <- obj4@Probabilities[["gamma"]]
gamma5 <- obj5@Probabilities[["gamma"]]
#gamma6 <- obj6@Probabilities[["gamma"]]




# expDesign is a list, each element is a data.frame (or a list of dfs)
# from your str: list [100 x 3] (S3: data.frame)
ed <- obj@info$expDesign
ed2 <- obj2@info$expDesign
ed3 <- obj3@info$expDesign
ed4 <- obj4@info$expDesign
ed5 <- obj5@info$expDesign


# If expDesign is a list of data.frames, bind them:
ed_df <- ed
ed_df2 <- ed2
ed_df3 <- ed3
ed_df4 <- ed4
ed_df5 <- ed5

# keep only what we need
ed_df <- ed_df[, c("name", "OName")]
ed_df2 <- ed_df2[, c("name", "OName")]
ed_df3 <- ed_df3[, c("name", "OName")]
ed_df4 <- ed_df4[, c("name", "OName")]
ed_df5 <- ed_df5[, c("name", "OName")]

# build mapping: Sample_# -> real sample ID
map <- setNames(ed_df$OName, ed_df$name)
map2 <- setNames(ed_df2$OName, ed_df2$name)
map3 <- setNames(ed_df3$OName, ed_df3$name)
map4 <- setNames(ed_df4$OName, ed_df4$name)
map5 <- setNames(ed_df5$OName, ed_df5$name)

# sanity check: do all gamma names exist in map?
missing <- setdiff(names(gamma), names(map))
if (length(missing) > 0) {
  stop("These gamma sample names are not in expDesign$name: ",
       paste(missing, collapse = ", "))
}

# rename safely using the mapping (order-independent)
names(gamma) <- unname(map[names(gamma)])
names(gamma2) <- unname(map2[names(gamma2)])
names(gamma3) <- unname(map3[names(gamma3)])
names(gamma4) <- unname(map4[names(gamma4)])
names(gamma5) <- unname(map5[names(gamma5)])

# assign back
obj@Probabilities[["gamma"]] <- gamma
obj2@Probabilities[["gamma"]] <- gamma2
obj3@Probabilities[["gamma"]] <- gamma3
obj4@Probabilities[["gamma"]] <- gamma4
obj5@Probabilities[["gamma"]] <- gamma5

# optional: also rename alpha/beta/psi consistently
for (slot in c("alpha","beta","psi")) {
  x <- obj@Probabilities[[slot]]
  if (!is.null(x)) {
    miss2 <- setdiff(names(x), names(map))
    if (length(miss2) == 0) {
      names(x) <- unname(map[names(x)])
      obj@Probabilities[[slot]] <- x
    }
  }
}

for (slot in c("alpha","beta","psi") ) {
  x2 <- obj2@Probabilities[[slot]]
  if (!is.null(x2)) {
    miss2 <- setdiff(names(x2), names(map2))
    if (length(miss2) == 0) {
      names(x2) <- unname(map2[names(x2)])
      obj2@Probabilities[[slot]] <- x2
    }
  }
}

for (slot in c("alpha","beta","psi") ) {
  x3 <- obj3@Probabilities[[slot]]
  if (!is.null(x3)) {
    miss3 <- setdiff(names(x3), names(map3))
    if (length(miss3) == 0) {
      names(x3) <- unname(map3[names(x3)])
      obj3@Probabilities[[slot]] <- x3
    }
  }
}

for (slot in c("alpha","beta","psi") ) {
  x4 <- obj4@Probabilities[[slot]]
  if (!is.null(x4)) {
    miss4 <- setdiff(names(x4), names(map4))
    if (length(miss4) == 0) {
      names(x4) <- unname(map4[names(x4)])
      obj4@Probabilities[[slot]] <- x4
    }
  }
}

for (slot in c("alpha","beta","psi") ) {
  x5 <- obj5@Probabilities[[slot]]
  if (!is.null(x5)) {
    miss5 <- setdiff(names(x5), names(map5))
    if (length(miss5) == 0) {
      names(x5) <- unname(map5[names(x5)])
      obj5@Probabilities[[slot]] <- x5
    }
  }
}



# quick check: show first few names
head(names(obj@Probabilities[["gamma"]]))
head(names(obj2@Probabilities[["gamma"]]))
head(names(obj3@Probabilities[["gamma"]]))
head(names(obj4@Probabilities[["gamma"]]))
head(names(obj5@Probabilities[["gamma"]]))


# save
saveRDS(obj, "Samples_chr1_2_RTIGER_results_renamed.rds")




# Get df 
get_file_map <- function(obj) {
  # returns named character vector: names = sample IDs (OName), values = TSV file paths
  
  if ("expDesign" %in% slotNames(obj)) {
    ed <- obj@expDesign
    ed_df <- if (is.data.frame(ed)) ed else do.call(rbind, ed)
    stopifnot(all(c("OName","files") %in% names(ed_df)))
    return(setNames(as.character(ed_df$files), as.character(ed_df$OName)))
  }
  
  if ("info" %in% slotNames(obj)) {
    inf <- obj@info
    
    # common patterns (based on your screenshot earlier)
    if (is.list(inf) && all(c("OName","files") %in% names(inf))) {
      return(setNames(as.character(inf$files), as.character(inf$OName)))
    }
    
    # sometimes info is a list of fields, not a data.frame
    if (is.list(inf) && !is.null(inf$expDesign)) {
      ed <- inf$expDesign
      ed_df <- if (is.data.frame(ed)) ed else do.call(rbind, ed)
      stopifnot(all(c("OName","files") %in% names(ed_df)))
      return(setNames(as.character(ed_df$files), as.character(ed_df$OName)))
    }
  }
  
  stop("Could not find TSV file paths in obj. I looked in obj@expDesign and obj@info.")
}

file_map <- get_file_map(obj)
file_map2 <- get_file_map(obj2)
file_map3 <- get_file_map(obj3)
file_map4 <- get_file_map(obj4)
file_map5 <- get_file_map(obj5)

head(file_map)
head(file_map2)
head(file_map3)

read_pos_one <- function(tsv_file, chr_key) {
  tb <- fread(tsv_file)
  nm <- names(tb); nml <- tolower(nm)
  
  chr_col <- nm[match(TRUE, nml %in% c("chr","chrom","chromosome"), nomatch = 0)]
  pos_col <- nm[match(TRUE, nml %in% c("pos","position","bp"), nomatch = 0)]
  if (length(chr_col) == 0) chr_col <- nm[1]
  if (length(pos_col) == 0) pos_col <- nm[2]
  
  chr <- as.character(tb[[chr_col[1]]])
  chr <- ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr))
  pos <- as.integer(tb[[pos_col[1]]])
  
  pos[chr == chr_key]
}

make_bin_map <- function(chr_key, chr_len, step_bp) {
  nb <- as.integer(ceiling(chr_len / step_bp))
  start <- (seq_len(nb) - 1L) * step_bp + 1L
  end   <- pmin(seq_len(nb) * step_bp, chr_len)
  mid   <- as.integer(floor((start + end) / 2))
  marker <- paste0(chr_key, ":", mid)
  data.frame(chr = chr_key, start = start, end = end, pos = mid, marker = marker)
}

bin_one_sample <- function(dos, pos, chr_len, step_bp) {
  nb <- as.integer(ceiling(chr_len / step_bp))
  bin_id <- as.integer((pos - 1L) %/% step_bp) + 1L
  bin_id <- pmin(pmax(bin_id, 1L), nb)
  
  out <- rep(NA_real_, nb)
  ok <- !is.na(dos) & !is.na(pos)
  if (!any(ok)) return(out)
  
  sumv <- tapply(dos[ok], bin_id[ok], sum)
  cntv <- tapply(rep(1L, sum(ok)), bin_id[ok], sum)
  
  idx <- as.integer(names(sumv))
  out[idx] <- as.numeric(sumv) / as.numeric(cntv)
  out
}

make_binned_dosage_chr <- function(obj, chr, chr_len, step_bp = 100000,
                                   iHET = 2, iTEO = 3) {
  
  gamma <- obj@Probabilities[["gamma"]]
  file_map <- get_file_map(obj)
  
  taxa <- intersect(names(gamma), names(file_map))
  stopifnot(length(taxa) > 0)
  
  bins <- make_bin_map(chr, chr_len, step_bp)
  
  Xbin <- matrix(NA_real_, nrow = length(taxa), ncol = nrow(bins),
                 dimnames = list(taxa, bins$marker))
  
  for (k in seq_along(taxa)) {
    id <- taxa[k]
    g  <- gamma[[id]][[chr]]
    if (is.null(g)) next
    
    dos <- as.numeric(g[iHET, ]) + 2 * as.numeric(g[iTEO, ])
    
    tsv <- file_map[[id]]
    pos <- read_pos_one(tsv, chr)
    
    if (length(pos) != length(dos)) {
      warning("Mismatch ", id, " ", chr,
              ": pos=", length(pos), " gamma=", length(dos),
              " (skipping)")
      next
    }
    
    Xbin[k, ] <- bin_one_sample(dos, pos, chr_len, step_bp)
  }
  
  list(X = Xbin, map = bins)
}


chr_len <- readRDS("/Users/nirwantandukar/Documents/Research/data/BZea/genotype/chr_len.rds")


# Single Chromosome
# out1 <- make_binned_dosage_chr(
#   obj     = obj,
#   chr     = "chr1",
#   chr_len = chr_len[["chr1"]],
#   step_bp = 100000,
#   iHET    = 2,
#   iTEO    = 3
# )
# geno_chr1 <- out1$X
# dim(geno_chr1)
# 
# geno_df <- data.frame(Line = rownames(geno_chr1), geno_chr1, check.names = FALSE)
# fwrite(geno_df, "RTiger_binned_dosage_chr1_step100kb.csv")



# All chromosomes
step_bp <- 10000
allX <- list()

for (cc in paste0("chr", 1:2)) {
  out <- make_binned_dosage_chr(
    obj     = obj,
    chr     = cc,
    chr_len = chr_len[[cc]],
    step_bp = step_bp,
    iHET    = 2,
    iTEO    = 3
  )
  allX[[cc]] <- out$X
}

geno_mat <- do.call(cbind, allX)   # rows=samples, cols=chr#:pos bins
geno_df  <- data.frame(Line = rownames(geno_mat), geno_mat, check.names = FALSE)

allX2 <- list()
for (cc in paste0("chr", 1:2)) {
  out2 <- make_binned_dosage_chr(
    obj     = obj2,
    chr     = cc,
    chr_len = chr_len[[cc]],
    step_bp = step_bp,
    iHET    = 2,
    iTEO    = 3
  )
  allX2[[cc]] <- out2$X
}

geno_mat2 <- do.call(cbind, allX2)   # rows=samples, cols=chr#:pos bins
geno_df2  <- data.frame(Line = rownames(geno_mat2), geno_mat2, check.names = FALSE)


allX3 <- list()
for (cc in paste0("chr", 1:2)) {
  out3 <- make_binned_dosage_chr(
    obj     = obj3,
    chr     = cc,
    chr_len = chr_len[[cc]],
    step_bp = step_bp,
    iHET    = 2,
    iTEO    = 3
  )
  allX3[[cc]] <- out3$X
}

geno_mat3 <- do.call(cbind, allX3)   # rows=samples, cols=chr#:pos bins
geno_df3  <- data.frame(Line = rownames(geno_mat3), geno_mat3, check.names = FALSE)


allX4 <- list()
for (cc in paste0("chr", 1:2)) {
  out4 <- make_binned_dosage_chr(
    obj     = obj4,
    chr     = cc,
    chr_len = chr_len[[cc]],
    step_bp = step_bp,
    iHET    = 2,
    iTEO    = 3
  )
  allX4[[cc]] <- out4$X
}
geno_mat4 <- do.call(cbind, allX4)   # rows=samples, cols=chr#:pos bins
geno_df4  <- data.frame(Line = rownames(geno_mat4), geno_mat4, check.names = FALSE)

allX5 <- list()
for (cc in paste0("chr", 1:2)) {
  out5 <- make_binned_dosage_chr(
    obj     = obj5,
    chr     = cc,
    chr_len = chr_len[[cc]],
    step_bp = step_bp,
    iHET    = 2,
    iTEO    = 3
  )
  allX5[[cc]] <- out5$X
}
geno_mat5 <- do.call(cbind, allX5)   # rows=samples, cols=chr#:pos bins
geno_df5  <- data.frame(Line = rownames(geno_mat5), geno_mat5, check.names = FALSE)

# Combine the geno_dfs 
geno_df <- bind_rows(geno_df, geno_df2, geno_df3, geno_df4, geno_df5)

fwrite(geno_df, sprintf("RTiger_binned_dosage_genome_step%s_chr1_2.csv", step_bp))

