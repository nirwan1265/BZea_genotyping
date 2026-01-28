suppressPackageStartupMessages({
  library(data.table)
})

# ---------------- SETTINGS ----------------
state_dir <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/statepaths"
map_dir   <- "/rsstu/users/r/rrellan/sara/ref/NAM_genetic_map"

out_dir   <- "/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/introgression_cmBridge"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Bridge any small gap (RR or RH) INSIDE donor if its total genetic length <= this
gap_cM <- 0.02   # try 0.02, 0.05, 0.10

# Define donor states (for continuity)
# If you want donor = HH only, set donor_set <- "HH"
donor_set <- c("HH","RH")

# Output files
out_summary <- file.path(out_dir, "all_samples.ref_alt_het_perc.cmBridge.tsv")
out_all_tracts <- file.path(out_dir, "all_samples.alt_het.tracts.cmBridge.tsv.gz")

# chromosome sizes (bp) you gave
chr_len <- c(308452471,243675191,238017767,250330460,226353449,
             181357234,185808916,182411202,163004744,152435371)
names(chr_len) <- paste0("chr", 1:10)
genome_bp <- sum(chr_len)

# ---------------- MAPS (LOAD ONCE) ----------------
maps <- lapply(paste0("chr",1:10), function(chr) {
  f <- file.path(map_dir, paste0(chr, ".map"))
  m <- fread(f, header=FALSE)
  # chrN.map = bp  1.0  cM
  setnames(m, c("bp","col2","cM"))
  setorder(m, bp)
  m
})
names(maps) <- paste0("chr",1:10)

interp_cm_vec <- function(chr, pos) {
  m <- maps[[chr]]
  # vectorized linear interpolation; rule=2 clamps outside range
  approx(x=m$bp, y=m$cM, xout=pos, rule=2)$y
}

# ---------------- RUN SEGMENTS ----------------
make_runs <- function(state, w_cM, pos, chr) {
  r <- rle(state)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  seg <- data.table(
    run_id = seq_along(r$values),
    state  = r$values,
    start_i = starts,
    end_i   = ends,
    n_sites = r$lengths
  )

  seg[, Chr := chr[start_i]]
  seg[, Start := pos[start_i]]
  seg[, End   := pos[end_i]]

  # run length in cM = sum of interval weights across SNPs in this run
  seg[, len_cM := vapply(seq_len(.N), function(i) {
    sum(w_cM[start_i[i]:end_i[i]], na.rm=TRUE)
  }, numeric(1))]

  # run length in bp (simple boundary span)
  seg[, len_bp := as.numeric(End - Start + 1L)]
  seg[]
}

# ---------------- BRIDGE GAPS INSIDE DONOR ----------------
# We build continuous donor tracts by allowing small "gap" runs (RR or RH)
# between donor runs, as long as total gap_cM <= threshold.
bridge_donor_runs <- function(seg, donor_set, gap_cM) {
  if (nrow(seg) == 0) return(data.table())

  out <- list()
  k <- 0L
  i <- 1L

  while (i <= nrow(seg)) {
    # start a donor tract only when we hit donor
    if (!(seg$state[i] %in% donor_set)) {
      i <- i + 1L
      next
    }

    chr <- seg$Chr[i]
    j <- i

    # extend through donor, and through small gaps that are flanked by donor
    while (TRUE) {
      # advance through consecutive donor runs
      while (j + 1L <= nrow(seg) && seg$Chr[j+1L] == chr && (seg$state[j+1L] %in% donor_set)) {
        j <- j + 1L
      }

      # now consider a gap block after j
      g_start <- j + 1L
      if (g_start > nrow(seg) || seg$Chr[g_start] != chr) break

      # compute total gap_cM across one or more non-donor runs until next donor
      g_end <- g_start
      gap_sum <- 0.0

      while (g_end <= nrow(seg) && seg$Chr[g_end] == chr && !(seg$state[g_end] %in% donor_set)) {
        gap_sum <- gap_sum + seg$len_cM[g_end]
        g_end <- g_end + 1L
        # early break if already too big
        if (gap_sum > gap_cM) break
      }

      # g_end is first run AFTER the gaps
      if (gap_sum <= gap_cM &&
          g_end <= nrow(seg) &&
          seg$Chr[g_end] == chr &&
          (seg$state[g_end] %in% donor_set)) {
        # bridge: include gaps + next donor block
        j <- g_end
        next
      } else {
        break
      }
    }

    # build tract i..j (includes bridged gaps)
    k <- k + 1L
    tract_runs <- seg[i:j]

    # label tract as Alt if ANY HH appears; else Het
    region <- if (any(tract_runs$state == "HH")) "Alt" else "Het"

    out[[k]] <- data.table(
      Chr = chr,
      Start = min(tract_runs$Start),
      End   = max(tract_runs$End),
      Region = region,
      n_sites = sum(tract_runs$n_sites),
      len_bp = as.numeric(max(tract_runs$End) - min(tract_runs$Start) + 1L),
      len_cM = sum(tract_runs$len_cM)
    )

    i <- j + 1L
  }

  rbindlist(out, fill=TRUE)
}

# gzip writer (portable)
write_tsv_gz <- function(dt, out_gz) {
  tmp <- tempfile(fileext = ".tsv")
  fwrite(dt, tmp, sep="\t")
  system2("gzip", c("-c", tmp), stdout = out_gz)
  unlink(tmp)
}

# ---------------- PER SAMPLE ----------------
files <- sort(list.files(state_dir, pattern="\\.statepath\\.tsv\\.gz$", full.names=TRUE))
if (length(files) == 0) stop("No *.statepath.tsv.gz in state_dir")

all_tracts <- vector("list", length(files))
all_stats  <- vector("list", length(files))

for (ii in seq_along(files)) {
  f <- files[ii]
  sample <- sub("\\.statepath\\.tsv\\.gz$", "", basename(f))

  dt <- fread(cmd = paste("zcat", shQuote(f)))
  # Expect at least: CHROM POS STATE
  # Make it robust to different column names
  setnames(dt, old=intersect(names(dt), c("CHROM","Chr")), new="CHROM", skip_absent=TRUE)
  setnames(dt, old=intersect(names(dt), c("POS","Pos")),   new="POS",   skip_absent=TRUE)
  setnames(dt, old=intersect(names(dt), c("STATE","State")), new="STATE", skip_absent=TRUE)

  dt <- dt[, .(CHROM, POS, STATE)]
  dt <- dt[CHROM %in% names(chr_len)]
  setorder(dt, CHROM, POS)

  # split per chr for speed and correct runs
  tr_list <- list()

  for (chr in names(chr_len)) {
    dchr <- dt[CHROM == chr]
    if (nrow(dchr) == 0) next

    # interpolate cM vector fast
    dchr[, cM := interp_cm_vec(chr, POS)]

    # interval weights in cM to next SNP within chr
    dchr[, next_cM := shift(cM, type="lead")]
    dchr[, w_cM := pmax(0, next_cM - cM)]
    dchr[is.na(w_cM), w_cM := 0]

    seg <- make_runs(dchr$STATE, dchr$w_cM, dchr$POS, dchr$CHROM)
    tr_chr <- bridge_donor_runs(seg, donor_set=donor_set, gap_cM=gap_cM)

    if (nrow(tr_chr) > 0) tr_list[[chr]] <- tr_chr
  }

  tr <- rbindlist(tr_list, fill=TRUE)
  if (nrow(tr) == 0) {
    tr <- data.table(Chr=character(), Start=integer(), End=integer(), Region=character(),
                     n_sites=integer(), len_bp=numeric(), len_cM=numeric())
  }

  tr[, Sample := sample]
  setcolorder(tr, c("Sample","Chr","Start","End","Region","n_sites","len_bp","len_cM"))

  # Per-sample output tracts
  out_f <- file.path(out_dir, paste0(sample, ".alt_het.tracts.cmBridge.tsv.gz"))
  write_tsv_gz(tr, out_f)

  # Percentages (bp-based) from bridged tracts
  Alt_bp <- sum(tr$len_bp[tr$Region=="Alt"], na.rm=TRUE)
  Het_bp <- sum(tr$len_bp[tr$Region=="Het"], na.rm=TRUE)
  Ref_bp <- genome_bp - Alt_bp - Het_bp
  if (Ref_bp < 0) Ref_bp <- 0

  st <- data.table(
    Sample = sample,
    Ref_perc = 100 * Ref_bp / genome_bp,
    Alt_perc = 100 * Alt_bp / genome_bp,
    Het_perc = 100 * Het_bp / genome_bp,
    Ref_bp = Ref_bp,
    Alt_bp = Alt_bp,
    Het_bp = Het_bp
  )

  all_tracts[[ii]] <- tr
  all_stats[[ii]]  <- st

  if (ii %% 50 == 0) message("Done ", ii, " / ", length(files))
}

tract_all <- rbindlist(all_tracts, fill=TRUE)
stats_all <- rbindlist(all_stats,  fill=TRUE)

# TOTAL row (pooled across samples)
nS <- nrow(stats_all)
tot_alt <- sum(stats_all$Alt_bp)
tot_het <- sum(stats_all$Het_bp)
tot_ref <- (genome_bp * nS) - tot_alt - tot_het
if (tot_ref < 0) tot_ref <- 0

total_row <- data.table(
  Sample="TOTAL",
  Ref_perc = 100 * tot_ref / (genome_bp * nS),
  Alt_perc = 100 * tot_alt / (genome_bp * nS),
  Het_perc = 100 * tot_het / (genome_bp * nS),
  Ref_bp = tot_ref, Alt_bp = tot_alt, Het_bp = tot_het
)

stats_out <- rbind(stats_all, total_row, fill=TRUE)
fwrite(stats_out, out_summary, sep="\t")

write_tsv_gz(tract_all, out_all_tracts)

cat("Wrote:\n")
cat("  Tracts per sample: ", out_dir, "\n", sep="")
cat("  Combined tracts:   ", out_all_tracts, "\n", sep="")
cat("  Summary perc:      ", out_summary, "\n", sep="")
