#!/usr/bin/env Rscript
## =============================================================================
## 33 — Environment at the collection site of every donor accession (climate, aridity, soil, mycorrhiza,
##      elevation), from rasters, for analyses of the teosinte k-mer profiles and QTL effects against the
##      donor's home environment. Adapted from Landadapt scripts/geno_env/08_extract_env_all_species.R.
##
## Run on the laptop that has the rasters, from the repo root:
##   Rscript final_script/33_env_extract.R
## Input:  data/metadata_BZea.csv (accession-id, taxa_code, race, latitude, longitude, elevation, country)
## Output: data/env/BZea_accession_env.csv   one row per accession: metadata + one column per layer
##         data/env/BZea_env_summary.tsv      per layer: accessions with a value, min, median, max
##
## Coordinates: all teosinte sites are in Mexico or Guatemala, west of Greenwich, so their longitude must be
## negative. In the metadata 75 of the 126 teosinte longitudes are positive (the minus sign is missing); every
## teosinte (and any Mexico/Guatemala row) gets -|longitude|, and the number fixed is printed. Landraces (Zm) are left out by default (KEEP_ZM=1
## keeps them; their coordinates are used as they are).
##
## Layers (value of the raster cell at the point; terra::extract, method "simple"):
##   WorldClim 2.1 BIO1-BIO19, 10 arc-min, raw units (BIO1 mean annual temperature in deg C, BIO12 annual
##     precipitation in mm, ...), and the WorldClim elevation (wc_elev)
##   aridity: Global Aridity Index v3 (annual), / 10000; 0 (sea, coast, no data) -> NA. et0: reference ET (mm/yr)
##   SoilGrids 0-5 cm, 1 km, raw units: nitrogen (cg/kg), phh2o (pH x 10), cec, soc, ocd, clay, sand, silt,
##     bdod, cfvo
##   AM / EcM roots colonized and relative abundance (Barcelo et al.)
##   elevation (ELEVATION.tif), NPP (MANPP.tif), mean annual precipitation (MAPRECIP.tif)
## Paths (environment, defaults as on Nirwan's laptop):
##   PRED_DIR  ~/Documents/Research/data/P_prediction/predictors
##   AI_DIR    ~/Documents/Research/data/P_prediction/Global-AI_ET0_v3_annual
##   WC_DIR    ~/Documents/Github/Landadapt/data/WorldClim/cache/climate/wc2.1_10m
##   EXTRA     more layers, "name=path,name=path" (e.g. other soil depths)
## Needs: R with terra and data.table.
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra) })
gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) path.expand(v) else path.expand(d) }
PRED <- gv("PRED_DIR", "~/Documents/Research/data/P_prediction/predictors")
AI   <- gv("AI_DIR", "~/Documents/Research/data/P_prediction/Global-AI_ET0_v3_annual")
WC   <- gv("WC_DIR", "~/Documents/Github/Landadapt/data/WorldClim/cache/climate/wc2.1_10m")
REPO <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
META <- file.path(REPO, "data", "metadata_BZea.csv")
OUT  <- file.path(REPO, "data", "env"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

soil <- c("nitrogen", "phh2o", "cec", "soc", "ocd", "clay", "sand", "silt", "bdod", "cfvo")
RAS <- c(setNames(file.path(WC, sprintf("wc2.1_10m_bio_%d.tif", 1:19)), paste0("BIO", 1:19)),
         wc_elev = file.path(WC, "wc2.1_10m_elev.tif"),
         aridity = file.path(AI, "ai_v3_yr.tif"),
         et0     = file.path(AI, "et0_v3_yr.tif"),
         setNames(file.path(PRED, sprintf("%s_0-5cm_mean_1000.tif", soil)), paste0(soil, "_0_5cm")),
         AM_roots  = file.path(PRED, "AM_roots_colonized.tif"),
         AM_rel    = file.path(PRED, "AM_rel.abundance_colonization.tif"),
         EcM_roots = file.path(PRED, "all_predictors", "EcM_roots_colonized.tif"),
         EcM_rel   = file.path(PRED, "all_predictors", "EcM_rel.abundance_colonization.tif"),
         elevation_raster = file.path(PRED, "ELEVATION.tif"),
         NPP       = file.path(PRED, "MANPP.tif"),
         MAP       = file.path(PRED, "MAPRECIP.tif"))
ex <- Sys.getenv("EXTRA")
if (nzchar(ex)) for (kv in strsplit(ex, ",")[[1]]) { p <- strsplit(kv, "=")[[1]]; RAS[p[1]] <- path.expand(p[2]) }

m <- fread(META, encoding = "UTF-8")
setnames(m, 1, sub("^﻿", "", names(m)[1]))
if (Sys.getenv("KEEP_ZM") != "1") m <- m[taxa_code != "Maiz"]
fix <- (m$country %in% c("México", "Mexico", "Guatemala") | m$taxa_code != "Maiz") & !is.na(m$longitude) & m$longitude > 0
m[fix, longitude := -longitude]
cat(sprintf("%d accessions; longitude sign fixed for %d (Mexico/Guatemala given as east of Greenwich)\n", nrow(m), sum(fix)))
ok <- !is.na(m$latitude) & !is.na(m$longitude)
pts <- vect(as.data.frame(m[ok, .(longitude, latitude)]), geom = c("longitude", "latitude"), crs = "EPSG:4326")

for (nm in names(RAS)) {
  f <- RAS[[nm]]
  if (!file.exists(f)) { warning("not found, column left empty: ", f); set(m, j = nm, value = NA_real_); next }
  r <- rast(f)
  p <- if (same.crs(r, pts)) pts else project(pts, crs(r))
  v <- extract(r, p, ID = FALSE)[[1]]
  col <- rep(NA_real_, nrow(m)); col[ok] <- as.numeric(v)
  set(m, j = nm, value = col)
  cat(sprintf("  %-18s %3d / %d with a value\n", nm, sum(!is.na(col[ok])), sum(ok)))
}
if ("aridity" %in% names(m)) { m[aridity == 0, aridity := NA_real_]; m[, aridity := aridity / 10000] }

fwrite(m, file.path(OUT, "BZea_accession_env.csv"))
sm <- rbindlist(lapply(names(RAS), function(nm) { x <- m[[nm]]; x <- x[!is.na(x)]
  data.table(layer = nm, accessions = length(x), min = if (length(x)) min(x) else NA,
             median = if (length(x)) median(x) else NA, max = if (length(x)) max(x) else NA) }))
fwrite(sm, file.path(OUT, "BZea_env_summary.tsv"), sep = "\t")
if (sum(complete.cases(m$elevation, m$wc_elev)) > 2)
  cat(sprintf("elevation check (metadata vs WorldClim): r = %.2f\n", cor(m$elevation, m$wc_elev, use = "complete.obs")))
cat("written to", OUT, "\n")
