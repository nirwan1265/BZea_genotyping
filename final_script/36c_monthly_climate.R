#!/usr/bin/env Rscript
## =============================================================================
## 36c — Monthly climate at the collection site of every donor accession: WorldClim 2.1 monthly precipitation
##       (prec, mm), mean temperature (tavg, deg C) and solar radiation (srad, kJ m-2 day-1), 10 arc-min, the same source and resolution as BIO1-19 in
##       data/env/BZea_accession_env.csv. For step 36d (which months of rain go with a QTL's teosinte effect).
##
## Run on the laptop, from the repo root:
##   Rscript final_script/36c_monthly_climate.R
## Coordinates: data/env/BZea_accession_env.csv (accession-id, latitude, longitude; longitudes already corrected
##   to west, see old_scripts/33_env_extract.R).
## Rasters: WC_DIR (default ~/Documents/Github/Landadapt/data/WorldClim/cache/climate/wc2.1_10m), files
##   wc2.1_10m_<var>_01.tif ... _12.tif for prec, tavg and srad. If they are not there, they are
##   downloaded once into that folder with geodata::worldclim_global() (about 10 MB per variable; needs internet).
## Output: data/env/BZea_accession_monthly.csv (accession-id, taxa_code, race, prec_01..12, tavg_01..12, srad_01..12)
## Needs: R with terra and data.table (and geodata if the monthly rasters must be downloaded).
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra) })
gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) path.expand(v) else path.expand(d) }
WC <- gv("WC_DIR", "~/Documents/Github/Landadapt/data/WorldClim/cache/climate/wc2.1_10m")
ENV <- "data/env/BZea_accession_env.csv"
OUT <- "data/env/BZea_accession_monthly.csv"
if (!file.exists(ENV)) stop("run from the repo root (", ENV, " not found)")

files <- function(v) file.path(WC, sprintf("wc2.1_10m_%s_%02d.tif", v, 1:12))
for (v in c("prec", "tavg", "srad")) {
  if (!all(file.exists(files(v)))) {
    if (!requireNamespace("geodata", quietly = TRUE)) stop("monthly ", v, " rasters not in ", WC,
                                                          "; install geodata (install.packages('geodata')) to download them")
    message("downloading WorldClim 2.1 monthly ", v, " (10 arc-min) into ", WC)
    geodata::worldclim_global(var = v, res = 10, path = dirname(dirname(WC)))   # geodata writes <path>/climate/wc2.1_10m
  }
  if (!all(file.exists(files(v)))) stop("monthly ", v, " rasters still missing in ", WC)
}

e <- fread(ENV)
pts <- vect(as.data.frame(e[, .(longitude, latitude)]), geom = c("longitude", "latitude"), crs = "EPSG:4326")
out <- e[, .(`accession-id`, taxa_code, race)]
for (v in c("prec", "tavg", "srad")) {
  r <- rast(files(v))
  x <- extract(r, pts, ID = FALSE)
  setnames(x, sprintf("%s_%02d", v, 1:12))
  out <- cbind(out, as.data.table(x))
}
fwrite(out, OUT)
cat("wrote", OUT, ":", nrow(out), "accessions;", sum(complete.cases(out)), "with all 36 values\n")
print(out[, lapply(.SD, function(z) round(median(z, na.rm = TRUE), 1)), .SDcols = sprintf("prec_%02d", 1:12)])
