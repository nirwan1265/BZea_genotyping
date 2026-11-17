#!/usr/bin/env Rscript
## =============================================================================
## 32d — The four mexicana races (Chalco, Mesa Central, Durango, Nobogame) in the teosinte k-mer profiles (32b).
##   A: where the accessions were collected (longitude / latitude; parviglumis in grey)
##   B: PCA of the mexicana pools' containment of each assembly (8 assemblies, centred log-ratio, prcomp)
##   C: per pool, log(containment of TIL18 / containment of TIL25), the two mexicana assemblies, by race
##   D: mexicana index f_Zx / (f_Zx + f_Zv) per accession against elevation (parviglumis in grey)
## Usage:  Rscript 32d_mexicana_races.R <kmer_taxon_profile.tsv> <metadata_BZea.csv> <out dir>
## Outputs: mexicana_races.png, mexicana_races_pools.tsv
## Teosinte longitudes are set to -|longitude| (Mexico and Guatemala; the sign is missing in the metadata).
## Needs: R with data.table, ggplot2, patchwork.
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); if (length(a) < 3) stop("usage: Rscript 32d_mexicana_races.R <kmer_taxon_profile.tsv> <metadata> <out dir>")
d <- fread(a[1]); out <- a[3]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
m <- fread(a[2], encoding = "UTF-8"); setnames(m, 1, sub("^﻿", "", names(m)[1]))
m[taxa_code != "Maiz" & longitude > 0, longitude := -longitude]
rc <- c(Chal = "Chalco", Mesa = "Mesa Central", Dura = "Durango", Nobo = "Nobogame", Bals = "parviglumis (Balsas)")
cols <- c(Chalco = "#08306B", `Mesa Central` = "#2171B5", Durango = "#6BAED6", Nobogame = "#9E9AC8",
          `parviglumis (Balsas)` = "grey70")
d <- merge(d[role == "donor" & label %in% c("Zx", "Zv")],
           m[, .(accession = `accession-id`, group = rc[taxa_code], elevation, latitude, longitude)], by = "accession")
d <- d[!is.na(group)]
d[, group := factor(group, levels = names(cols))]
x <- d[label == "Zx"]

## A: sites
st <- unique(d[, .(accession, group, latitude, longitude, elevation)])
pA <- ggplot(st[order(group != "parviglumis (Balsas)")], aes(longitude, latitude, colour = group)) +
  geom_point(size = 2.2, alpha = 0.9) + scale_colour_manual(values = cols, name = NULL, drop = FALSE) + coord_quickmap() +
  labs(x = "Longitude", y = "Latitude", title = "A. Collection sites") + theme_classic(base_size = 11)
## B: PCA of per-assembly containment, mexicana pools
cc <- grep("^c_", names(x), value = TRUE)
X <- as.matrix(x[, ..cc]); X <- X / rowSums(X); L <- log(pmax(X, min(X[X > 0]) / 2)); L <- L - rowMeans(L)
pc <- prcomp(L); ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
x[, `:=`(PC1 = pc$x[, 1], PC2 = pc$x[, 2])]
pB <- ggplot(x, aes(PC1, PC2, colour = group)) + geom_point(size = 1.3, alpha = 0.75) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = sprintf("B. PCA of %d mexicana pools (containment of 8 assemblies)", nrow(x))) + theme_classic(base_size = 11)
## C: TIL18 vs TIL25
x[, til := log(c_Zx_TIL18 / c_Zx_TIL25)]
pC <- ggplot(x, aes(group, til, colour = group)) + geom_hline(yintercept = 0, colour = "grey70") +
  geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_colour_manual(values = cols, guide = "none") + scale_x_discrete(drop = TRUE) +
  labs(x = NULL, y = "log(containment TIL18 / TIL25)", title = "C. Closer to mexicana assembly TIL18 (+) or TIL25 (-)") +
  theme_classic(base_size = 11)
## D: mexicana index along elevation
d[, mx := f_Zx / (f_Zx + f_Zv)]
ac <- d[is.finite(mx), .(mx = median(mx), pools = .N), by = .(accession, group, elevation)]
pD <- ggplot(ac[order(group != "parviglumis (Balsas)")], aes(elevation, mx, colour = group, size = pools)) +
  geom_point(alpha = 0.85) + scale_colour_manual(values = cols, guide = "none") + scale_size_area(max_size = 4, name = "Pools") +
  labs(x = "Elevation of the accession (m)", y = "f_Zx / (f_Zx + f_Zv), median per accession",
       title = "D. mexicana index along elevation") + theme_classic(base_size = 11)
p <- (pA | pB) / (pC | pD) + plot_layout(guides = "collect")
ggsave(file.path(out, "mexicana_races.png"), p, width = 13, height = 9, dpi = 200)
fwrite(x[, .(sample, accession, race = group, elevation, PC1, PC2, log_TIL18_TIL25 = til, f_Zx, f_Zv)],
       file.path(out, "mexicana_races_pools.tsv"), sep = "\t")
print(x[, .(pools = .N, accessions = uniqueN(accession), median_elevation = as.numeric(median(elevation)),
            median_log_TIL18_TIL25 = round(median(til), 3), median_PC1 = round(median(PC1), 2)), by = group][order(group)])
message("WROTE ", out, "/mexicana_races.png, mexicana_races_pools.tsv")
