#!/usr/bin/env Rscript
## =============================================================================
## 34 — Does the donor's home environment explain (A) the teosinte k-mer profile and (B) the QTL effects?
##
## A. Per accession (parviglumis and mexicana separately, 36 and 37 accessions): the mexicana index
##    f_Zx / (f_Zx + f_Zv) and the first two within-taxon PCs of the per-assembly k-mer containment
##    (centred log-ratio, prcomp) against every environment layer of step 33: Spearman rho and p, Benjamini-
##    Hochberg q over all layers x responses, and the partial rho given elevation (Spearman correlation of the
##    rank residuals after regressing both on the rank of elevation), to see whether a layer adds anything
##    beyond elevation.
## B. For every QTL of the joint model (step 29.2, taxon run), whether the effect of the teosinte allele changes
##    with the donor's environment: lmer(y ~ taxon * dosage + env + dosage:env + (1 | accession)) (lme4 +
##    lmerTest), where dosage = teosinte copies at the QTL window from the k-mer calls (B73 0, SEG 1,
##    DONOR 2), env = the donor accession's value (z-score). dosage:taxon is in the model, so dosage:env asks
##    about the environment *within* taxa, not the taxon differences already tested by "taxa differ".
##    BH q over all QTL x layers. A simple mixed model without the LOCO kinship: a screen; hits should be
##    confirmed with R/qtl2 (scan1 with intcovar = env).
## Usage (repo root): Rscript final_script/34_env_tests.R [out dir, default data/env/tests]
## Inputs: data/env/BZea_accession_env.csv (33), data/introgression_spread/sourmash/kmer_taxon_profile.tsv (32b),
##         data/introgression_spread/tracts_all.tsv.gz + pool_summary.tsv (21),
##         data/phenotypes/BZea_spats_BLUEs_combined.csv (28), data/introgression_spread/joint_qtl/noperm (29.2)
## Env:    LAYERS_B (default "elevation,BIO1,BIO6,BIO12,BIO15,aridity,phh2o_0_5cm") — layers tested in part B
## Needs:  R with data.table, ggplot2, lme4, lmerTest.
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(lme4); library(lmerTest) })
out <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "data/env/tests"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
D <- "data/introgression_spread"
env <- fread("data/env/BZea_accession_env.csv")
lab <- c(BIO1 = "annual mean temperature", BIO2 = "mean diurnal range", BIO3 = "isothermality", BIO4 = "temperature seasonality",
         BIO5 = "max temp. warmest month", BIO6 = "min temp. coldest month", BIO7 = "temperature annual range",
         BIO8 = "mean temp. wettest quarter", BIO9 = "mean temp. driest quarter", BIO10 = "mean temp. warmest quarter",
         BIO11 = "mean temp. coldest quarter", BIO12 = "annual precipitation", BIO13 = "precip. wettest month",
         BIO14 = "precip. driest month", BIO15 = "precipitation seasonality", BIO16 = "precip. wettest quarter",
         BIO17 = "precip. driest quarter", BIO18 = "precip. warmest quarter", BIO19 = "precip. coldest quarter",
         elevation = "elevation (metadata)", aridity = "aridity index", et0 = "reference ET", MAP = "mean annual precip.",
         NPP = "net primary productivity", nitrogen_0_5cm = "soil N", phh2o_0_5cm = "soil pH", cec_0_5cm = "soil CEC",
         soc_0_5cm = "soil organic C", ocd_0_5cm = "organic C density", clay_0_5cm = "clay", sand_0_5cm = "sand",
         silt_0_5cm = "silt", bdod_0_5cm = "bulk density", cfvo_0_5cm = "coarse fragments", AM_roots = "AM roots colonized",
         EcM_roots = "EcM roots colonized")
layers <- intersect(names(lab), names(env))
layers <- layers[sapply(layers, function(v) sum(!is.na(env[[v]])) >= 0.8 * nrow(env))]     # drop sparse layers
E <- env[, c("accession-id", layers), with = FALSE]; setnames(E, 1, "accession")

## ---------------- A: k-mer profile vs environment ----------------
k <- fread(file.path(D, "sourmash", "kmer_taxon_profile.tsv"))[role == "donor" & label %in% c("Zv", "Zx")]
cc <- grep("^c_", names(k), value = TRUE)
k[, mx := f_Zx / (f_Zx + f_Zv)]
resA <- list(); acc_all <- list()
for (tx in c("Zv", "Zx")) {
  x <- k[label == tx]
  X <- as.matrix(x[, ..cc]); X <- X / rowSums(X); L <- log(pmax(X, min(X[X > 0]) / 2)); L <- L - rowMeans(L)
  pc <- prcomp(L); x[, `:=`(PC1 = pc$x[, 1], PC2 = pc$x[, 2])]
  a <- x[, .(mexicana_index = median(mx), kmer_PC1 = median(PC1), kmer_PC2 = median(PC2), pools = .N), by = accession]
  a <- merge(a, E, by = "accession"); a[, taxon := tx]; acc_all[[tx]] <- a
  for (resp in c("mexicana_index", "kmer_PC1", "kmer_PC2")) for (v in layers) {
    ok <- complete.cases(a[[resp]], a[[v]], a$elevation); if (sum(ok) < 10) next
    ct <- suppressWarnings(cor.test(a[[resp]][ok], a[[v]][ok], method = "spearman"))
    pr <- NA_real_
    if (v != "elevation") {
      re <- rank(a$elevation[ok])
      pr <- cor(resid(lm(rank(a[[resp]][ok]) ~ re)), resid(lm(rank(a[[v]][ok]) ~ re)))
    }
    resA[[length(resA) + 1]] <- data.table(taxon = tx, response = resp, layer = v, what = lab[v], n = sum(ok),
                                           rho = round(ct$estimate, 3), p = ct$p.value, partial_rho_given_elevation = round(pr, 3),
                                           rho_layer_elevation = round(cor(a[[v]][ok], a$elevation[ok], method = "spearman"), 2))
  }
}
A <- rbindlist(resA); A[, q := p.adjust(p, "BH")]; A <- A[order(taxon, response, p)]
fwrite(A, file.path(out, "A_kmer_profile_vs_env.tsv"), sep = "\t")
fwrite(rbindlist(acc_all), file.path(out, "A_accession_values.tsv"), sep = "\t")
cat("\n== A: k-mer profile vs environment, best layers (q < 0.05) ==\n")
print(A[q < 0.05, .(taxon, response, layer, what, n, rho, q = signif(q, 2), partial_rho_given_elevation, rho_layer_elevation)][order(taxon, response, q)][, head(.SD, 8), by = .(taxon, response)])
## geography: within parviglumis, is the mexicana index a west-east gradient, an elevation gradient, or both?
geo <- fread("data/env/BZea_accession_env.csv")[, .(accession = `accession-id`, latitude, longitude)]
zv <- merge(acc_all[["Zv"]], geo, by = "accession")
prank <- function(y, x, zs) { Z <- cbind(1, sapply(zs, rank)); f <- function(w) resid(lm.fit(Z, rank(w))); cor.test(f(y), f(x)) }
G <- rbindlist(list(
  data.table(test = "index ~ longitude", rho = cor(zv$mexicana_index, zv$longitude, method = "spearman"),
             p = cor.test(zv$mexicana_index, zv$longitude, method = "spearman")$p.value),
  data.table(test = "index ~ elevation", rho = cor(zv$mexicana_index, zv$elevation, method = "spearman"),
             p = cor.test(zv$mexicana_index, zv$elevation, method = "spearman")$p.value),
  { t <- prank(zv$mexicana_index, zv$longitude, list(zv$elevation)); data.table(test = "index ~ longitude | elevation", rho = t$estimate, p = t$p.value) },
  { t <- prank(zv$mexicana_index, zv$elevation, list(zv$longitude)); data.table(test = "index ~ elevation | longitude", rho = t$estimate, p = t$p.value) },
  { t <- prank(zv$mexicana_index, zv$BIO4, list(zv$longitude, zv$elevation)); data.table(test = "index ~ BIO4 | longitude, elevation", rho = t$estimate, p = t$p.value) }))
G[, rho := round(rho, 3)]; fwrite(G, file.path(out, "A_parviglumis_geography.tsv"), sep = "\t")
cat("\n== A: parviglumis mexicana index, geography vs elevation (Spearman; | = partial) ==\n"); print(G)
pa <- ggplot(A[response == "mexicana_index"], aes(factor(layer, levels = layers), taxon, fill = rho)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(q < 0.05, sprintf("%.2f", rho), "")), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", high = "#B2182B", mid = "white", limits = c(-1, 1), name = "Spearman rho") +
  scale_x_discrete(labels = function(z) sprintf("%s (%s)", z, lab[z])) +
  labs(x = NULL, y = NULL, title = "mexicana index of the accession vs its environment (numbers: BH q < 0.05)") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 55, hjust = 1))
ggsave(file.path(out, "A_mexicana_index_vs_env.png"), pa, width = 13, height = 4.2, dpi = 200)

## ---------------- B: QTL effect x donor environment ----------------
LB <- strsplit(Sys.getenv("LAYERS_B", "elevation,BIO1,BIO6,BIO12,BIO15,aridity,phh2o_0_5cm"), ",")[[1]]
LB <- intersect(LB, layers)
ps <- fread(file.path(D, "pool_summary.tsv"))[ROLE == "donor" & FLAG == "ok"]
tr <- fread(file.path(D, "tracts_all.tsv.gz"))[SAMPLE %in% ps$SAMPLE]
ph <- fread("data/phenotypes/BZea_spats_BLUEs_combined.csv")
ph <- merge(ph, ps[, .(Sequencing_ID = SAMPLE, taxon = GROUP)], by = "Sequencing_ID")
ph[, accession := sub("_.*", "", new_genotype)]
ph <- merge(ph, E, by = "accession")
dose_at <- function(chr, pos) {
  w <- tr[CHROM == chr & START <= pos & END > pos, .(Sequencing_ID = SAMPLE, d = c(B73 = 0, SEG = 1, DONOR = 2)[STATE])]
  w[!duplicated(Sequencing_ID)]
}
traits <- c("DTS", "DTA", "ASI", "PH", "EH", "EN", "Prolif", "NBR", "SL", "BL", "BW")
resB <- list()
for (t in traits) {
  qf <- file.path(D, "joint_qtl", "noperm", sprintf("qtl2_%s_qtl.tsv", t)); if (!file.exists(qf)) next
  Q <- fread(qf)
  for (i in seq_len(nrow(Q))) {
    pos <- as.numeric(sub(".*:", "", Q$qtl[i])); dd <- dose_at(Q$chr[i], pos)
    x <- merge(ph[!is.na(get(t))], dd, by = "Sequencing_ID")
    x[, y := get(t)]
    for (v in LB) {
      z <- x[!is.na(get(v))]; if (sum(z$d > 0) < 20) next
      z[, e := as.numeric(scale(get(v)))]
      f <- tryCatch(suppressMessages(lmer(y ~ taxon * d + e + d:e + (1 | accession), data = z)), error = function(err) NULL)
      if (is.null(f)) next
      co <- summary(f)$coefficients
      if (!"d:e" %in% rownames(co)) next
      resB[[length(resB) + 1]] <- data.table(trait = t, qtl = Q$qtl[i], layer = v, what = lab[v], plots = nrow(z),
                                             carriers = sum(z$d > 0), effect_per_copy = round(Q$effect_per_copy_joint[i], 3),
                                             change_per_copy_per_SD_env = round(co["d:e", "Estimate"], 3),
                                             se = round(co["d:e", "Std. Error"], 3), p = co["d:e", "Pr(>|t|)"])
    }
  }
}
B <- rbindlist(resB); B[, q := p.adjust(p, "BH")]; B <- B[order(p)]
fwrite(B, file.path(out, "B_qtl_effect_x_env.tsv"), sep = "\t")
cat(sprintf("\n== B: QTL x environment, %d tests (%d QTL x %d layers); smallest p ==\n", nrow(B), uniqueN(B$qtl), length(LB)))
print(head(B[, .(trait, qtl, layer, carriers, effect_per_copy, change_per_copy_per_SD_env, se, p = signif(p, 2), q = signif(q, 2))], 15))
message("WROTE ", out, "/A_kmer_profile_vs_env.tsv, A_accession_values.tsv, A_mexicana_index_vs_env.png, B_qtl_effect_x_env.tsv")
