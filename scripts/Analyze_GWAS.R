################################################################################
### PACKAGES LOAD
################################################################################
rm(list=ls())

# Load packages
library(vroom)
library(bigmemory)
library(biganalytics)
library(compiler)
library(dplyr)
library(qtl2)
library(GAPIT)

keep <- c("package:stats","package:graphics","package:grDevices","package:utils",
          "package:datasets","package:methods","package:base")
for (pkg in setdiff(search()[grepl("^package:", search())], keep)) {
  detach(pkg, character.only = TRUE, unload = TRUE)
}

library(GAPIT)


################################################################################
### LOAD THE DATA AND PCA
################################################################################

# Read phenotype
myY_FT <- read.csv("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv", header = TRUE)
myY_FT <- myY_FT[,-c(2,4)]
# Remove the PNs
myY_FT <- myY_FT[!grepl("PN", myY_FT$new_genotype), ]
samples_to_keep <- myY_FT$new_genotype
str(myY_FT)

# Read PCA and subset to common individuals
pca <- read.table("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/BZea.beagle.imputed.allchr.renamed.pca.eigenvec") %>% dplyr::select(c(1,3:4))
colnames(pca) <- c("sample.id","EV1","EV2")
pca <- pca[pca$sample.id %in% samples_to_keep, ]
str(pca)

# get samples names and filter for myY_FT
samples_to_keep <- pca$sample.id
myY_FT <- myY_FT[myY_FT$new_genotype %in% samples_to_keep, ]

str(myY_FT)

myY_FT <- as.data.frame(myY_FT, stringsAsFactors = FALSE)
pca    <- as.data.frame(pca,    stringsAsFactors = FALSE)

colnames(myY_FT)[1] <- "taxa"
colnames(pca)[1]    <- "taxa"

# (optional but recommended) make covariate names sane
colnames(pca)[2:3]  <- c("PC1","PC2")

# remove duplicates just in case
myY_FT <- myY_FT[!duplicated(myY_FT$taxa), , drop = FALSE]
pca    <- pca[!duplicated(pca$taxa),    , drop = FALSE]

myY_FT$taxa <- trimws(myY_FT$taxa)
pca$taxa    <- trimws(pca$taxa)

# RUN GWAS
# Loop through chromosomes for GAPIT models
basewd <- "/Users/nirwantandukar/Documents/Github/BZea_genotyping"
folder_hapmap <- "/Users/nirwantandukar/Documents/Research/data/BZea/genotype/"
for (chr in 1:10) {
  out_dir <- file.path(basewd, "results", paste0("chr", chr))
  #dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  hapmap_file <- file.path(folder_hapmap,
                           paste0("BZea.beagle.imputed.allchr.renamed.chr", chr, ".hmp.txt"))
  
  myG_N <- read.delim(
    hapmap_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  
  
  myGAPIT <- GAPIT(
    Y = myY_FT,
    G = myG_N,
    CV = pca,
    model = "FarmCPU",
    file.output = TRUE
    #Geno.View.output = FALSE,   # <<< prevents the layout.matrix crash
    #Phenotype.View = FALSE
  )
  
  saveRDS(myGAPIT, file = paste0("GAPIT_FarmCPU_DTA_chr", chr, ".rds"))
}
setwd(basewd)



# --- enforce GAPIT ID column name ---

# quick overlap check (should be ~897)
cat("Y ∩ CV:", length(intersect(myY_FT$taxa, pca$taxa)), "\n")



# 1. Check sample counts
cat("Phenotype samples:", nrow(myY_FT), "\n")
cat("PCA samples:", nrow(pca), "\n")

# 2. Check for exact matches
pheno_samples <- myY_FT[,1]
pca_samples <- pca[,1]

# Read a small chunk of genotype to check samples
genotype_samples <- colnames(myG_N)

cat("Sample match check:\n")
cat("Pheno-PCA match:", sum(pheno_samples %in% pca_samples), "\n")
cat("Pheno-Geno match:", sum(pheno_samples %in% genotype_samples), "\n")
cat("PCA-Geno match:", sum(pca_samples %in% genotype_samples), "\n")




myY_N
pca
myG_N

