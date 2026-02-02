# Packages
#BiocManager::install(version = "3.14")
#BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "TailRank", "IRanges", "Gviz"))
#install.packages("JuliaCall")
library(JuliaCall)
library(dplyr)
#install.packages("RTIGER")
library(RTIGER)


### SETUP
# Done once
#setupJulia(JULIA_HOME="/Applications/Julia-1.0.app/Contents/Resources/julia/bin")
setupJulia(JULIA_HOME="/Applications/Julia-1.10.app/Contents/Resources/julia/bin")
# Needs to be run everytime we load RTIGER
sourceJulia()


library(doParallel)
library(foreach)
library(doParallel)
library(foreach)
library(tools)
library(RTIGER)


# ----- RTIGER constants -----
chr_len <- c(308452471,243675191,238017767,250330460,226353449,
             181357234,185808916,182411202,163004744,152435371)
names(chr_len) <- paste0("chr", 1:10)


chr_len <- c(308452471)
names(chr_len) <- paste0("chr", c(1))

post_post.processing <- TRUE


in_dir  <- "/Users/nirwantandukar/Documents/Research/data/BZea/graphtyper/allele_counts/"
files   <- list.files(in_dir, pattern="\\.tsv$", full.names=TRUE)
files <- files[c(1)]

expDesign <- data.frame(
  files = files,
  name  = sub("\\.rtiger.chr1-10.ref_refc_alt_altc.DP2plus.tsv$", "", basename(files)),
  stringsAsFactors = FALSE
)

outdir <- "/Users/nirwantandukar/Documents/Research/data/BZea/rtiger_results/rtiger_results_graphtyper"

myres <- RTIGER(
  expDesign = expDesign,
  outputdir = outdir,
  seqlengths = chr_len,
  
  # START CONSERVATIVE for your pruned marker density
  rigidity = 40,
  # Best rigidity = 512
  
  # Let RTIGER tune (but give realistic depth)
  autotune = F,
  average_coverage = 0.08,
  #crossovers_per_megabase = 0.05,
  crossovers_per_megabase = 0.01,
  
  nstates = 3,
  post.processing = post_post.processing,
  save.results = TRUE,
  verbose = TRUE
)

optimize_R(myres, average_coverage = 0.8)

# [1] 323

# after you have a myres
co_tbl <- calcCOnumber(myres)  # CO counts per sample/chrom (RTIGER helper)

# genome size in Mb
genome_mb <- sum(as.numeric(chr_len)) / 1e6

# per-sample CO/Mb (rough)
co_per_mb <- rowSums(co_tbl, na.rm = TRUE) / genome_mb

summary(co_per_mb)
quantile(co_per_mb, c(0.9, 0.95, 0.99), na.rm = TRUE)





saveRDS(myres,"Samples_chr1_2_201_300_RTIGER_results.rds")

str(myres@Probabilities$gamma)
str(myres)
myres@Viterbi

myres@matobs
myres@params
