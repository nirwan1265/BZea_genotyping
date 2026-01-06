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

post_post.processing <- TRUE


in_dir  <- "/Users/nirwantandukar/Documents/Research/data/BZea/BZea_allele_counts_25k/"
files   <- list.files(in_dir, pattern="\\.tsv$", full.names=TRUE)
files <- files[300:325]

expDesign <- data.frame(
  files = files,
  name  = sub("\\.tsv$", "", basename(files)),
  stringsAsFactors = FALSE
)

outdir <- "/Users/nirwantandukar/Documents/Research/data/BZea/rtiger_results_121_200"

myres <- RTIGER(
  expDesign = expDesign,
  outputdir = outdir,
  seqlengths = chr_len,
  
  # START CONSERVATIVE for your pruned marker density
  rigidity = 20,
  # Best rigidity = 512
  
  # Let RTIGER tune (but give realistic depth)
  autotune = FALSE,
  average_coverage = 0.8,
  crossovers_per_megabase = 0.05,
  
  nstates = 3,
  post.processing = post_post.processing,
  save.results = TRUE,
  verbose = TRUE
)

saveRDS(myres,"Samples_121_200_RTIGER_results.rds")

str(myres@Probabilities$gamma)
str(myres)
myres@Viterbi

myres@matobs
myres@params
