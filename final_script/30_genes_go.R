#!/usr/bin/env Rscript
## =============================================================================
## 30 — Genes in the significant regions (step 27) joined to their GO annotation, with trait-relevant
##      genes flagged by keyword.
##
## Usage:  Rscript 30_genes_go.R <genes dir (step 27 output)> <GO file> <out dir>
##   genes dir: regions_genes_<TRAIT>_ALL.tsv (preferred: has the region's -log10 p and effect) or peak_genes_<TRAIT>_ALL.tsv
##   GO file:   tab-separated, columns GeneID GeneName Family_Subfamily ProteinClass GO_MF GO_BP GO_CC
## Outputs (<out dir>):
##   genes_GO_all.tsv            every gene in every region, with GO (one row per trait x gene)
##   genes_GO_flagged.tsv        genes whose name / family / GO matches the keywords for that trait
##   genes_GO_flagged_<TRAIT>.tsv
## Needs: R with dplyr, readr, stringr, purrr (tidyverse).
## Example (laptop):
##   Rscript final_script/30_genes_go.R ~/Documents/Research/data/BZea/introgression_spread/assoc_gridlmm_spats/genes \
##     ~/Documents/Research/data/BZea/introgression_spread/GO/GO_terms_all_candidates.txt \
##     ~/Documents/Research/data/BZea/introgression_spread/GO
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(purrr) })
args <- commandArgs(TRUE)
if (length(args) < 3) stop("usage: Rscript 30_genes_go.R <genes dir> <GO file> <out dir>")
genes_dir <- args[1]; go_file <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- keywords per trait group (regular expressions, case-insensitive) ----
kw <- list(
  flowering = "flower|heading|photoperiod|circadian|CONSTANS|CCT|GHD7|PEBP|ZCN|FLOWERING LOCUS|TERMINAL FLOWER|MADS|phytochrome|GIGANTEA|vernaliz|EARLY FLOWERING|LATE ELONGATED|pseudo-response|PRR\\d|TIMING OF CAB|SQUAMOSA|short-day|long-day|ADAGIO|LNK\\d|Rap2\\.7|HD3A|tga1",
  height    = "gibberellin|brassinosteroid|DELLA|GRAS|internode|stem elongation|auxin efflux|PIN-FORMED|cell elongation|expansin|dwarf|GROWTH-REGULATING|GRF\\d|ent-kaurene|BZR|BSK\\d|CYP734|GA2|GA20|GA3|THICK TASSEL|ramosa|ra2|CCT|GHD7",
  branching = "branch|tiller|strigolactone|axillary|meristem|TEOSINTE BRANCHED|TCP|cytokinin|CKX|KNOTTED|inflorescence|spikelet|BARREN|RAMOSA|ra2|shoot apical|carotenoid cleavage|CCD\\d|D14|MAX2|LAX PANICLE|FEA\\d|PERIANTHIA|liguleless|LG\\d|lg4|ZCN|PEBP|SQUAMOSA|tsh|Rap2\\.7"
)
group_of <- c(DTS = "flowering", DTA = "flowering", ASI = "flowering", PH = "height", EH = "height",
              EN = "branching", Prolif = "branching", NBR = "branching", SL = "branching", BL = "branching",
              BW = "branching", LAE = "branching", Kinki = "branching", SPAD = "flowering", SPAD2 = "flowering",
              StPu = "branching", StPi = "branching")

## ---- genes in the regions ----
files <- list.files(genes_dir, pattern = "^regions_genes_.*_ALL\\.tsv$", full.names = TRUE)
use_regions <- length(files) > 0
if (!use_regions) files <- list.files(genes_dir, pattern = "^peak_genes_.*_ALL\\.tsv$", full.names = TRUE)
if (!length(files)) stop("no regions_genes_*_ALL.tsv or peak_genes_*_ALL.tsv in ", genes_dir)
genes <- map_dfr(files, function(f) {
  trait <- str_match(basename(f), "_genes_(.+)_ALL\\.tsv$")[, 2]
  d <- read_tsv(f, col_types = cols(.default = "c"))
  if (use_regions) {
    d %>% filter(GENE != "") %>%
      transmute(TRAIT = trait, REGION = as.integer(REGION), CHR, INTERVAL_START = as.numeric(INTERVAL_START),
                INTERVAL_END = as.numeric(INTERVAL_END), TOP_WINDOW, TOP_MLOG10P = as.numeric(TOP_MLOG10P),
                BETA_PER_COPY = as.numeric(BETA_PER_COPY), GENE, GENE_START = as.numeric(GENE_START),
                GENE_END = as.numeric(GENE_END), NAME, DESCRIPTION)
  } else {
    d %>% transmute(TRAIT = trait, REGION = as.integer(PEAK), CHR, GENE, GENE_START = as.numeric(START),
                    GENE_END = as.numeric(END), NAME, DESCRIPTION)
  }
})

## ---- join GO ----
go <- read_tsv(go_file, col_types = cols(.default = "c")) %>% distinct(GeneID, .keep_all = TRUE)
all <- genes %>% left_join(go, by = c("GENE" = "GeneID")) %>%
  mutate(across(c(NAME, DESCRIPTION, GeneName, Family_Subfamily, ProteinClass, GO_MF, GO_BP, GO_CC), ~ coalesce(.x, "")),
         group = unname(group_of[TRAIT]),
         text = paste(NAME, DESCRIPTION, GeneName, Family_Subfamily, GO_BP, GO_MF, sep = " | "))
all$flag <- map2_lgl(all$text, all$group, ~ !is.na(.y) && str_detect(.x, regex(kw[[.y]], ignore_case = TRUE)))
all$keywords <- map2_chr(all$text, all$group, function(t, g) {
  if (is.na(g)) return("")
  hits <- str_extract_all(t, regex(kw[[g]], ignore_case = TRUE))[[1]]
  paste(unique(tolower(hits)), collapse = ";")
})
all <- all %>% select(-text) %>% arrange(TRAIT, REGION, GENE_START)
write_tsv(all, file.path(out_dir, "genes_GO_all.tsv"), na = "")
fl <- all %>% filter(flag)
write_tsv(fl, file.path(out_dir, "genes_GO_flagged.tsv"), na = "")
fl %>% group_split(TRAIT) %>% walk(~ write_tsv(.x, file.path(out_dir, sprintf("genes_GO_flagged_%s.tsv", .x$TRAIT[1])), na = ""))

## ---- summary ----
summ <- all %>% group_by(TRAIT) %>%
  summarise(regions = n_distinct(REGION), genes = n_distinct(GENE), with_GO = n_distinct(GENE[GeneName != ""]),
            flagged = n_distinct(GENE[flag]), .groups = "drop")
print(as.data.frame(summ), row.names = FALSE)
message("WROTE ", out_dir, "/genes_GO_all.tsv, genes_GO_flagged.tsv, genes_GO_flagged_<TRAIT>.tsv")
