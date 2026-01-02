# -----------------------------
# BZea SpATS: clean build + namechange-aware join + 2-block fitting
# -----------------------------
library(dplyr)
library(tidyr)
library(readr)
library(SpATS)

# -----------------------------
# 0) Read inputs
# -----------------------------
fieldmap  <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/fieldmap.csv",
                      header = FALSE, stringsAsFactors = FALSE)

fieldID   <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/fieldID.csv",
                      header = TRUE, stringsAsFactors = FALSE)

phenotype <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/field_phenotype_filtered.csv",
                      header = TRUE, stringsAsFactors = FALSE)

borders   <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/field_borders.csv",
                      header = TRUE, stringsAsFactors = FALSE)

block1    <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/block1.csv",
                      header = FALSE, stringsAsFactors = FALSE)

irri      <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/irrigation.csv",
                      header = FALSE, stringsAsFactors = FALSE)

namechange <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/name_change.csv",
                       header = TRUE, stringsAsFactors = FALSE)

namechange_new <- read.csv("/Users/nirwantandukar/Documents/Research/data/BZea/name_change_new.csv",
                           header = TRUE, stringsAsFactors = FALSE)

str(namechange)
str(namechange_new)

namechange2 <- namechange %>%
  mutate(Female_genotype = as.character(Female_genotype)) %>%
  left_join(
    namechange_new %>%
      mutate(oldold_genotype = as.character(oldold_genotype)) %>%
      select(oldold_genotype, new_genotype),
    by = c("Female_genotype" = "oldold_genotype")
  ) %>%
  select(Female_genotype, FieldID, FieldID_Species, Sequencing_ID, new_genotype)

# quick sanity checks
cat("Rows in namechange2:", nrow(namechange2), "\n")
cat("Missing new_genotype:", sum(is.na(namechange2$new_genotype) | namechange2$new_genotype == ""), "\n")
cat("Missing Sequencing_ID:", sum(is.na(namechange2$Sequencing_ID) | namechange2$Sequencing_ID == ""), "\n")

# optional: see which ones failed to map
namechange2 %>%
  filter(is.na(new_genotype) | new_genotype == "") %>%
  select(Female_genotype, FieldID_Species, Sequencing_ID) %>%
  head(20)

length(namechange2$Sequencing_ID)

# -----------------------------
# 1) Layout: fieldmap matrix -> (FieldID, row, col)
#    Remove "X" holes.
# -----------------------------
layout_df <- fieldmap %>%
  mutate(row = row_number()) %>%
  pivot_longer(cols = -row, names_to = "col_name", values_to = "FieldID") %>%
  mutate(
    col = readr::parse_number(col_name),
    FieldID = as.character(FieldID)
  ) %>%
  filter(!is.na(FieldID), FieldID != "X") %>%
  select(FieldID, row, col)

# -----------------------------
# 2) Plotbook: FieldID, Rep, LineID + Block1 membership from block1.csv
# -----------------------------
block1_ids <- as.character(as.vector(as.matrix(block1)))

plotbook <- fieldID %>%
  mutate(
    FieldID = as.character(FieldID),
    Rep     = factor(Rep),
    LineID  = as.character(LineID)
  ) %>%
  mutate(
    Block1 = factor(ifelse(FieldID %in% block1_ids, 0, 1))  # 0=in block1_ids, 1=other field
  ) %>%
  select(FieldID, Rep, LineID, Block1)

# -----------------------------
# 3) Add border + irrigation flags at plot level
# -----------------------------
borders_vec <- as.character(borders$Borders)
irri_vec    <- as.character(irri[[1]])

plot_flags <- layout_df %>%
  mutate(
    Border = factor(ifelse(FieldID %in% borders_vec, 1, 0)),
    Irri   = factor(ifelse(FieldID %in% irri_vec,    1, 0))
  ) %>%
  select(FieldID, row, col, Border, Irri)

# -----------------------------
# 4) Merge layout + plotbook
# -----------------------------
dat0 <- plot_flags %>%
  left_join(plotbook, by = "FieldID")

# sanity: how many FieldIDs in map that have no plotbook record?
cat("Missing plotbook rows:", sum(is.na(dat0$LineID) | is.na(dat0$Rep)), "\n")

# -----------------------------
# 5) Phenotypes: wide (DTS_1..3, DTA_1..3) -> long (Rep1..Rep3)
#    This is the correct pivot (no cross-product).
# -----------------------------
phen_long <- phenotype %>%
  mutate(Female_genotype = as.character(Female_genotype)) %>%
  pivot_longer(
    cols = matches("^(DTS|DTA)_[0-9]+$"),
    names_to = c(".value", "RepNum"),
    names_pattern = "^(DTS|DTA)_([0-9]+)$"
  ) %>%
  mutate(Rep = factor(paste0("Rep", RepNum))) %>%
  select(Female_genotype, Rep, DTS, DTA)

# -----------------------------
# 6) NAMECHANGE: harmonize phenotype genotype IDs to match plotbook LineID
#    Your namechange columns are genotype-name variants (NOT plot FieldID).
#    We automatically choose which column best matches plotbook LineID and map everything to that.
# -----------------------------
plot_ids <- unique(plotbook$LineID)

candidate_cols <- c("Female_genotype", "FieldID_Species", "FieldID")
candidate_cols <- candidate_cols[candidate_cols %in% names(namechange)]

# choose the namechange column with the best overlap to plotbook LineID
overlaps <- sapply(candidate_cols, function(cc) sum(namechange[[cc]] %in% plot_ids, na.rm = TRUE))
target_col <- candidate_cols[which.max(overlaps)]

cat("Namechange target column used to match plotbook LineID:", target_col, "\n")
cat("Overlap counts:", paste(candidate_cols, overlaps, sep="=", collapse="; "), "\n")

# build a single mapping table: any variant -> target_col
maps <- bind_rows(
  if ("Female_genotype" %in% names(namechange))
    transmute(namechange, from = Female_genotype, to = .data[[target_col]]),
  if ("FieldID_Species" %in% names(namechange))
    transmute(namechange, from = FieldID_Species, to = .data[[target_col]]),
  if ("FieldID" %in% names(namechange))
    transmute(namechange, from = FieldID, to = .data[[target_col]])
) %>%
  filter(!is.na(from), from != "", !is.na(to), to != "") %>%
  distinct(from, .keep_all = TRUE)

recode_to_plotbook <- function(x) {
  x <- as.character(x)
  # if already matches plotbook ids, keep
  out <- x
  need <- !(out %in% plot_ids)
  if (any(need)) {
    hit <- match(out[need], maps$from)
    repl <- maps$to[hit]
    out[need & !is.na(hit)] <- repl[!is.na(hit)]
  }
  out
}

phen_long <- phen_long %>%
  mutate(LineID_std = recode_to_plotbook(Female_genotype))

# report mapping success
cat("Phenotypes matching plotbook BEFORE mapping:",
    sum(phen_long$Female_genotype %in% plot_ids), "of", nrow(phen_long), "\n")

cat("Phenotypes matching plotbook AFTER mapping:",
    sum(phen_long$LineID_std %in% plot_ids), "of", nrow(phen_long), "\n")

# show top unmapped names (if any)
unmapped <- phen_long %>%
  filter(!(LineID_std %in% plot_ids)) %>%
  count(Female_genotype, sort = TRUE)

if (nrow(unmapped) > 0) {
  cat("Top unmapped phenotype IDs (need fix in namechange or phenotype file):\n")
  print(head(unmapped, 20))
}

# collapse to ONE value per (LineID_std, Rep)
phen_keyed <- phen_long %>%
  group_by(LineID_std, Rep) %>%
  summarise(
    DTS = if (all(is.na(DTS))) NA_real_ else mean(DTS, na.rm = TRUE),
    DTA = if (all(is.na(DTA))) NA_real_ else mean(DTA, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 7) Join phenotypes onto plots by (LineID, Rep)
#    This is ONE-TO-MANY (good): one genotype mean per rep -> many plots with same genotype.
# -----------------------------
dat1 <- dat0 %>%
  left_join(phen_keyed, by = c("LineID" = "LineID_std", "Rep" = "Rep"))

# -----------------------------
# 8) Irrigation handling (choose ONE)
#    A) If irrigation plots are compromised: set traits to NA (recommended)
#    B) If not compromised: comment this out and just keep Irri as a covariate later
# -----------------------------
dat1 <- dat1 %>%
  mutate(
    DTS = ifelse(Irri == "1", NA, DTS),
    DTA = ifelse(Irri == "1", NA, DTA)
  )

# -----------------------------
# 9) Final SpATS-ready data + integrity checks
# -----------------------------
dat_fit <- dat1 %>%
  filter(!is.na(Rep), !is.na(LineID)) %>%
  mutate(
    row = as.numeric(row),
    col = as.numeric(col),
    RowF = factor(row),
    ColF = factor(col),
    LineID = factor(LineID),
    Rep    = factor(Rep),
    Block1 = factor(Block1),
    Border = factor(Border),
    Irri   = factor(Irri),
    FieldID = factor(FieldID)
  )

# MUST be one row per plot FieldID
if (any(duplicated(dat_fit$FieldID))) stop("FATAL: duplicated FieldID rows. Your joins made >1 row per plot.")

cat("Plots (rows) in dat_fit:", nrow(dat_fit), "\n")
cat("Non-missing DTS:", sum(!is.na(dat_fit$DTS)), "\n")
cat("Non-missing DTA:", sum(!is.na(dat_fit$DTA)), "\n")

# -----------------------------
# 10) Fit SpATS separately for each physical block (recommended with two separate fields)
# -----------------------------
pick_nseg <- function(df) {
  nx <- length(unique(df$col))
  ny <- length(unique(df$row))
  c(max(8, min(20, round(nx / 3))),
    max(8, min(20, round(ny / 3))))
}

fit_spats_block <- function(df, response) {
  seg <- pick_nseg(df)
  SpATS(
    response = response,
    spatial  = ~ SAP(col, row, nseg = seg, degree = 3, pord = 2),
    genotype = "LineID",
    genotype.as.random = FALSE, # TRUE = BLUP; FALSE = BLUE
    fixed   = ~ Rep + Border,     # add + Irri if you did NOT set irrigation traits to NA
    random  = ~ RowF + ColF,      # optional; helps soak up row/col noise beyond smooth
    data    = df,
    control = list(tolerance = 1e-3, monitoring = 0)
  )
}

dat_b0 <- dat_fit %>% filter(Block1 == "0")
dat_b1 <- dat_fit %>% filter(Block1 == "1")

m_DTS_b0 <- fit_spats_block(dat_b0 %>% filter(!is.na(DTS)), "DTS")
m_DTS_b1 <- fit_spats_block(dat_b1 %>% filter(!is.na(DTS)), "DTS")

m_DTA_b0 <- fit_spats_block(dat_b0 %>% filter(!is.na(DTA)), "DTA")
m_DTA_b1 <- fit_spats_block(dat_b1 %>% filter(!is.na(DTA)), "DTA")

# -----------------------------
# 11) Extract genotype predictions (BLUP-like) per block
# -----------------------------
pred_DTS_b0_BLUE <- predict(m_DTS_b0, which = "LineID", predFixed = "marginal")
pred_DTS_b1_BLUE <- predict(m_DTS_b1, which = "LineID", predFixed = "marginal")

pred_DTA_b0_BLUE <- predict(m_DTA_b0, which = "LineID", predFixed = "marginal")
pred_DTA_b1_BLUE <- predict(m_DTA_b1, which = "LineID", predFixed = "marginal")


str(pred_DTS_b0)
str(pred_DTS_b1)
str(pred_DTA_b0)
str(pred_DTA_b1)

# -----------------------------
# 12) Extract genotype predictions (BLU3-like) per block
# -----------------------------
DTS0_BLUE <- pred_DTS_b0_BLUE %>%
  transmute(LineID = as.character(LineID),
            DTS_b0_BLUE = predicted.values,
            SE_DTS_b0_BLUE = standard.errors)

DTS1_BLUE <- pred_DTS_b1_BLUE %>%
  transmute(LineID = as.character(LineID),
            DTS_b1_BLUE = predicted.values,
            SE_DTS_b1_BLUE = standard.errors)

DTS_BLUE <- DTS0_BLUE %>%
  full_join(DTS1_BLUE, by = "LineID") %>%
  mutate(DTS_BLUE = rowMeans(cbind(DTS_b0_BLUE, DTS_b1_BLUE), na.rm = TRUE))



DTA0_BLUE <- pred_DTA_b0_BLUE %>%
  transmute(LineID = as.character(LineID),
            DTA_b0_BLUE = predicted.values,
            SE_DTA_b0_BLUE = standard.errors)

DTA1_BLUE <- pred_DTA_b1_BLUE %>%
  transmute(LineID = as.character(LineID),
            DTA_b1_BLUE = predicted.values,
            SE_DTA_b1_BLUE = standard.errors)

DTA_BLUE <- DTA0_BLUE %>%
  full_join(DTA1_BLUE, by = "LineID") %>%
  mutate(DTA_BLUE = rowMeans(cbind(DTA_b0_BLUE, DTA_b1_BLUE), na.rm = TRUE))


head(DTS_BLUE)
head(DTA_BLUE)

# -----------------------------
# 13) Weighted BLUES across blocks
# -----------------------------

wmean2 <- function(x1, se1, x2, se2) {
  w1 <- ifelse(is.na(x1) | is.na(se1) | se1 <= 0, NA_real_, 1 / se1^2)
  w2 <- ifelse(is.na(x2) | is.na(se2) | se2 <= 0, NA_real_, 1 / se2^2)
  num <- x1*w1 + x2*w2
  den <- w1 + w2
  out <- num / den
  out[is.nan(out)] <- NA_real_
  out
}

DTS_BLUE <- DTS_BLUE %>%
  mutate(
    DTS_BLUE_w = wmean2(DTS_b0_BLUE, SE_DTS_b0_BLUE, DTS_b1_BLUE, SE_DTS_b1_BLUE)
  )

DTA_BLUE <- DTA_BLUE %>%
  mutate(
    DTA_BLUE_w = wmean2(DTA_b0_BLUE, SE_DTA_b0_BLUE, DTA_b1_BLUE, SE_DTA_b1_BLUE)
  )

head(DTS_BLUE)
head(DTA_BLUE)



# -----------------------------
# 14) Export for GWAS
# -----------------------------

# Combine the tables
pheno_gwas <- full_join(
  DTS_BLUE %>% select(LineID, DTS = DTS_BLUE_w),
  DTA_BLUE %>% select(LineID, DTA = DTA_BLUE_w),
  by = "LineID"
)
pheno_gwas <- pheno_gwas[complete.cases(pheno_gwas), ]

# Add the sequencing ID
# --- 0) Ensure Sequencing_ID exists in pheno_gwas (even if empty) ---
if (!"Sequencing_ID" %in% names(pheno_gwas)) {
  pheno_gwas <- pheno_gwas %>% mutate(Sequencing_ID = NA_character_)
}

# --- 1) Pick which namechange column matches your current LineID space ---
plot_ids <- unique(pheno_gwas$LineID)

candidate_cols <- c("Female_genotype", "FieldID_Species", "FieldID")
candidate_cols <- candidate_cols[candidate_cols %in% names(namechange)]

overlaps <- sapply(candidate_cols, function(cc) sum(namechange[[cc]] %in% plot_ids, na.rm = TRUE))
target_col <- candidate_cols[which.max(overlaps)]

cat("Mapping pheno_gwas$LineID using namechange column:", target_col, "\n")
cat("Overlap counts:", paste(candidate_cols, overlaps, sep="=", collapse="; "), "\n")

# --- 2) Build mapping: LineID -> oldold_genotype + Sequencing_ID (from namechange) ---
id_map <- namechange %>%
  transmute(
    LineID = .data[[target_col]],
    oldold_genotype = Female_genotype,
    Sequencing_ID = Sequencing_ID
  ) %>%
  filter(!is.na(LineID), LineID != "") %>%
  distinct(LineID, .keep_all = TRUE)

# --- 3) Build mapping: oldold_genotype -> new_genotype ---
new_map <- namechange_new %>%
  transmute(oldold_genotype = oldold_genotype,
            new_genotype    = new_genotype) %>%
  filter(!is.na(oldold_genotype), oldold_genotype != "",
         !is.na(new_genotype), new_genotype != "") %>%
  distinct(oldold_genotype, .keep_all = TRUE)

# --- 4) Join + resolve Sequencing_ID cleanly + set new_genotype ---
pheno_gwas_new <- pheno_gwas %>%
  left_join(id_map,  by = "LineID", suffix = c("", ".nc")) %>%   # makes Sequencing_ID.nc if needed
  mutate(
    Sequencing_ID = dplyr::coalesce(Sequencing_ID, Sequencing_ID.nc)
  ) %>%
  select(-Sequencing_ID.nc) %>%
  left_join(new_map, by = "oldold_genotype") %>%
  mutate(
    new_genotype = ifelse(is.na(new_genotype) | new_genotype == "", LineID, new_genotype)
  ) %>%
  select(-LineID, -oldold_genotype) %>%
  relocate(new_genotype, .before = 1) %>%
  relocate(Sequencing_ID, .after = last_col())

pheno_gwas_new <- pheno_gwas_new[complete.cases(pheno_gwas_new),]

# --- 5) Sanity checks ---
cat("Missing new_genotype:", sum(is.na(pheno_gwas_new$new_genotype)), "of", nrow(pheno_gwas_new), "\n")
cat("Missing Sequencing_ID:", sum(is.na(pheno_gwas_new$Sequencing_ID)), "of", nrow(pheno_gwas_new), "\n")
cat("Duplicate new_genotype:", any(duplicated(pheno_gwas_new$new_genotype)), "\n")

# --- 6) Write ---
write.csv(pheno_gwas_new,
          "data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv",
          row.names = FALSE)


