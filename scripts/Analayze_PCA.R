library(tidyverse)
library(plotly)
library(htmlwidgets)
library(ggplot2)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 24,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 24,      # X‐axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 24,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),

    # ---- grids: light grey major grids for x + y, no minor grids ----
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),

    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 16),

    plot.margin    = margin(15, 15, 15, 15)

  )

read_plink_eigenvec <- function(path) {
  df <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  
  n_pcs <- ncol(df) - 2
  stopifnot(n_pcs >= 2)
  
  colnames(df) <- c("FID", "IID", paste0("PC", seq_len(n_pcs)))
  
  df %>%
    mutate(
      sample = IID,
      Group  = case_when(
        str_detect(sample, "^PN") ~ NA_character_,              # drop PN*
        str_detect(sample, "\\.") ~ str_remove(sample, "\\..*"),# prefix before first "."
        TRUE ~ sample                                           # fallback (e.g., B73)
      )
      # If you REALLY want the dot in group labels (Zv., Zl.), use:
      # Group = if_else(!is.na(Group) & !str_ends(Group, "\\."), paste0(Group, "."), Group)
    ) %>%
    filter(!is.na(Group))
}

plot_pca <- function(df, title = NULL, level = 0.95) {
  ggplot(df, aes(PC2, PC3, color = Group)) +
    # 95% confidence ellipse per group
    stat_ellipse(
      aes(group = Group),
      type  = "t",
      level = level,
      linewidth = 0.8
    ) +
    geom_point(alpha = 0.75, size = 1.6) +
    theme_bw() +
    labs(title = title, x = "PC2", y = "PC3") +
    theme(
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    ) +
    plot_theme
}

# ---- Load ----
filtered_imputed_pca <- read_plink_eigenvec(
  "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/BZea.beagle.imputed.allchr.renamed.pca.eigenvec"
)

filtered_unimputed_pca <- read_plink_eigenvec(
  "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/50perc_missing/BZea.DP2.MAF005.MISS50.allchr.renamed.pca.eigenvec"
)



# ---- Plot ----
p_imp  <- plot_pca(filtered_imputed_pca, "Imputed (unfiltered) PCA: PC1 vs PC2")
p_unim <- plot_pca(filtered_unimputed_pca, "Unimputed (filtered) PCA: PC1 vs PC2")

quartz()
p_imp
quartz()
p_unim

# ---- Optional: save ----
ggsave("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/PCA_filtered_imputed_PC1_PC2.png",  p_imp,  width = 7, height = 5, dpi = 300, bg = "white")
ggsave("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/PCA_filtered_unimputed_PC1_PC2.png", p_unim, width = 7, height = 5, dpi = 300, bg = "white")




### 3D pca
plot_pca_3d <- function(df,
                        title = "PCA: PC1 vs PC2 vs PC3",
                        point_size = 3,
                        show_labels = FALSE) {
  df <- df %>%
    dplyr::filter(!is.na(Group)) %>%
    dplyr::mutate(
      Group  = as.factor(Group),
      label  = if ("sample" %in% names(df)) df$sample else df$IID
    )
  
  p <- plot_ly(
    data   = df,
    x      = ~PC1,
    y      = ~PC2,
    z      = ~PC3,
    color  = ~Group,
    colors = "Set2",
    type   = "scatter3d",
    mode   = if (show_labels) "markers+text" else "markers",
    marker = list(size = point_size),
    text   = if (show_labels) ~label else NULL,
    hoverinfo = "text",
    hovertext = ~paste0(
      "Sample: ", label,
      "<br>Group: ", Group,
      "<br>PC1: ", signif(PC1, 4),
      "<br>PC2: ", signif(PC2, 4),
      "<br>PC3: ", signif(PC3, 4)
    )
  ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      )
    )
  
  p
}

# ---- examples ----
p3d_imp  <- plot_pca_3d(unfiltered_imputed_pca,  title = "Imputed (unfiltered): PC1–PC3")
p3d_unim <- plot_pca_3d(filtered_unimputed_pca, title = "Unimputed (filtered): PC1–PC3")

quartz()
p3d_imp
p3d_unim

# ---- Save interactive HTML ----
saveWidget(p3d_imp,  "data/pca/PCA3D_imputed_PC1_PC2_PC3.html",  selfcontained = TRUE)
saveWidget(p3d_unim, "data/pca/PCA3D_unimputed_PC1_PC2_PC3.html", selfcontained = TRUE)




# Scree Plot

# Variance explained (%) and plot a scree.
plot_scree_plink <- function(eigenval_file,
                             title = "Scree plot",
                             n_show = 20) {
  ev <- scan(eigenval_file, quiet = TRUE)
  stopifnot(length(ev) > 0)
  
  df <- tibble(
    PC = seq_along(ev),
    eigenvalue = ev,
    var_frac = ev / sum(ev),
    var_pct  = 100 * var_frac,
    cum_pct  = 100 * cumsum(var_frac)
  ) %>%
    slice_head(n = min(n_show, nrow(.)))
  
  ggplot(df, aes(x = PC, y = var_pct)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(
      title = title,
      x = "Principal component",
      y = "Variance explained (%)"
    ) + plot_theme
}

ev <- scan("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/BZea.beagle.imputed.allchr.renamed.pca.eigenval")

# ---- examples ----
scree_imp <- plot_scree_plink(
  "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/BZea.beagle.imputed.allchr.renamed.pca.eigenval",
  title = "Imputed (unfiltered): Scree plot",
  n_show = 15
)

scree_unim <- plot_scree_plink(
  "/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/BZea.DP2.MAF005.MISS50.allchr.renamed.pca.eigenval",
  title = "Unimputed (filtered): Scree plot",
  n_show = 15
)

quartz()
scree_imp
quartz()
scree_unim

# ---- optional save ----
ggsave("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/Scree_imputed.png",  scree_imp,  width = 7, height = 5, dpi = 300)
ggsave("/Users/nirwantandukar/Documents/Github/BZea_genotyping/data/pca/Scree_unimputed.png", scree_unim, width = 7, height = 5, dpi = 300)

