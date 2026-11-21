#!/usr/bin/env Rscript
## =============================================================================
## 31 — Side-by-side figure per trait: (A) k-mer additive scan (step 26, GridLMM) and (B) the joint QTL model
##      with per-taxon effects (step 29 / 29.2, R/qtl2) at the same QTL.
##
## Usage:  Rscript 31_joint_qtl_figures.R <scan dir> <joint QTL dir> <out dir> [TRAIT ...]
##   scan dir:      scan_<TRAIT>.tsv and summary_<TRAIT>.tsv (step 26; e.g. introgression_spread/assoc_gridlmm_spats)
##   joint QTL dir: qtl2_<TRAIT>_qtl.tsv and qtl2_<TRAIT>_by_taxon.tsv (step 29 / 29.2; e.g. introgression_spread/joint_QTL)
##   TRAITs:        default = every trait with a qtl2_<TRAIT>_qtl.tsv
## Panel A: -log10 p per 100-kb window, chromosomes in black/grey, dotted line = 5% permutation threshold (when the
##          step-26 summary has one), every QTL of the scan numbered at its top window: filled triangle = it holds in
##          the joint model (joint -log10 p above the Bonferroni cut over the trait's QTL, -log10(0.05 / number of QTL)),
##          open triangle = it does not. "*" after the number = added by the conditional scans of step 29 (a second
##          QTL inside a significant region).
## Panel B: per QTL, effect per teosinte copy (trait units) with 95% CI: black = joint model (all QTL together, all taxa),
##          colours = within each taxon (taxa with >= 5 carrier plots). Right: -log10 p of the joint test and of
##          "taxa differ" (two columns; these are p-values, not effects).
## Outputs: joint_qtl_<TRAIT>.png (one per trait) and joint_qtl_all_traits.pdf (one page per trait).
## Needs:   R with ggplot2, patchwork, data.table.
## =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table) })
args <- commandArgs(TRUE)
if (length(args) < 3) stop("usage: Rscript 31_joint_qtl_figures.R <scan dir> <joint QTL dir> <out dir> [TRAIT ...]")
scan_dir <- args[1]; jq_dir <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
traits <- if (length(args) > 3) args[-(1:3)] else
  sub("^qtl2_(.+)_qtl\\.tsv$", "\\1", list.files(jq_dir, pattern = "^qtl2_.+_qtl\\.tsv$"))
if (!length(traits)) stop("no qtl2_<TRAIT>_qtl.tsv in ", jq_dir)
taxon_col <- c(Zd = "#E69F00", Zh = "#6A3D9A", Zl = "#009E73", Zv = "#CC79A7", Zx = "#0072B2",
               Zx_Chal = "#08306B", Zx_Mesa = "#2171B5", Zx_Dura = "#6BAED6", Zx_Nobo = "#9E9AC8", other = "#999999")
taxon_lab <- c(Zd = "Zd diploperennis", Zh = "Zh huehuetenangensis", Zl = "Zl luxurians", Zv = "Zv parviglumis",
               Zx = "Zx mexicana", Zx_Chal = "Zx Chalco", Zx_Mesa = "Zx Mesa Central", Zx_Dura = "Zx Durango",
               Zx_Nobo = "Zx Nobogame", other = "other")

make_fig <- function(trait) {
  sf <- file.path(scan_dir, sprintf("scan_%s.tsv", trait))
  qf <- file.path(jq_dir, sprintf("qtl2_%s_qtl.tsv", trait)); tf <- file.path(jq_dir, sprintf("qtl2_%s_by_taxon.tsv", trait))
  if (!all(file.exists(c(sf, qf, tf)))) { message(trait, ": missing ", paste(basename(c(sf, qf, tf))[!file.exists(c(sf, qf, tf))], collapse = ", ")); return(NULL) }
  sc <- fread(sf); q <- fread(qf); tb <- fread(tf)
  thr <- NA
  smf <- file.path(scan_dir, sprintf("summary_%s.tsv", trait))
  if (file.exists(smf)) { s <- fread(smf); if ("threshold" %in% names(s)) thr <- suppressWarnings(as.numeric(s$threshold[s$scope == "ALL"][1])) }

  ## ---- A: Manhattan ----
  sc[, chrn := as.integer(sub("^chr", "", chr))]
  sc <- sc[order(chrn, pos)]
  len <- sc[, .(len = max(pos)), by = chrn][order(chrn)]
  len[, off := c(0, cumsum(as.numeric(len))[-.N])]; len[, mid := off + len / 2]
  sc <- merge(sc, len[, .(chrn, off)], by = "chrn"); sc[, x := pos + off]
  sc[, y := if ("mlog10p_add" %in% names(sc)) mlog10p_add else -log10(p_add)]
  q[, chrn := as.integer(sub("^chr", "", chr))]
  q <- q[order(chrn, pos_Mb)]; q[, id := paste0("Q", seq_len(.N))]
  q <- merge(q, len[, .(chrn, off)], by = "chrn"); q[, x := pos_Mb * 1e6 + off]
  q[, cond := if ("found_by" %in% names(q)) grepl("^conditional", found_by) else FALSE]
  cut <- -log10(0.05 / nrow(q)); q[, holds := mlog10p_joint > cut]
  q[, id_lab := paste0(id, ifelse(cond, "*", ""))]
  q[, ytop := sapply(seq_len(.N), function(i) { w <- sc[chrn == q$chrn[i] & abs(pos - q$pos_Mb[i] * 1e6) < 1e5]; if (nrow(w)) max(w$y) else NA })]
  pA <- ggplot(sc, aes(x, y)) +
    geom_point(aes(colour = factor(chrn %% 2)), size = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = c("0" = "grey55", "1" = "black")) +
    scale_x_continuous(breaks = len$mid, labels = len$chrn, expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    geom_point(data = q, aes(x, ytop, fill = holds), shape = 25, colour = "#D55E00", size = 2.6, stroke = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white")) +
    geom_text(data = q, aes(x, ytop, label = id_lab), vjust = -1.2, size = 3.2, colour = "#D55E00") +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))),
         title = sprintf("%s: k-mer additive scan (GridLMM, taxon + LOCO K, %s plots)", trait,
                         if ("n_plots" %in% names(sc)) max(sc$n_plots, na.rm = TRUE) else ""),
         caption = sprintf("Filled triangle: the QTL holds in the joint model (joint -log10 p > %.2f, Bonferroni over %d QTL); open: it does not.%s",
                           cut, nrow(q), if (any(q$cond)) "  * found by a conditional scan." else "")) +
    theme_classic(base_size = 11) + theme(plot.title = element_text(size = 11),
                                          plot.caption = element_text(size = 8, colour = "grey30", hjust = 0))
  if (!is.na(thr)) pA <- pA + geom_hline(yintercept = thr, linetype = "dotted") +
    annotate("text", x = min(sc$x), y = thr, label = "5% permutation", hjust = 0, vjust = -0.4, size = 3)

  ## ---- B: joint and per-taxon effects ----
  q[, lab := sprintf("%s  chr%d %.1f Mb", id, chrn, pos_Mb)]
  lev <- rev(q$lab)
  jt <- q[, .(lab, est = effect_per_copy_joint, se, group = "joint (all taxa)")]
  tb[, tx := ifelse(grepl("^other", taxon), "other", taxon)]
  tt <- merge(tb[!is.na(effect_per_copy)], q[, .(qtl, lab)], by = "qtl")[, .(lab, est = effect_per_copy, se, group = tx)]
  d <- rbind(jt, tt); d[, lab := factor(lab, levels = lev)]
  d[, group := factor(group, levels = c("joint (all taxa)", names(taxon_col)))]
  cols <- c("joint (all taxa)" = "black", taxon_col); labs_ <- c("joint (all taxa)" = "joint model (all taxa)", taxon_lab)
  pd <- position_dodge(width = 0.7)
  pB <- ggplot(d, aes(est, lab, colour = group, group = group)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_linerange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), position = pd, linewidth = 0.5) +
    geom_point(aes(shape = group == "joint (all taxa)", size = group == "joint (all taxa)"), position = pd) +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 3.2, `FALSE` = 1.8), guide = "none") +
    scale_colour_manual(values = cols, labels = labs_, drop = TRUE, name = NULL) +
    labs(x = "Effect per teosinte copy (95% CI)", y = NULL, title = "Joint QTL model (R/qtl2) and per-taxon effects") +
    theme_classic(base_size = 11) + theme(plot.title = element_text(size = 11), legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))
  ## -log10 p of the joint test and of "taxa differ", as two labelled columns
  tx <- rbind(q[, .(lab = factor(lab, levels = lev), col = "Joint\ntest", v = sprintf("%.1f", mlog10p_joint), b = holds)],
              q[, .(lab = factor(lab, levels = lev), col = "Taxa\ndiffer",
                    v = ifelse(is.na(mlog10p_taxa_differ), "NA", sprintf("%.1f", mlog10p_taxa_differ)),
                    b = !is.na(mlog10p_taxa_differ) & mlog10p_taxa_differ > cut)])
  pC <- ggplot(tx, aes(col, lab, label = v)) + geom_text(aes(fontface = ifelse(b, "bold", "plain")), size = 3.2) +
    scale_x_discrete(position = "top") +
    labs(title = expression(-log[10](italic(p))), x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 10, hjust = 0.5), axis.text.x.top = element_text(size = 9, lineheight = 0.9))
  pB <- pB + labs(caption = sprintf("Taxa with < %d carrier plots at a QTL have no per-taxon estimate.", 5)) +
    theme(plot.caption = element_text(size = 8, colour = "grey30", hjust = 0))
  (pA | pB | pC) + plot_layout(widths = c(2.3, 1.3, 0.35))
}

figs <- list()
for (t in traits) {
  f <- make_fig(t); if (is.null(f)) next
  h <- max(4.5, 1.2 + 0.45 * nrow(fread(file.path(jq_dir, sprintf("qtl2_%s_qtl.tsv", t)))))
  ggsave(file.path(out_dir, sprintf("joint_qtl_%s.png", t)), f, width = 17, height = h, dpi = 200)
  figs[[t]] <- list(f = f, h = h); message("WROTE ", file.path(out_dir, sprintf("joint_qtl_%s.png", t)))
}
if (length(figs)) {
  pdf(file.path(out_dir, "joint_qtl_all_traits.pdf"), width = 17, height = max(sapply(figs, `[[`, "h")))
  for (x in figs) print(x$f)
  invisible(dev.off()); message("WROTE ", file.path(out_dir, "joint_qtl_all_traits.pdf"))
}
