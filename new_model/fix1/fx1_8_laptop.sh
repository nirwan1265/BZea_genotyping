#!/bin/bash
# =============================================================================
# fix1 / 8 — the laptop steps on the new model, after fx1_7 has finished and its outputs were copied here
# (rsync commands in README.md of this folder). Run from the repo root:  bash new_model/fix1/fx1_8_laptop.sh
#   30   genes + GO in the GWAS regions (taxon and mexicana races)
#   33a  genes per joint-model QTL      33  SuSiE fine-mapping (PIP)      33b  PIPs into the gene tables
#   31   GWAS + joint-model figures     31b  DTS by year      24b  introgression landscape      34e  targeted-seq figure
#   32c  k-mer PCA vs the new calls     checks: tract-length simulation, ra2 x tsh4, NBR Q8 vs flowering
# Inputs: data/introgression_fix1/ (calls, scans, joint QTL, step-25 inputs in assoc/, targseq/).
# Outputs: data/introgression_fix1/{GO,GO_subtaxon,qtl_genes,finemap_susie,checks,sourmash/pca}, figures/fix1/.
# The current model's files (data/introgression_spread, qtl_genes/, figures/) are not touched.
# Needs: R (susieR, data.table, tidyverse, ggplot2, patchwork); python (numpy, pandas, matplotlib, statsmodels).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../.."
D=data/introgression_fix1; OLD=data/introgression_spread; F=figures/fix1
TRAITS="DTS DTA ASI NBR SPAD2 PH EH EN Prolif SL BL BW"
for x in tracts_all.tsv.gz pool_summary.tsv assoc/geno_dosage.csv assoc_gridlmm_spats/genes joint_qtl/noperm targseq/sites_used.tsv.gz; do
  [[ -e $D/$x ]] || { echo "missing $D/$x: copy the cluster outputs first (README of new_model/fix1)"; exit 1; }
done
mkdir -p $D/GO $D/GO_subtaxon $D/checks $F/subtaxon $F/results
cp $OLD/GO/GO_terms_all_candidates.txt $D/GO/                                   # GO annotation (input, unchanged)
echo "== 30 genes + GO"
Rscript final_script/30_genes_go.R $D/assoc_gridlmm_spats/genes $D/GO/GO_terms_all_candidates.txt $D/GO
Rscript final_script/30_genes_go.R $D/assoc_gridlmm_spats_subtaxon/genes $D/GO/GO_terms_all_candidates.txt $D/GO_subtaxon
echo "== 33a genes per QTL"
python final_script/33a_qtl_gene_tables.py --qtl_dir $D/joint_qtl/noperm --scan_dir $D/assoc_gridlmm_spats \
  --go_all $D/GO/genes_GO_all.tsv $D/GO_subtaxon/genes_GO_all.tsv --out $D/qtl_genes --traits "$TRAITS"
echo "== 33 SuSiE"
for T in $TRAITS; do
  [[ -s $D/joint_qtl/noperm/qtl2_${T}_qtl.tsv ]] || continue
  Rscript final_script/33_susie_finemap.R $D/assoc $T $D/joint_qtl/noperm/qtl2_${T}_qtl.tsv $D/finemap_susie $D/qtl_genes/${T}_qtl_genes.tsv
done
python final_script/33b_pip_to_gene_tables.py --susie $D/finemap_susie --genes $D/qtl_genes
echo "== 31 figures"
Rscript final_script/31_joint_qtl_figures.R $D/assoc_gridlmm_spats $D/joint_qtl/noperm $F
Rscript final_script/31_joint_qtl_figures.R $D/assoc_gridlmm_spats_subtaxon $D/joint_qtl/noperm_subtaxon $F/subtaxon
python final_script/31b_year_compare.py --trait DTS --base $D --out $F/DTS_by_year.png
echo "== 24b landscape, 34e targeted sequencing, 32c PCA"
python final_script/24b_introgression_landscape.py --tracts $D/tracts_all.tsv.gz --pool_summary $D/pool_summary.tsv --out $F/results
python final_script/34e_targseq_figure.py --sites $D/targseq/sites_used.tsv.gz --tracts $D/targseq/tracts_targseq.tsv --out $F/targseq_validation.png
mkdir -p $D/sourmash/pca; POOL_SUMMARY=$D/pool_summary.tsv Rscript final_script/32c_kmer_pca.R $OLD/sourmash/kmer_taxon_profile.tsv $D/sourmash/pca data/metadata_BZea.csv
echo "== checks"
python final_script/tract_length_sim.py --map_dir final_script/genetic_maps --tracts $D/tracts_all.tsv.gz --pool_summary $D/pool_summary.tsv > $D/checks/tract_length_sim.txt
python final_script/ra2_tsh4_interaction.py --tracts $D/tracts_all.tsv.gz --pool_summary $D/pool_summary.tsv --pheno data/phenotypes/BZea_spats_BLUEs_combined.csv > $D/checks/ra2_tsh4.txt
python final_script/nbr_q8_flowering_check.py --tracts $D/tracts_all.tsv.gz --pool_summary $D/pool_summary.tsv --pheno data/phenotypes/BZea_spats_BLUEs_combined.csv > $D/checks/nbr_q8.txt
echo "done: $D and $F"
