#!/bin/bash
#SBATCH --job-name=family_donor
#SBATCH --array=1-10
#SBATCH --output=logs/13c_family_donor.%A_%a.out
#SBATCH --error=logs/13c_family_donor.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=4G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 13c — Donor allele of every BC1 family, from the sister-based donor haplotypes
#       (DONOR_ALLELES_DIR/<chr>.donor_alleles.vcf.gz; one column per BC1 family, e.g.
#       Zd.0010_P2, with FORMAT GT:LLR:PRIOR:PALT:SRC). One job per chromosome.
#
#   PALT = probability the family's donor carries the ALT allele at the site.
#   SRC  = own (from that family's own sister data) or gap (filled in).
#   Output: ${D13C}/<FAMILY>/<chr>.tsv.gz    POS REF ALT PALT SRC
#   (REF must be the B73 base; the HMM checks REF/ALT against each marker)
#
#   Submit (independent of 13/13b):  sbatch 13c_family_donor_alleles.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
IN=$(donor_alleles "$CHR")
[[ -s $IN ]] || die "missing ${IN}"
mkdir -p "$D13C"
NAMES=${D13C}/.names_${CHR}
bcftools query -l "$IN" > "$NAMES"
while read -r fam; do mkdir -p "${D13C}/${fam}"; done < "$NAMES"
ulimit -n "$(ulimit -Hn)" 2> /dev/null || true
echo "${CHR}: $(wc -l < "$NAMES") families"

bcftools view -m2 -M2 -v snps -Ou "$IN" \
  | bcftools query -f '%POS\t%REF\t%ALT[\t%PALT\t%SRC]\n' \
  | awk -F'\t' -v OFS='\t' -v chr="$CHR" -v dir="$D13C" -v names="$NAMES" '
      BEGIN { while ((getline s < names) > 0) { n++; f[n] = dir "/" s "/" chr ".tsv" } }
      { for (i = 1; i <= n; i++) {
          p = $(2 + 2 * i); if (p == "." || p == "") continue
          print $1, $2, $3, p, $(3 + 2 * i) > f[i] } }
      END { for (i = 1; i <= n; i++) close(f[i]) }'
while read -r fam; do
  if [[ -s ${D13C}/${fam}/${CHR}.tsv ]]; then gzip -f "${D13C}/${fam}/${CHR}.tsv"; else echo "  ${fam}: no donor alleles on ${CHR}"; fi
done < "$NAMES"
rm -f "$NAMES"
echo "WROTE: ${D13C}/<FAMILY>/${CHR}.tsv.gz"
