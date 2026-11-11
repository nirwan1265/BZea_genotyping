#!/bin/bash
#SBATCH --job-name=sites
#SBATCH --array=1-10
#SBATCH --output=logs/07_sites.%A_%a.out
#SBATCH --error=logs/07_sites.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 07 — Build the site list for read counting, one job per chromosome.
#
#   Takes ALL biallelic SNPs ANGSD discovered (the old 0.6 filter is NOT applied:
#   its population-wide AF >= 0.05 cut removes donor-specific SNPs; site
#   filtering now happens in step 2 using read counts, B73 controls and donor
#   groups). Every SNP is re-polarized to the B73 base (ANGSD wrote major,minor).
#
#   Per chromosome:
#     BZea_<chr>.alleles.tsv.gz   CHROM POS B73 ALT SWAP   (no header)
#     BZea_<chr>.sites (+ .bin/.idx)  ANGSD -sites file
#     BZea_<chr>.polarize.log     kept / swapped / dropped counts
#   Shared (identical in every job):
#     samples.txt         PN*_SID* in ANGSD column order (checked against the BCF)
#     b73_controls.idx    1-based positions of the B73 control pools in samples.txt
#     ad_header.vcf       VCF header used by step 08
#
#   Submit: sbatch 07_make_site_tables.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
BCF=$(discovery_bcf "$CHR")
[[ -s $BCF ]] || die "missing discovery BCF ${BCF} (check DISCOVERY_BCF_DIR in config.sh)"
[[ -s ${REF}.fai ]] || die "missing ${REF}.fai"
[[ -s $SAMPLE_META ]] || die "missing ${SAMPLE_META}"
mkdir -p "$D07"

# write stdin to $1 atomically (every array job writes the same shared files)
write_atomic() { local tmp="$1.tmp.$$"; cat > "$tmp"; mv -f "$tmp" "$1"; }

# ---- samples: ANGSD column order = BAM-list order --------------------------
sample_names_from_bamlist | write_atomic "${D07}/samples.txt"
bcftools query -l "$BCF" | cmp -s - "${D07}/samples.txt" \
  || die "sample names/order in ${BCF} differ from the BAM list"
NIND=$(wc -l < "${D07}/samples.txt")

# ---- B73 control pools (metadata column 2 == "B73"; file has CRLF endings) --
tr -d '\r' < "$SAMPLE_META" | awk -F'\t' 'NR > 1 && $2 == "B73" {print $1}' | sort -u > "${D07}/.b73_ids.$$"
awk 'NR == FNR {c[$1] = 1; next} ($1 in c) {print FNR}' "${D07}/.b73_ids.$$" "${D07}/samples.txt" \
  | write_atomic "${D07}/b73_controls.idx"
NCTRL_META=$(wc -l < "${D07}/.b73_ids.$$")
NCTRL=$(wc -l < "${D07}/b73_controls.idx")
rm -f "${D07}/.b73_ids.$$"
echo "B73 control pools: ${NCTRL} found in BAM list (${NCTRL_META} in metadata)"
(( NCTRL > 0 )) || die "no B73 control pools found"

# ---- VCF header for step 08 --------------------------------------------------
{
  echo '##fileformat=VCFv4.2'
  echo "##reference=file://${REF}"
  awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${REF}.fai"
  echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Reads at this site over all pools (A+C+G+T, after ANGSD read filters)">'
  echo '##INFO=<ID=OTH,Number=1,Type=Integer,Description="Reads with a base that is neither B73 nor ALT, all pools">'
  echo '##INFO=<ID=NCOV,Number=1,Type=Integer,Description="Pools with >=1 B73 or ALT read">'
  echo '##INFO=<ID=NALT,Number=1,Type=Integer,Description="Pools with >=1 ALT read">'
  echo '##INFO=<ID=SWAP,Number=0,Type=Flag,Description="ANGSD major allele was not the B73 base; alleles re-polarized">'
  echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read counts: B73 (REF) allele, ALT allele">'
  echo '##source=BZea_genotyping/final_script/08 (ANGSD -doCounts 1 -dumpCounts 4)'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
  awk '{printf "\t%s", $0} END {printf "\n"}' "${D07}/samples.txt"
} | write_atomic "${D07}/ad_header.vcf"

# ---- per-chromosome site table ----------------------------------------------
ALLELES=${D07}/BZea_${CHR}.alleles.tsv.gz
SITES=${D07}/BZea_${CHR}.sites

bcftools view -m2 -M2 -v snps "$BCF" -Ou \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  | python3 "${FINAL_DIR}/07_polarize_to_b73.py" --fasta "$REF" --log "${D07}/BZea_${CHR}.polarize.log" \
  | bgzip -c > "$ALLELES"

gzip -dc "$ALLELES" | cut -f1,2 > "$SITES"
angsd sites index "$SITES"

echo "---- ${CHR} ----"
cat "${D07}/BZea_${CHR}.polarize.log"
echo "pools: ${NIND}"
echo "WROTE: ${ALLELES} ${SITES}"
