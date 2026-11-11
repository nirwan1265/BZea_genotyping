# =============================================================================
# 08_counts_to_ad.awk — convert ANGSD -dumpCounts 4 output to VCF records with
#                       FORMAT/AD = (B73 reads, ALT reads) per pool.
# Called by 08_angsd_counts_by_part.sh. POSIX awk (works with mawk and gawk).
#
# stdin: ANGSD .pos.gz and .counts.gz pasted side by side, headers removed:
#   chr  pos  totDepth  ind0_A ind0_C ind0_G ind0_T  ind1_A ...  (ANGSD leaves a trailing tab)
#
# -v alleles=  07 allele table (gz): CHROM POS B73 ALT SWAP
# -v start= -v end=  part interval (only these allele rows are loaded)
# -v nind=     number of pools (columns) ; -v samples= sample names, one per line
# -v ctrl=     1-based indices of B73 control pools
# -v sampleqc= per-pool totals for this part (tsv)
# -v sitesum=  per-site summary for this part (tsv.gz)
#
# stdout: VCF body lines (header is written by the calling script)
# =============================================================================
BEGIN {
  FS = "\t"; OFS = "\t"
  code["A"] = 0; code["C"] = 1; code["G"] = 2; code["T"] = 3
  base[0] = "A"; base[1] = "C"; base[2] = "G"; base[3] = "T"

  cmd = "gzip -dc " alleles
  while ((cmd | getline line) > 0) {
    split(line, f, "\t")
    p = f[2] + 0
    if (p < start || p > end) continue
    rb[p] = code[f[3]]; ab[p] = code[f[4]]; sw[p] = f[5] + 0
  }
  close(cmd)

  n = 0
  while ((getline line < samples) > 0) name[n++] = line
  close(samples)
  if (n != nind) { print "ERROR: " n " sample names but nind=" nind > "/dev/stderr"; bad = 1; exit 2 }

  for (i = 0; i < nind; i++) isctrl[i] = 0
  while ((getline line < ctrl) > 0) if (line != "") isctrl[line - 1] = 1
  close(ctrl)

  need = 3 + 4 * nind
  sitecmd = "gzip -c > " sitesum
  print "CHROM", "POS", "B73", "ALT", "DP", "OTH", "REF_READS", "ALT_READS", "NCOV", "NALT", \
        "CTRL_REF", "CTRL_ALT", "CTRL_NCOV" | sitecmd
}

{
  if (NF < need) {
    print "ERROR: input line " NR " has " NF " fields, expected >= " need > "/dev/stderr"
    bad = 1; exit 3
  }
  p = $2 + 0
  if (!(p in rb)) { nmiss++; next }
  r = rb[p]; a = ab[p]

  dp = 0; ot = 0; tr = 0; ta = 0; nc = 0; na = 0; cr = 0; ca = 0; cn = 0
  for (i = 0; i < nind; i++) {
    k = 4 + 4 * i                       # field of ind<i>_A
    x = $(k + r) + 0; y = $(k + a) + 0
    t = $k + $(k + 1) + $(k + 2) + $(k + 3)
    o = t - x - y
    R[i] = x; A[i] = y
    dp += t; ot += o; tr += x; ta += y
    if (x + y > 0) { nc++; qc[i]++ }
    if (y > 0) na++
    qr[i] += x; qa[i] += y; qo[i] += o
    if (isctrl[i]) { cr += x; ca += y; if (x + y > 0) cn++ }
  }

  printf "%s\t%d\t.\t%s\t%s\t.\t.\tDP=%d;OTH=%d;NCOV=%d;NALT=%d%s\tAD", \
         $1, p, base[r], base[a], dp, ot, nc, na, (sw[p] == 1 ? ";SWAP" : "")
  for (i = 0; i < nind; i++) printf "\t%d,%d", R[i], A[i]
  printf "\n"

  print $1, p, base[r], base[a], dp, ot, tr, ta, nc, na, cr, ca, cn | sitecmd
  nout++
}

END {
  if (bad) exit 3
  close(sitecmd)
  print "SAMPLE", "SITES_COVERED", "REF_READS", "ALT_READS", "OTHER_READS" > sampleqc
  for (i = 0; i < nind; i++) print name[i], qc[i] + 0, qr[i] + 0, qa[i] + 0, qo[i] + 0 > sampleqc
  close(sampleqc)
  printf "sites written: %d ; ANGSD sites missing from allele table: %d\n", nout, nmiss + 0 > "/dev/stderr"
}
