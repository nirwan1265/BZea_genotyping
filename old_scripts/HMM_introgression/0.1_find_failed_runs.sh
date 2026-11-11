cd /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs

# 1) List the 30 failing stderr files (the ones that contain the lsbatch "Killed" line)
grep -l "/home/ntanduk/.lsbatch/" *.err > failed.err.list

# sanity check (should be 30)
wc -l failed.err.list

# 2) Also make a list based on the explicit word "Killed" (sometimes nicer)
grep -l "Killed" *.err > failed.killed.list
wc -l failed.killed.list


awk -F'.' '{print $(NF-1), $0}' failed.err.list | sort -n > failed_indices_and_files.tsv
head failed_indices_and_files.tsv


while read -r errf; do
  base="${errf%.err}"
  outf="${base}.out"
  printf "\n=== %s ===\n" "$errf"
  if [[ -f "$outf" ]]; then
    grep -E "CHR=|INTERVAL=|PART=|OUT=" "$outf" | head -n 20
  else
    echo "Missing matching stdout: $outf"
  fi
done < failed.err.list > failed_jobs_debug.txt





cd /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs

# Input: failed_indices_and_files.tsv (you already made this)
# Output: failed_jobs_chr_interval.tsv
out_tsv="failed_jobs_chr_interval.tsv"

{
  echo -e "index\tjobid\terr_file\tout_file\tCHR\tPART\tINTERVAL\tOUT"
  while read -r idx errf; do
    # errf like: angsd.52149.1.err
    jobid=$(echo "$errf" | awk -F'.' '{print $2}')
    out="${errf%.err}.out"

    CHR="NA"; PART="NA"; INTERVAL="NA"; OUTP="NA"

    if [[ -f "$out" ]]; then
      CHR=$(grep -m1 '^CHR=' "$out" | sed 's/^CHR=//' || true)
      PART=$(grep -m1 '^PART=' "$out" | sed 's/^PART=//' || true)
      INTERVAL=$(grep -m1 'INTERVAL=' "$out" | sed -E 's/.*INTERVAL=([^[:space:]]+).*/\1/' || true)
      OUTP=$(grep -m1 '^OUT=' "$out" | sed 's/^OUT=//' || true)
      [[ -z "$CHR" ]] && CHR="NA"
      [[ -z "$PART" ]] && PART="NA"
      [[ -z "$INTERVAL" ]] && INTERVAL="NA"
      [[ -z "$OUTP" ]] && OUTP="NA"
    fi

    echo -e "${idx}\t${jobid}\t${errf}\t${out}\t${CHR}\t${PART}\t${INTERVAL}\t${OUTP}"
  done < failed_indices_and_files.tsv
} > "$out_tsv"

echo "Wrote: $out_tsv"
wc -l "$out_tsv"
