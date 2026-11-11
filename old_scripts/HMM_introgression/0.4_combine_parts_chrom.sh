#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara 
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J chr10
#BSUB -o stdout.%J
#BSUB -e stderr.%J

set -euo pipefail

ROOT="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps"
OUTDIR="${ROOT}/per_chr"

mkdir -p "${OUTDIR}"
cd "${ROOT}"

for c in {10..10}; do
  chr="chr${c}"
  out="${OUTDIR}/BZea_${chr}.full.bcf"

  echo "==== ${chr} ===="

  # Grab all parts for this chr
  parts=( BZea_${chr}.part*.bcf )

  if [[ ${#parts[@]} -eq 0 ]]; then
    echo "WARNING: no parts found for ${chr} (pattern BZea_${chr}.part*.bcf)" >&2
    continue
  fi

  # Sort by numeric START from filename: ...partN.START_END.bcf
  # Output: "START<TAB>filename" -> sort -> filenames
  mapfile -t sorted_parts < <(
    printf '%s\n' "${parts[@]}" \
    | awk -F'.' '
        {
          # field with START_END is second-to-last: ... .START_END .bcf
          se=$(NF-1)
          split(se,a,"_")
          start=a[1]+0
          print start "\t" $0
        }' \
    | sort -n -k1,1 \
    | cut -f2-
  )

  echo "Found ${#sorted_parts[@]} parts for ${chr}"
  echo "First: ${sorted_parts[0]}"
  echo "Last : ${sorted_parts[-1]}"

  # Build a list file for bcftools concat (safer than long argv)
  listfile="$(mktemp)"
  printf "%s\n" "${sorted_parts[@]}" > "${listfile}"

  # Concatenate: parts are non-overlapping, ordered by coordinate
  bcftools concat -f "${listfile}" -Ob -o "${out}"
  bcftools index -f "${out}"

  rm -f "${listfile}"

  echo "WROTE: ${out}"
done

echo "DONE all chromosomes."
