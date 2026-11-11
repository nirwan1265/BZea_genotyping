#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara 
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J Combine_subparts
#BSUB -o stdout.%J
#BSUB -e stderr.%J

set -euo pipefail

SRC="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/remaining"
DEST="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps"

#mkdir -p "${DEST}"

cd "${SRC}"

# Loop over all sub1 BCFs (each should have a matching sub2)
shopt -s nullglob
for f1 in BZea_chr*.part*.sub1.*_*.bcf; do
  # Parse: BZea_chr1.part4.sub1.46267873_53979184.bcf
  if [[ "${f1}" =~ ^(BZea_chr[^.]+\.part[0-9]+)\.sub1\.([0-9]+)_([0-9]+)\.bcf$ ]]; then
    prefix="${BASH_REMATCH[1]}"
    s1="${BASH_REMATCH[2]}"
    e1="${BASH_REMATCH[3]}"

    # Find matching sub2
    f2_glob="${prefix}.sub2.*_*.bcf"
    f2=( ${f2_glob} )
    if [[ ${#f2[@]} -ne 1 ]]; then
      echo "ERROR: expected exactly 1 match for ${f2_glob}, found ${#f2[@]}" >&2
      printf '%s\n' "${f2[@]:-}" >&2
      exit 2
    fi
    f2="${f2[0]}"

    if [[ "${f2}" =~ ^${prefix}\.sub2\.([0-9]+)_([0-9]+)\.bcf$ ]]; then
      s2="${BASH_REMATCH[1]}"
      e2="${BASH_REMATCH[2]}"
    else
      echo "ERROR: could not parse sub2 filename: ${f2}" >&2
      exit 2
    fi

    # Sanity checks: sub1 should start before sub2; intervals should not overlap weirdly
    if (( s2 <= s1 )); then
      echo "ERROR: unexpected ordering: ${f1} and ${f2}" >&2
      exit 2
    fi

    out="${prefix}.${s1}_${e2}.bcf"

    echo "----"
    echo "Merging:"
    echo "  ${f1}"
    echo "  ${f2}"
    echo "=> ${out}"

    # Concatenate (intervals are disjoint/adjacent, so concat is correct)
    bcftools concat -a -Ob -o "${out}" "${f1}" "${f2}"

    # Index (bcftools will usually create .csi for BCF)
    bcftools index -f "${out}"

    # Move merged BCF + index to DEST
    if [[ -f "${out}.csi" ]]; then
      mv -v "${out}" "${out}.csi" "${DEST}/"
    elif [[ -f "${out}.tbi" ]]; then
      mv -v "${out}" "${out}.tbi" "${DEST}/"
    else
      echo "WARNING: index not found for ${out} (expected .csi or .tbi)" >&2
      mv -v "${out}" "${DEST}/"
    fi

  else
    echo "SKIP (unexpected name): ${f1}" >&2
  fi
done

echo "DONE."
