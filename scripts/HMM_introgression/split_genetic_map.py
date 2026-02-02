#!/usr/bin/env python3
"""
Split NAM_genetic_map.txt into per-chromosome .map files for HMM.

Usage: python split_genetic_map.py NAM_genetic_map.txt output_dir
"""

import os
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python split_genetic_map.py <input_map> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)

    # Read and parse the input file
    chr_data = {f"chr{i}": [] for i in range(1, 11)}

    with open(input_file, "r") as f:
        header = next(f)  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            try:
                bp = int(parts[0])
                chrom_num = int(parts[1])
                cm = float(parts[2])

                chrom = f"chr{chrom_num}"
                if chrom in chr_data:
                    chr_data[chrom].append((bp, cm))
            except ValueError:
                continue

    # Write per-chromosome map files
    # Format expected by HMM: first column = bp, last column = cM
    for chrom, data in chr_data.items():
        if not data:
            print(f"Warning: No data for {chrom}")
            continue

        # Sort by bp position
        data.sort(key=lambda x: x[0])

        out_file = os.path.join(output_dir, f"{chrom}.map")
        with open(out_file, "w") as f:
            for bp, cm in data:
                f.write(f"{bp}\t{cm}\n")

        print(f"Wrote {len(data)} markers to {out_file}")

    print(f"\nDone! Map files written to {output_dir}")

if __name__ == "__main__":
    main()
