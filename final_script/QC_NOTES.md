# QC notes and decisions (read counts, steps 07–10)

Run date: 2026-11-13. Source: `10_qc/sample_qc.tsv`, `10_qc/taxon_summary.tsv`, `10_qc/site_qc.tsv`.
The sample lists referred to below are in [`sample_lists/`](sample_lists/).

## Overall

| | |
|---|---|
| SNP sites with read counts | 9,675,593 (chr1–10) |
| Pools in the run | 1,535 (plates PN2–PN17; PN1 and PN18 were not in the BAM list) |
| Reads per site, all pools together | ~1,050 (≈ 0.6–0.7 per pool) |
| Share of sites with ≥1 read, per pool | ~38% (median) |
| Third-base reads (error proxy) | ~0.18% |
| ALT fraction, clean B73 controls | 0.39% (background) |
| ALT fraction, donor pools (median) | Zd 1.1%, Zl 1.7%, Zh 1.8%, Zv 1.9%, Zx 2.1% |

A donor pool's ALT fraction minus the B73 background, divided by the expected donor share (12.5%), shows how many of the 9.7 M SNPs a donor haplotype actually carries. That is about 11–13% for Zv, Zx, Zl and Zh, and about 5.5% for Zd. For any single line, most sites tell us nothing about its donor, which is why step 2 builds site sets per donor group.

## Decisions

1. **Two B73 controls are not clean B73 → removed from the controls.**
   PN3_SID236 (ALT 1.76%) and PN5_SID468 (1.84%) look like donor lines, not B73, so they are mislabeled or contaminated. The other 10 controls (`b73_controls_clean.txt`) are at 0.35–0.53%. Only those 10 are used as controls. List: `exclude_b73_not_clean.txt`.

2. **PN13_SID1225 is PurpleCheck, not Zd.0020 → excluded.**
   The metadata lists it twice, as `Zd.0020_P3_P1_P2.5.1.1` and as `PurpleCheck`. Its ALT fraction is 20.2%, the same as the 18 PurpleCheck pools (median 20.4%).
   It cannot be a Zd line. Zd pools have a median ALT of 1.1% and none is above 6.4%. Even a line that was 100% diploperennis would only reach about 0.4% + 5.5% ≈ 6%, so 20% is out of reach for any Zd BC2S3. It also has low coverage (1.27 M reads). List: `exclude_mislabeled_purplecheck.txt`.

3. **65 donor pools look like pure B73 → kept, watched.**
   These have ALT ≤ 0.64%, the clean-B73 maximum + 20%: 31 Zd, 22 Zx, 11 Zv and 1 Zh. Zd.0010 alone accounts for 13 of its 87 pools.
   There are three possible explanations:
   - the line lost its donor segments
   - a failed cross or a selfed B73
   - (for Zd) diploperennis reads mapping too poorly to show up

   The step-3 caller should tell these apart. List: `watch_donor_looks_b73.txt`.

4. **9 donor pools are unusually high (6–11% ALT) → kept, watched.**
   Eight are Zx and one is Zd. They may carry more donor genome than expected, or be outcrossed or mixed seed. List: `watch_donor_high_alt.txt`.

5. **69 pools have no metadata → excluded for now.**
   They are on plates PN2 (62), PN3 (6) and PN13 (1). Their ALT fractions (1.6–6%) look like donor lines. Add them back once their pedigrees are found. List: `exclude_no_metadata.txt`.

6. **`Zm` (29 pools) = probably maize landrace donors → their own donor group.**
   They have the highest ALT fraction (4.7%) and a different pedigree format (`Zm.0950_P1_P1_1`).

7. **25 low-coverage pools (< 20% of the median read count, down to 5,541 reads) → excluded.** List: `exclude_low_coverage.txt`.

8. **Sites with ALT reads in the B73 controls → removed in step 2.**
   About 3.5% of sites showed ALT reads in the original 12 controls. Sequencing errors alone would explain about 0.4%, so most of these are likely paralogs, repeats or B73 reference problems. Step 2 recomputes this with the 10 clean controls only.

`exclude_all.txt` is the union of the lists from points 1, 2, 5 and 7: 93 pools, because some pools are on more than one list. Those pools are left out of steps 2 and 3.

## Teosinte reference panel

- **SNPversity 2.0** (`/rsstu/users/r/rrellan/sara/ref/SNPversity`, 1,498 accessions) is almost all maize inbreds, with only a few teosinte (e.g. TIL05, TIL09). **Not used.**
- **Schnable 2023 teosinte subset** (`/rsstu/users/r/rrellan/sara/ref/Schnable_SNPset/teosinte_only`) is 227 teosinte samples on B73 v5 with no missing calls (AN = 454). **Used** as the teosinte panel, with all samples pooled into one group, because no subspecies labels are available yet.
- **Overlap check on chr10:**
  - The panel has 3.21 M biallelic SNPs; our site set has 735,015.
  - 251,075 of our sites (34%) are in the panel.
  - 250,271 of those (99.7%) have the same B73/ALT alleles, which independently confirms our sites and their orientation to B73.
