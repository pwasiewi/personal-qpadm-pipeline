# ══════════════════════════════════════════════════════════════════════════════
# OPTIONAL: merge RAW data from multiple platforms before step 01
# Script: aadr/00_merge_arrays.sh <plink_A> <plink_B> <outdir>
# Input:  two pw_plink filesets after step 01 (same person, different platforms)
# Output: <outdir>/pw_plink.{bed,bim,fam} — SNP union, then normal pipeline
#
# Supported CSV formats (step 01):
#   FTDNA Family Finder : RSID,CHROMOSOME,POSITION,RESULT  (no ## header)
#   MyHeritage          : ## reference=build37 + RSID,CHROMOSOME,POSITION,RESULT
#   Ancestry            : rsid,chromosome,position,allele1,allele2  (two allele columns!)
#                         → step 01 auto-detects format from header
#
# Two arrays of the same person → merge at step 00:
#   ./aadr/01_prepare_files.sh /path/to/FTDNA_PW.csv     PERSON_FTDNA
#   ./aadr/01_prepare_files.sh /path/to/MyHeritage_PW.csv PERSON_MH
#   SAMPLE_ID=PW ./aadr/00_merge_arrays.sh \
#       PERSON_FTDNA/pw_plink PERSON_MH/pw_plink PERSON_MERGED
#   # continue pipeline from step 02 with ./PERSON_MERGED as outdir
#
# Three arrays (FTDNA + MyHeritage + Ancestry) — iteratively:
#   SAMPLE_ID=PW ./aadr/00_merge_arrays.sh \
#       PERSON_FTDNA/pw_plink PERSON_MH/pw_plink PERSON_MH_FF
#   SAMPLE_ID=PW ./aadr/00_merge_arrays.sh \
#       PERSON_MH_FF/pw_plink PERSON_ANC/pw_plink PERSON_MERGED
#
# Practical notes (from experience merging MH_PW + PW_FF):
#   - Overlap MH ∩ FTDNA: ~556K / ~600K SNPs (92%) — arrays overlap heavily
#   - Union MH+FTDNA: ~612K SNPs, but after extract_f2(maxmiss=0) FEWER than either alone
#     (~49,218 vs ~50,114 / ~49,883) — edge-of-overlap missing pairs reduce yield
#   - Ancestry DNA (Illumina OmniExpress v2.5, ~730K SNPs) has different content than FTDNA/MH;
#     three-platform union should yield more unique SNPs overlapping AADR
#   - Ancestry format: two allele columns (allele1, allele2) — step 01 must concatenate them
#     into a single RESULT string before conversion; TODO: add support to 01_prepare_files.sh
#   - Allele conflicts (missnp): typically 1–10 SNPs — removed automatically by step 00
#   - "same-position variants" warnings (rs8192284 / rs2228145 etc.) are normal —
#     different rsIDs at the same position; do not block the merge
# ══════════════════════════════════════════════════════════════════════════════

# ── PW (FTDNA Family Finder) ─────────────────────────────────────────────────
./aadr/01_prepare_files.sh /usr/local/share/aadr/FF_PW_Chrom_Autoso_20260223.csv PW_FF
Rscript aadr/02_build_aadr_subset.R /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB ./PW_FF
./aadr/03_merge_and_test.sh ./PW_FF
Rscript aadr/04_run_qpadm_test.R ./PW_FF/merged_pw_aadr_final_ready PW
Rscript aadr/05_robustness_outgroups.R ./PW_FF/merged_pw_aadr_final_ready /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB PW
Rscript aadr/06_ancestry_models.R ./PW_FF/merged_pw_aadr_extended /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB PW
Rscript aadr/07_extended_models.R ./PW_FF/merged_pw_aadr_final_ready /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB PW
Rscript aadr/08_two_source.R     ./PW_FF/merged_pw_aadr_final_ready /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB PW

# ── JW (MyHeritage) — 5 differences vs PW ───────────────────────────────────
# 1. SAMPLE_ID=JW   (env var for step 1; arg 2 for steps 3-8)
# 2. different CSV  (MyHeritage_JW_raw_dna_data.csv)
# 3. different outdir (JW_MH instead of PW_FF)
# 4. step 3 requires 2nd argument JW  (default "PW" → grep fails in step 3g)
# 5. steps 4-8 use JW as target_id
SAMPLE_ID=JW ./aadr/01_prepare_files.sh /usr/local/share/aadr/MyHeritage_JW_raw_dna_data.csv JW_MH
Rscript aadr/02_build_aadr_subset.R /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB ./JW_MH
./aadr/03_merge_and_test.sh ./JW_MH JW
Rscript aadr/04_run_qpadm_test.R ./JW_MH/merged_pw_aadr_final_ready JW
Rscript aadr/05_robustness_outgroups.R ./JW_MH/merged_pw_aadr_final_ready /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB JW
Rscript aadr/06_ancestry_models.R ./JW_MH/merged_pw_aadr_extended /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB JW
Rscript aadr/07_extended_models.R ./JW_MH/merged_pw_aadr_final_ready /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB JW
Rscript aadr/08_two_source.R     ./JW_MH/merged_pw_aadr_final_ready /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB JW

# ══════════════════════════════════════════════════════════════════════════════
# STEP 07 (extended 2026-06-28): 23 single-source populations — Slavic/Baltic/
#   Lusatian + Pannonian Avars/Huns/Sarmatians/Steppe Bulgars. Verification of
#   G25 "Avar_Kecskemet": Pannonian Avars REJECTED as single source for PW.
# STEP 08 (2026-06-28): rotating 2-source scan — 2 base sources × 23 diverse
#   second sources (Uralic/Siberian/Steppe/Hunnic/Middle Eastern/Caucasian).
#   SIGNIFICANT-2src verdict detects a real additional component;
#   UNIDENTIFIABLE flags weight blowup (Poland_EarlySlav base).
#   Result: 0 SIGNIFICANT-2src for both targets → no hidden Avar/steppe admixture;
#   Avar_Kecskemet = artefact same as Germany_IA. f2caches and TSVs are per-target
#   (PW/JW do not collide). options(width=200) in scripts 06/07/08 — tables do not
#   wrap columns.
# Full results: PW_FF_slavic_aadr.md §5A/5B/§10, JW_MH_slavic_aadr.md §4A/4B
# ══════════════════════════════════════════════════════════════════════════════
