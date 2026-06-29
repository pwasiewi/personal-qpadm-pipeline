#!/usr/bin/env bash
# ==============================================================================
# 00_merge_arrays.sh -- merge two consumer-array pw_plink filesets for the
#                       same individual into one combined PLINK fileset.
#
# Use before step 02 when you have both MyHeritage and FTDNA raw data for the
# same person and want the union of SNPs for a higher-power qpAdm run.
#
# Usage:
#   ./aadr/00_merge_arrays.sh <plink_prefix_A> <plink_prefix_B> <outdir>
#
# Example:
#   ./aadr/00_merge_arrays.sh MH_PW/pw_plink FF_PW/pw_plink PW_MH_FF
#
# Both inputs must have the same FID/IID (same individual).
# The script handles strand conflicts with the same probe+exclude logic as
# step 3d (--flip does NOT fix one-sided missingness conflicts).
#
# Output:
#   <outdir>/pw_plink.{bed,bim,fam}   -- ready as input for step 02
# ==============================================================================
set -euo pipefail

A_PREFIX="${1:?Usage: $0 <plink_A> <plink_B> <outdir>}"
B_PREFIX="${2:?Usage: $0 <plink_A> <plink_B> <outdir>}"
OUTDIR="${3:?Usage: $0 <plink_A> <plink_B> <outdir>}"

mkdir -p "$OUTDIR"
WORK="$OUTDIR/_merge_work"
mkdir -p "$WORK"

echo "================================================================================"
echo "  STEP 0: Merge two consumer arrays for the same individual"
echo "================================================================================"
echo "  Array A : $A_PREFIX  ($(wc -l < "${A_PREFIX}.bim") SNPs)"
echo "  Array B : $B_PREFIX  ($(wc -l < "${B_PREFIX}.bim") SNPs)"
echo "  Outdir  : $OUTDIR"
echo

# ---- Idempotency ----
if [[ -f "$OUTDIR/pw_plink.bed" ]]; then
  echo "  $OUTDIR/pw_plink.bed already exists -- skipping merge."
  echo "  Delete to force re-run."
  echo "================================================================================"
  exit 0
fi

# ---- 0a. Sort both filesets ----
echo "[0a] Sorting variants in both arrays..."
plink2 --bfile "$A_PREFIX" --sort-vars --make-pgen --out "$WORK/a_pgen"
plink2 --pfile "$WORK/a_pgen" --make-bed --out "$WORK/a_sorted"

plink2 --bfile "$B_PREFIX" --sort-vars --make-pgen --out "$WORK/b_pgen"
plink2 --pfile "$WORK/b_pgen" --make-bed --out "$WORK/b_sorted"

# ---- 0b. Probe merge to find strand/allele conflicts ----
echo
echo "[0b] Probe merge to discover allele conflicts..."
echo "     (A/T and C/G SNPs on mismatched strands appear as 3+-allele conflicts)"
set +e
p-link --bfile "$WORK/a_sorted" \
       --bmerge "$WORK/b_sorted.bed" "$WORK/b_sorted.bim" "$WORK/b_sorted.fam" \
       --make-bed --out "$WORK/probe" 2>&1 | tee "$WORK/probe_attempt.log"
set -e

if [[ -f "$WORK/probe-merge.missnp" ]]; then
  n_conflict="$(wc -l < "$WORK/probe-merge.missnp")"
  echo "  Found $n_conflict conflict SNPs. Excluding from both arrays."
  plink2 --bfile "$WORK/a_sorted" --exclude "$WORK/probe-merge.missnp" \
         --make-bed --out "$WORK/a_final"
  plink2 --bfile "$WORK/b_sorted" --exclude "$WORK/probe-merge.missnp" \
         --make-bed --out "$WORK/b_final"
else
  echo "  No conflicts -- using sorted filesets directly."
  cp "$WORK/a_sorted.bed" "$WORK/a_final.bed"
  cp "$WORK/a_sorted.bim" "$WORK/a_final.bim"
  cp "$WORK/a_sorted.fam" "$WORK/a_final.fam"
  cp "$WORK/b_sorted.bed" "$WORK/b_final.bed"
  cp "$WORK/b_sorted.bim" "$WORK/b_final.bim"
  cp "$WORK/b_sorted.fam" "$WORK/b_final.fam"
fi

# ---- 0c. Real merge ----
echo
echo "[0c] Merging arrays with p-link --bmerge..."
set +e
p-link --bfile "$WORK/a_final" \
       --bmerge "$WORK/b_final.bed" "$WORK/b_final.bim" "$WORK/b_final.fam" \
       --make-bed --out "$WORK/merged" 2>&1 | tee "$WORK/merge_final.log"
merge_exit=$?
set -e

if [[ $merge_exit -eq 0 && -f "$WORK/merged.bim" ]]; then
  echo "  Merge completed (main --make-bed succeeded)."
  cp "$WORK/merged.bed" "$OUTDIR/pw_plink.bed"
  cp "$WORK/merged.bim" "$OUTDIR/pw_plink.bim"
  cp "$WORK/merged.fam" "$OUTDIR/pw_plink.fam"
elif [[ -f "$WORK/merged-merge.bim" ]]; then
  echo "  --make-bed exit $merge_exit but intermediate found -- converting via plink2."
  plink2 --bfile "$WORK/merged-merge" --make-bed --out "$OUTDIR/pw_plink"
else
  echo "ERROR: merge failed and no intermediate found." >&2
  echo "       Check $WORK/merge_final.log" >&2
  exit 1
fi

n_a="$(wc -l < "${A_PREFIX}.bim")"
n_b="$(wc -l < "${B_PREFIX}.bim")"
n_out="$(wc -l < "$OUTDIR/pw_plink.bim")"
echo
echo "  Array A    : $n_a SNPs"
echo "  Array B    : $n_b SNPs"
echo "  Union      : $n_out SNPs"
echo "  Duplicates : $(( n_a + n_b - n_out )) SNPs shared between arrays"
echo

echo "================================================================================"
echo "  STEP 0 DONE."
echo "  Combined PLINK: $OUTDIR/pw_plink.{bed,bim,fam}"
echo "  Next: Rscript aadr/02_build_aadr_subset.R <aadr_prefix> ./$OUTDIR"
echo "================================================================================"
