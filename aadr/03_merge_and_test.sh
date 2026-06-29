#!/usr/bin/env bash
# ==============================================================================
# 03_merge_and_test.sh -- Germany_IA artefact test, step 3/4
#
# Merges your converted personal genome (pw_plink, from step 1) with the AADR
# population subset (aadr_subset, from step 2), fixes the FID/population
# labelling, and produces a final PLINK fileset ready for extract_f2()/qpAdm
# in R (see 04_run_qpadm_test.R).
#
# IDEMPOTENCY: if merged_pw_aadr_final_ready.{bed,bim,fam} already exist,
# all merge steps are skipped -- safe to re-run after a partial failure.
#
# PITFALL FIXES (all discovered the hard way):
#   1. plink2's --pmerge / --pmerge-list cannot handle "non-concatenating"
#      merges (samples with largely non-overlapping SNP sets) in this plink2
#      version -- it errors "Non-concatenating --pmerge[-list] is under
#      development." FIX: use plink1.9 (binary name "p-link" on this system,
#      NOT "plink") for the actual --bmerge, which has handled this case for
#      a decade.
#   2. plink1.9's --sort-vars requirement: --make-bed cannot sort directly
#      ("Fixed-width .bed/.pgen output doesn't support sorting yet"). FIX:
#      sort via --make-pgen first, then convert that sorted pgen to bed.
#   3. SNPs with A1=A2 in the .bim (all-missing SNPs from step 2, or any
#      other source) make plink1.9 abort the merge with "Identical A1 and A2
#      alleles". FIX: exclude these BEFORE merging.
#   4. Differing platforms (AADR capture vs consumer array) produce a few %
#      of SNPs where one side has a real allele and the other has "." --
#      plink1.9 sometimes mis-flags these as "3+ alleles present" even though
#      it's not a true strand conflict. --flip does NOT fix this (verified:
#      identical error count before/after). FIX: exclude the conflicting SNP
#      IDs from BOTH filesets before merging.
#   5. [RESOLVED in p-link v1.90b7.11.d (2026-06-06)] The 2014 build
#      (v1.90b2p) reliably segfaulted at ~99% during the FINAL --make-bed write
#      after --bmerge, even though the intermediate <out>-merge.* was correct.
#      The 2026 build from the current plink-ng master does NOT exhibit this
#      bug (tested on 1,046,397 variants x 342 samples). The intermediate
#      fileset (-merge.*) is still used as a fallback if the main --make-bed
#      fails, but no longer expected to be needed.
#   6. extract_f2()'s auto_only=TRUE chromosome filter can fail to parse a
#      mixed numeric/X/Y chromosome column. FIX: strip non-autosomal variants
#      via plink2 --chr 1-22 BEFORE calling extract_f2().
#   7. After any merge, .fam FID defaults to sample ID, not population --
#      extract_f2(pops=...) matches against FID, so without a fix every
#      population-based query fails with "Populations missing in indfile".
#      FIX: rewrite FID from the population lookup table written in step 2,
#      keeping your own sample's FID as its own ID.
#   8. aadr_subset_pop_lookup.tsv may be absent if the AADR subset was built
#      manually (not via 02_build_aadr_subset.R). FIX: set AADR_PREFIX env
#      var to the AADR .ind prefix so the lookup can be reconstructed.
#
# Usage:
#   ./03_merge_and_test.sh <run_outdir> [your_sample_id]
#
#   AADR_PREFIX env var (optional):
#     If aadr_subset_pop_lookup.tsv is missing, set this to the AADR .ind
#     prefix (e.g. /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB) and the
#     lookup will be reconstructed from the .ind file automatically.
#
#   <run_outdir>     the SAME directory steps 1-2 wrote to.
#                    Must contain pw_plink.*, aadr_subset.*, and
#                    aadr_subset_pop_lookup.tsv (or AADR_PREFIX must be set).
#   [your_sample_id] the sample ID used in step 1 (default: PW)
#
# Output (all inside <run_outdir>):
#   merged_pw_aadr_final_ready.{bed,bim,fam}  -- ready for R / qpAdm
# ==============================================================================
set -euo pipefail

RUN_DIR="${1:?Usage: $0 <run_outdir> [sample_id]}"
SAMPLE_ID="${2:-PW}"
AADR_PREFIX="${AADR_PREFIX:-}"

cd "$RUN_DIR"

require_file() {
  [[ -f "$1" ]] || { echo "ERROR: required file not found: $1" >&2; exit 1; }
}

echo "================================================================================"
echo "  STEP 3: Merge personal genome with AADR subset"
echo "================================================================================"

# ---- Idempotency: skip everything if final output already exists ----
if [[ -f merged_pw_aadr_final_ready.bed ]]; then
  echo
  echo "  merged_pw_aadr_final_ready.{bed,bim,fam} already exist -- skipping merge."
  echo "  Delete them to force a re-run from scratch."
  echo "================================================================================"
  echo "  STEP 3 ALREADY DONE."
  echo "  Final fileset: $(pwd)/merged_pw_aadr_final_ready.{bed,bim,fam}"
  echo "================================================================================"
  exit 0
fi

require_file pw_plink.bed
require_file aadr_subset.map
require_file aadr_subset.ped

# ---- Pop lookup: reconstruct from AADR .ind if missing ----
if [[ ! -f aadr_subset_pop_lookup.tsv ]]; then
  if [[ -n "${AADR_PREFIX}" && -f "${AADR_PREFIX}.ind" ]]; then
    echo
    echo "  aadr_subset_pop_lookup.tsv missing -- reconstructing from ${AADR_PREFIX}.ind"
    { printf "ID\tPop\n"; awk 'NR==FNR{iid[$1]=1;next} $1 in iid{print $1"\t"$3}' \
        aadr_subset.fam "${AADR_PREFIX}.ind"; } > aadr_subset_pop_lookup.tsv
    echo "  Wrote $(( $(wc -l < aadr_subset_pop_lookup.tsv) - 1 )) ID→population mappings."
  else
    echo "ERROR: aadr_subset_pop_lookup.tsv not found." >&2
    echo "       Set AADR_PREFIX=/path/to/v66.p1_HO.aadr.patch.PUB to auto-reconstruct." >&2
    exit 1
  fi
fi

# ---- 3a. Convert aadr_subset PED/MAP -> binary PLINK ----
echo
echo "[3a] Converting aadr_subset.{map,ped} -> binary PLINK..."
plink2 --pedmap aadr_subset --make-bed --out aadr_subset_raw

# ---- 3b. Exclude all-missing SNPs ----
echo
echo "[3b] Excluding SNPs with identical A1/A2 (all-missing in this subset)..."
awk '$5 == $6 {print $2}' aadr_subset_raw.bim > bad_snps_aadr.txt
n_bad="$(wc -l < bad_snps_aadr.txt)"
echo "  Found $n_bad such SNPs."
if [[ "$n_bad" -gt 0 ]]; then
  plink2 --bfile aadr_subset_raw --exclude bad_snps_aadr.txt --make-bed --out aadr_subset_clean
else
  cp aadr_subset_raw.bed aadr_subset_clean.bed
  cp aadr_subset_raw.bim aadr_subset_clean.bim
  cp aadr_subset_raw.fam aadr_subset_clean.fam
fi

awk '$5 == $6 {print $2}' pw_plink.bim > bad_snps_pw.txt
n_bad_pw="$(wc -l < bad_snps_pw.txt)"
echo "  Found $n_bad_pw such SNPs in your own sample."
if [[ "$n_bad_pw" -gt 0 ]]; then
  plink2 --bfile pw_plink --exclude bad_snps_pw.txt --make-bed --out pw_plink_clean
else
  cp pw_plink.bed pw_plink_clean.bed
  cp pw_plink.bim pw_plink_clean.bim
  cp pw_plink.fam pw_plink_clean.fam
fi

# ---- 3c. Sort both filesets (required by --bmerge) ----
echo
echo "[3c] Sorting variants (--sort-vars needs pgen->bed roundtrip in plink2)..."
plink2 --bfile aadr_subset_clean --sort-vars --make-pgen --out aadr_subset_pgen
plink2 --pfile aadr_subset_pgen --make-bed --out aadr_subset_sorted

plink2 --bfile pw_plink_clean --sort-vars --make-pgen --out pw_plink_pgen
plink2 --pfile pw_plink_pgen --make-bed --out pw_plink_sorted

# ---- 3d. First merge attempt to discover 3+-allele conflicts ----
echo
echo "[3d] First merge attempt to discover cross-platform allele conflicts..."
echo "     (expected to fail if any exist -- --flip does NOT resolve these)"
set +e
p-link --bfile aadr_subset_sorted \
       --bmerge pw_plink_sorted.bed pw_plink_sorted.bim pw_plink_sorted.fam \
       --make-bed --out merge_probe 2>&1 | tee merge_probe_attempt.log
set -e

if [[ -f merge_probe-merge.missnp ]]; then
  n_conflict="$(wc -l < merge_probe-merge.missnp)"
  echo "  Found $n_conflict cross-platform allele conflicts. Excluding from both sides."
  plink2 --bfile aadr_subset_sorted --exclude merge_probe-merge.missnp \
         --make-bed --out aadr_subset_final
  plink2 --bfile pw_plink_sorted --exclude merge_probe-merge.missnp \
         --make-bed --out pw_plink_final
else
  echo "  No conflicts reported -- proceeding with sorted filesets directly."
  cp aadr_subset_sorted.bed aadr_subset_final.bed
  cp aadr_subset_sorted.bim aadr_subset_final.bim
  cp aadr_subset_sorted.fam aadr_subset_final.fam
  cp pw_plink_sorted.bed pw_plink_final.bed
  cp pw_plink_sorted.bim pw_plink_final.bim
  cp pw_plink_sorted.fam pw_plink_final.fam
fi

# ---- 3e. Real merge via plink1.9 ----
# p-link v1.90b7.11.d (2026-06-06) no longer segfaults at --make-bed after
# --bmerge (tested on 1,046,397 variants x 342 samples). If --make-bed
# fails for any reason but the intermediate <out>-merge.* was written
# (behaviour of the 2014 build), we fall back to that intermediate.
echo
echo "[3e] Merging via p-link 1.9 (handles non-concatenating merges plink2 cannot)..."
set +e
p-link --bfile aadr_subset_final \
       --bmerge pw_plink_final.bed pw_plink_final.bim pw_plink_final.fam \
       --make-bed --out merged_pw_aadr 2>&1 | tee merge_final_attempt.log
merge_exit=$?
set -e

if [[ $merge_exit -eq 0 && -f merged_pw_aadr.bim ]]; then
  echo "  Merge completed successfully."
  plink2 --bfile merged_pw_aadr --make-bed --out merged_pw_aadr_final
elif [[ -f merged_pw_aadr-merge.bim ]]; then
  echo "  --make-bed exit $merge_exit but intermediate merge fileset found."
  echo "  Using merged_pw_aadr-merge.* (correct merge result) via plink2..."
  plink2 --bfile merged_pw_aadr-merge --make-bed --out merged_pw_aadr_final
else
  echo "ERROR: merge failed and no intermediate fileset found." >&2
  echo "       Check merge_final_attempt.log for the actual plink1.9 error." >&2
  exit 1
fi

# ---- 3f. Strip non-autosomal variants ----
echo
echo "[3f] Restricting to chromosomes 1-22..."
plink2 --bfile merged_pw_aadr_final --chr 1-22 --make-bed --out merged_pw_aadr_auto

# ---- 3g. Fix FID to carry population labels ----
echo
echo "[3g] Restoring population labels as FID (required for extract_f2(pops=))..."
python3 - "$SAMPLE_ID" <<'PYEOF'
import sys

sample_id = sys.argv[1]

pop_map = {}
with open("aadr_subset_pop_lookup.tsv") as f:
    next(f)  # header
    for line in f:
        iid, pop = line.rstrip("\n").split("\t")
        pop_map[iid] = pop

fixed_rows = []
with open("merged_pw_aadr_auto.fam") as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        iid = fields[1]
        if iid == sample_id:
            fields[0] = sample_id
        elif iid in pop_map:
            fields[0] = pop_map[iid]
        fixed_rows.append(fields)

with open("merged_pw_aadr_auto.fam", "w") as f:
    for fields in fixed_rows:
        f.write("\t".join(fields) + "\n")

print(f"  Rewrote FID for {len(fixed_rows)} samples using population lookup table.")
PYEOF

echo
echo "  Verification:"
grep -w "$SAMPLE_ID" merged_pw_aadr_auto.fam | head -2
echo "  ..."
awk '{print $1}' merged_pw_aadr_auto.fam | sort -u | wc -l
echo "  unique population (FID) labels found above."

mv merged_pw_aadr_auto.bed merged_pw_aadr_final_ready.bed
mv merged_pw_aadr_auto.bim merged_pw_aadr_final_ready.bim
mv merged_pw_aadr_auto.fam merged_pw_aadr_final_ready.fam

echo
echo "================================================================================"
echo "  STEP 3 DONE."
echo "  Final fileset: $(pwd)/merged_pw_aadr_final_ready.{bed,bim,fam}"
echo "  Next: run 04_run_qpadm_test.R against this prefix."
echo "================================================================================"
