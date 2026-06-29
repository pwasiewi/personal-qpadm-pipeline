#!/usr/bin/env bash
# ==============================================================================
# 01_prepare_files.sh -- Germany_IA artefact test, step 1/4
#
# Converts a raw consumer DNA CSV (FTDNA or MyHeritage format) to PLINK.
#
# OUTPUT LOCATION: a run-specific directory derived from the CSV filename,
# created relative to wherever you invoke this script from -- NOT inside the
# (large, shared, read-only-in-spirit) AADR data directory. This keeps each
# person's conversion run self-contained and lets you run multiple people's
# files against the same AADR install without them colliding.
#
#   Example: running from /home/pwas with
#somewhere/MyHeritage_PW_raw_dna_data.csv
#   creates:  /home/pwas/MyHeritage_PW/  (outdir = CSV basename, minus
#             extension and common raw-data suffixes, under CWD)
#
# KNOWN PITFALLS THIS SCRIPT AVOIDS (discovered the hard way in the original
# session -- see comments inline at each fix):
#   - plink2 needs --pedmap (not bare --file) to read PED/MAP
#   - plink2 cannot read AADR's native packedancestrymap format at all
#   - sample-major .bed write from --pedmap can choke on indel/missing alleles
#
# Usage:
#   ./01_prepare_files.sh /path/to/raw_myheritage_or_ftdna.csv [outdir]
#
#   [outdir]  optional override. Default: derived from CSV filename, under CWD.
#
# Output:
#   <outdir>/pw_plink.{bed,bim,fam}   -- your sample, converted, ready for step 3
#   <outdir>/pw.map, pw.ped           -- intermediate PED/MAP (kept for inspection)
# ==============================================================================
set -euo pipefail

RAW_CSV="${1:?Usage: $0 <raw_dna_csv> [outdir]}"
RAW_CSV="$(realpath "$RAW_CSV")"
SAMPLE_ID="${SAMPLE_ID:-PW}"   # override with: SAMPLE_ID=yourname ./01_...

# ---- Derive OUTDIR from CSV filename, created relative to CWD ----
# Strips extension and common raw-export suffixes so
# "MyHeritage_PW_raw_dna_data.csv" -> outdir "MyHeritage_PW", and
# "FF_PW_Chrom_Autoso_20260223.csv" -> outdir
# "FF_PW_Chrom_Autoso_20260223" (no generic suffix to strip there,
# left as-is -- pass an explicit [outdir] argument if you want a shorter name).
_basename="$(basename "$RAW_CSV")"
_stem="${_basename%.*}"
_stem="${_stem%_raw_dna_data}"
_stem="${_stem%_raw}"

OUTDIR="${2:-$(pwd)/$_stem}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "================================================================================"
echo "  STEP 1: Convert raw DNA CSV -> PLINK"
echo "================================================================================"
echo "  Input CSV : $RAW_CSV"
echo "  Output dir: $OUTDIR"
echo

# ---- Detect format (FTDNA vs MyHeritage) by header inspection ----
# FTDNA Family Finder: header "RSID,CHROMOSOME,POSITION,RESULT", no build comment.
# MyHeritage: leading "##" comment block including "##reference=buildNN".
BUILD_LINE="$(grep -m1 '^##reference=' "$RAW_CSV" 2>/dev/null || true)"
if [[ -n "$BUILD_LINE" ]]; then
  echo "  Detected MyHeritage-style file. $BUILD_LINE"
  if [[ "$BUILD_LINE" != *build37* ]]; then
    echo "  WARNING: build is not build37/hg19 -- AADR is hg19. A liftover would"
    echo "           be needed before merging; this script does not do that."
  fi
else
  echo "  No '##reference=' header found -- assuming FTDNA-style file (build hg19"
  echo "  typical for Family Finder, but NOT verified by this script)."
fi

# ---- Convert CSV -> PED/MAP ----
# Handles both formats: lines may or may not be quoted, may or may not have a
# leading "#" comment block. Skips rows with chromosome 0/X/Y/MT or position 0
# (matches AADR's auto_only=TRUE filtering downstream, and avoids "0,0,--"
# placeholder rows seen in some FTDNA exports for unmapped probes).
python3 - "$RAW_CSV" "$SAMPLE_ID" <<'PYEOF'
import csv, sys

raw_path, sample_id = sys.argv[1], sys.argv[2]

snps = []
with open(raw_path, encoding="utf-8-sig") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        if line.upper().startswith("RSID"):
            continue
        # Strip quotes (MyHeritage quotes every field; FTDNA does not)
        parts = [p.strip('"') for p in line.split(",")]
        if len(parts) != 4:
            continue
        rsid, chrom, pos, result = parts
        if chrom in ("0", "X", "Y", "MT", "0.0") or pos in ("0", "0.0", ""):
            continue
        if not rsid.startswith("rs"):
            # Skip non-standard probe IDs (e.g. FTDNA "2010-08-Y-1221",
            # "GSA-1:115030388") -- these have no AADR-comparable rsID anyway.
            continue
        # Normalise missing/indel codes to PLINK's "0" missing-allele symbol
        result = result.replace("-", "0").replace("I", "0").replace("D", "0")
        a1 = result[0] if len(result) > 0 else "0"
        a2 = result[1] if len(result) > 1 else "0"
        snps.append((rsid, chrom, pos, a1, a2))

with open("pw.map", "w") as f:
    for rsid, chrom, pos, _, _ in snps:
        f.write(f"{chrom}\t{rsid}\t0\t{pos}\n")

with open("pw.ped", "w") as f:
    alleles = "\t".join(f"{a1}\t{a2}" for _, _, _, a1, a2 in snps)
    f.write(f"{sample_id}\t{sample_id}\t0\t0\t0\t-9\t{alleles}\n")

print(f"  Wrote pw.map / pw.ped ({len(snps)} SNPs, autosomes only, sample={sample_id})")
PYEOF

# ---- PED/MAP -> binary PLINK ----
# PITFALL FIX: plink2 needs --pedmap explicitly; bare --file (plink1.9 syntax)
# is not recognised ("unknown option").
plink2 --pedmap pw --make-bed --out pw_plink

echo
echo "  Sample genotyping rate check (should be 0 missing -- this is YOUR own"
echo "  file before any merge; non-zero here means the CSV parser above has a"
echo "  bug, not a downstream merge issue):"
plink2 --bfile pw_plink --missing --out pw_plink_selfcheck 2>&1 | tail -3
cat pw_plink_selfcheck.smiss

echo
echo "================================================================================"
echo "  STEP 1 DONE. Output: $OUTDIR/pw_plink.{bed,bim,fam}"
echo "  Next: run 02_build_aadr_subset.R, pointing its outdir at:"
echo "    $OUTDIR"
echo "================================================================================"
