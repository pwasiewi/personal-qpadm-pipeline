# AADR / qpAdm Ancestry Analysis Pipeline (EU variant) — Beginner's Guide

> **`aadr_eu` — West/Central-European-expanded variant of `aadr/`.**
> Same preparation and f-statistics machinery as `../aadr/`, but it runs only the
> shared prep (steps **00–03**) plus steps **07** (1-source) and **08** (2-source
> rotating scan) — the two steps that carry an additional batch of **West/Central
> European** reference populations (British / Anglo-Saxon, continental Germanic,
> Gaulish/French, Italic/Roman, Iberian + modern English/French/Italian/Spanish
> anchors). All Polish/Slavic/Baltic/steppe populations from `aadr/` are
> **retained** — the EU pops are *added*, not substituted. The Polish-target
> probes (steps 04–06, 09) are **not** duplicated here; use `../aadr/` for those.
>
> Use this variant for **admixed / non-Polish targets** — e.g. US samples from
> the public PGP-HMS / Harvard dataset (`familytreedna_am/n01/nw_autosomal.csv`),
> whose dominant ancestry is Western European with minor Asian/other components.
> The pool is grown **gradually** (each added pop costs RAM + disk during
> `extract_f2`); the n≥3 gate in steps 07/08 silently skips any label absent
> from the installed AADR version.

This pipeline analyses your personal DNA (from a consumer test kit) against
thousands of ancient genomes from the AADR database using **qpAdm** — a
formal statistical method based on f-statistics (allele frequency correlations).
The result tells you which ancient populations are statistically compatible
with your ancestry, and which are not.

This is a different, more rigorous approach than G25/nMonte (geometric
distance). qpAdm works on raw allele data and can distinguish genuine shared
ancestry from geometric coincidence.

---

## 0. Scientific purpose

This repository contains an experimental research pipeline for analysing
consumer autosomal DNA using qpAdm and the Allen Ancient DNA Resource (AADR).

The project was developed as an independent scientific and educational study
to evaluate the applicability of formal f-statistics to personal ancestry
analysis. It is intended for reproducible research, methodological evaluation,
and learning.

No personal genomic data are included in this repository. Users must analyse
their own legally obtained genotype files (e.g. FamilyTreeDNA, MyHeritage or
AncestryDNA) and must obtain the AADR reference dataset separately according
to its licence terms.

The pipeline is provided for research purposes only. The results should not be
interpreted as medical, forensic, or genealogical proof without independent
validation.

---

## 1. Buy and download your raw DNA data

You need a **raw DNA file** (CSV format) from one of these providers:

### Family Tree DNA — Family Finder
1. Order at **familytreedna.com** → "Family Finder" autosomal test
2. After results arrive: **myFTDNA** → top-right menu → **"Download Raw Data"**
   → "Chromosome Browser Raw Data" → download the zip
3. Unzip: `unzip <filename>.zip` → you get a `.csv` file

### MyHeritage DNA
1. Order at **myheritage.com/dna**
2. After results arrive: **DNA** tab → **"Manage DNA kits"** → your kit →
   **"Download"** → "Download raw DNA data" → confirm by email link
3. Unzip: `unzip MyHeritage_raw_dna_data.zip` → you get a `.csv` file

### Ancestry DNA  *(partial support — step 01 auto-detects)*
Download via: **DNA** → **"Settings"** → **"Download DNA Raw Data"**

> **Tip:** If you have results from multiple providers for the same person,
> run step 01 separately for each, then merge with step 00
> (`00_merge_arrays.sh`) before continuing.

---

## 2. Prerequisites — software and data

Everything below must already be installed. On Gentoo with overlay `pwr`:

```bash
# R packages
Rscript -e 'remotes::install_github("uqrmaie1/admixtools")'

# Tools
sudo emerge sci-biology/plink2-bin   # binary: plink2
sudo emerge sci-biology/plink        # binary: p-link  (v1.90b7+, overlay pwr)

# AADR v66.p1 HO dataset — ~30 GB download, stored read-only
# Download from Harvard Dataverse: doi:10.7910/DVN/FFIDCW
# Place at: /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB.{geno,snp,ind,anno}
ls /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB.geno   # should exist
```

---

## 3. Run the pipeline

**Replace `NN` with your own short identifier throughout** (2–4 letters,
no spaces, e.g. `AB`, `ANNA`, `JOHN`). This becomes your sample label in all
output tables.

Run all commands from `~/Claude` (or wherever this repo lives):

```bash
cd ~/Claude
```

### Step 01 — Convert raw CSV to PLINK format

```bash
# FTDNA Family Finder:
SAMPLE_ID=NN ./aadr_eu/01_prepare_files.sh /path/to/your_ftdna_rawdata.csv NN_FF

# MyHeritage:
SAMPLE_ID=NN ./aadr_eu/01_prepare_files.sh /path/to/MyHeritage_raw_dna_data.csv NN_MH
```

- `SAMPLE_ID=NN` — sets your label (used in .fam file and all later outputs)
- Second argument (`NN_FF` / `NN_MH`) — name of the output directory
- Output: `./NN_FF/pw_plink.{bed,bim,fam}` — your genome in PLINK format

### Step 02 — Build AADR population subset

```bash
Rscript aadr_eu/02_build_aadr_subset.R \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
    ./NN_FF
```

- Reads AADR, extracts the ~23 reference populations needed for analysis
- Memory-safe: does NOT load the full 129 GB panel — uses `pops=` filter
- Output: `./NN_FF/aadr_subset.{bed,bim,fam}` + `aadr_subset_pop_lookup.tsv`
- Runtime: ~5–10 minutes

### Step 03 — Merge personal genome with AADR subset

```bash
./aadr_eu/03_merge_and_test.sh ./NN_FF NN
```

- Second argument must match your `SAMPLE_ID` from step 01
- Merges your ~520K SNPs with AADR HO panel → ~41–50K overlapping SNPs
  *(this is expected and normal — not a bug)*
- Output: `./NN_FF/merged_pw_aadr_final_ready.{bed,bim,fam}`
- Runtime: ~5–15 minutes

> **Steps 04–06 and 09 are not part of the EU variant.** They are Polish-target
> probes — the Germany_IA artefact test (04–05), Poland-only ancestry models (06),
> and the Przeworsk/Wielbark Iron Age channel test (09). For an admixed / non-Polish
> target they only re-confirm closed Polish findings, so the EU pipeline goes
> straight from **step 03 to steps 07 and 08**. If you do want those Polish-specific
> analyses, run the canonical copies in **`../aadr/`** (steps `04`–`06`, `09`).

### Step 07 — Extended single-source scan (EU-expanded)

```bash
Rscript aadr_eu/07_extended_models.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Tests single-source models across a wide pool: Slavic, Baltic, Lusatian,
  steppe (Avar/Hun/Sarmatian/Bulgar/Xiongnu) **plus the EU-expanded batch** —
  British/Anglo-Saxon, continental Germanic, Gaulish/French, Italic/Roman,
  Iberian + modern English/French/Italian/Spanish/Orcadian/Icelandic anchors
- Identifies the best single-source ancient match for your ancestry
- Reference result (Polish targets): **Ukraine_IA_Lusatian** is the best fit
  (p ≈ 0.86–0.96). For an admixed EU target, expect a broad spread of
  non-rejected European sources (the consumer-array ceiling — see CLAUDE.md)
- Output TSV: `./NN_FF/extended_models_results.tsv`

### Step 08 — Two-source rotating scan (hidden components)

```bash
Rscript aadr_eu/08_two_source.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Rotates a panel of well-powered base anchors against a wide pool of second
  sources (Uralic, Siberian, Steppe, Hunnic, Middle Eastern, Caucasian +
  the EU-expanded batch) to detect a hidden second ancestry component
- Verdict `SIGNIFICANT-2src` = real additional component found. The detector is
  deliberately strict (hardened 2026-06-29): a base must have **n ≥ 15**; the
  second source must **improve** on the base-alone fit, carry a non-vestigial
  weight (both weights in 0.1–0.9), reach **|z| > 3**, and survive a
  Benjamini-Hochberg correction (`q_second < 0.05`) across all ~400 models. This
  rejects weak-anchor artefacts (a low-n base floats the whole model onto the
  second source and manufactures fake "signals")
- Verdict `UNIDENTIFIABLE` = model numerically unstable (not a finding)
- Reference result: **0 significant two-source models** for either the Polish
  targets or the admixed US target `AM` → no hidden Avar/steppe/Asian/EU
  minority component resolvable at the consumer-array SNP ceiling
- Output TSV: `./NN_FF/two_source_results_NN.tsv` (column `q_second` = the
  multiple-testing-corrected significance)

---

## 4. Full command sequence at a glance

```bash
cd ~/Claude

# -- replace NN with your label, and update the CSV path --
SAMPLE_ID=NN ./aadr_eu/01_prepare_files.sh /path/to/rawdata.csv NN_FF

Rscript aadr_eu/02_build_aadr_subset.R \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB ./NN_FF

./aadr_eu/03_merge_and_test.sh ./NN_FF NN

Rscript aadr_eu/07_extended_models.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN

Rscript aadr_eu/08_two_source.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

> Steps 04–06 and 09 (Polish-target probes) are not part of this variant — see
> `../aadr/` if you need them.

### Worked example — US public-dataset sample (`AM`)

For an admixed US genome from the public Harvard/PGP-HMS data, the EU-expanded
candidate pools in steps 07–08 are exactly what's needed — run **01–03 then
07–08**:

```bash
cd ~/Claude
AADR=/usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB

SAMPLE_ID=AM ./aadr_eu/01_prepare_files.sh \
    /usr/local/share/aadr/familytreedna_am_autosomal.csv AM_FF
Rscript aadr_eu/02_build_aadr_subset.R  "$AADR" ./AM_FF
./aadr_eu/03_merge_and_test.sh          ./AM_FF AM
Rscript aadr_eu/07_extended_models.R    ./AM_FF/merged_pw_aadr_final_ready  "$AADR" AM
Rscript aadr_eu/08_two_source.R         ./AM_FF/merged_pw_aadr_final_ready  "$AADR" AM
```

The other two public genomes use the same recipe — swap `am`→`n01`/`nw` and the
label `AM`→`N01`/`NW` (and the output dir `AM_FF`→`N01_FF`/`NW_FF`). Run them one
at a time: each `extract_f2` pass with the widened EU pool is the RAM/disk peak.

---

## 5. What to look at in the results

| File | What it shows |
|------|---------------|
| `NN_FF/extended_models_results.tsv` | Best single-source match across the full pool; sort by `p` descending |
| `NN_FF/two_source_results_NN.tsv` | Whether a second hidden ancestry component exists (`q_second` = corrected significance) |

**Reading the p-value:** A model is compatible with your data when `p > 0.05`.
Higher p = better fit (the model cannot be statistically rejected). The best
single-source fit typically has `p > 0.5` for Polish/Slavic ancestry.

**Reading the verdict column (steps 07–08):**
- `OK-*` — model fits, p > 0.05 (suffix says how the weight splits)
- `REJECTED` — model rejected, p < 0.05
- `SIGNIFICANT-2src` — real second component detected (step 08 only; passes the
  strict n≥15 / improvement / |z|>3 / BH-corrected gate)
- `UNIDENTIFIABLE` — weights numerically unstable; not a meaningful finding

---

## 6. Example data and reference runs

- **Public genomes to try this on** — `examples/site.txt` points to the PGP-HMS
  open repository (<https://my.pgp-hms.org/public_genetic_data>), where you can
  download consenting individuals' raw autosomal files from various testing
  companies. Worked-example outputs for these will be added under `examples/`.
- **Canonical Polish reference runs** (the methodology baseline) live with the
  Polish pipeline: `../aadr/PW_FF_slavic_aadr.md` (FTDNA, western Mazovia) and
  `../aadr/JW_MH_slavic_aadr.md` (MyHeritage, Podlasie).

---

## 7. Expected runtimes and SNP counts

| Step | Runtime | Notes |
|------|---------|-------|
| 01 | < 1 min | CSV conversion |
| 02 | 5–10 min | AADR subset extraction |
| 03 | 5–15 min | Merge; yields ~41–50K SNPs (normal for consumer arrays) |
| 07 | 10–20 min | f2cache built on first run, reused on re-run |
| 08 | 20–40 min | ~400 two-source models; f2 pass split into 6 chunk-pairs (RAM-capped) |

**~41–50K SNPs after merge is correct and expected.** Consumer arrays cover
~520–560K positions; AADR HO has 584K; physical overlap ~57K; after
`maxmiss=0` across all ancient samples → 41–50K. Not a bug.

---

## 8. Developer reference

Pitfall documentation and conventions for this variant are in `aadr_eu/CLAUDE.md`.
For the full original command history and the platform-merge workflow
(FTDNA + MyHeritage union), see the Polish pipeline in `../aadr/`.

## 9. Acknowledgements

Parts of this project were developed with the assistance of Claude Code
(Anthropic) as an AI programming assistant. All methodological decisions,
implementation, testing, and validation were performed by the repository
author.
