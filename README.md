# AADR / qpAdm Ancestry Analysis Pipeline — Beginner's Guide

> **Which variant do I want?** This `aadr/` pipeline is tuned for **Poland and
> its neighbours** (Slavic / Baltic / Central-East European ancestry) — its
> candidate populations are Polish, Slavic, Baltic, Lusatian and steppe sources.
> For **Western/Central European or US (admixed) genomes**, use the sibling
> `../aadr_eu/` instead: it adds British/Anglo-Saxon, continental Germanic,
> Gaulish/French, Italic/Roman and Iberian reference populations. Note that
> `aadr_eu/` is at an **early stage** (only the prep steps 00–03 and the
> EU-expanded scans 07–08), whereas `aadr/` is the mature, full pipeline.

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
> run step 01 separately for each, then merge with step 00 before continuing.
> See `README_commands.md` for the merge workflow.

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
SAMPLE_ID=NN ./aadr/01_prepare_files.sh /path/to/your_ftdna_rawdata.csv NN_FF

# MyHeritage:
SAMPLE_ID=NN ./aadr/01_prepare_files.sh /path/to/MyHeritage_raw_dna_data.csv NN_MH
```

- `SAMPLE_ID=NN` — sets your label (used in .fam file and all later outputs)
- Second argument (`NN_FF` / `NN_MH`) — name of the output directory
- Output: `./NN_FF/pw_plink.{bed,bim,fam}` — your genome in PLINK format

### Step 02 — Build AADR population subset

```bash
Rscript aadr/02_build_aadr_subset.R \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
    ./NN_FF
```

- Reads AADR, extracts the ~23 reference populations needed for analysis
- Memory-safe: does NOT load the full 129 GB panel — uses `pops=` filter
- Output: `./NN_FF/aadr_subset.{bed,bim,fam}` + `aadr_subset_pop_lookup.tsv`
- Runtime: ~5–10 minutes

### Step 03 — Merge personal genome with AADR subset

```bash
./aadr/03_merge_and_test.sh ./NN_FF NN
```

- Second argument must match your `SAMPLE_ID` from step 01
- Merges your ~520K SNPs with AADR HO panel → ~41–50K overlapping SNPs
  *(this is expected and normal — not a bug)*
- Output: `./NN_FF/merged_pw_aadr_final_ready.{bed,bim,fam}`
- Runtime: ~5–15 minutes

### Step 04 — Germany_IA artefact test

```bash
Rscript aadr/04_run_qpadm_test.R \
    ./NN_FF/merged_pw_aadr_final_ready NN
```

- Tests whether Germany_IA appears as a genuine ancestry source or a proxy artefact
- Established result for Slavic/Polish ancestry: Germany_IA is **always an artefact**
  (shared Corded Ware substrate, not Iron Age Germanic ancestry)

### Step 05 — Robustness check across outgroup combinations

```bash
Rscript aadr/05_robustness_outgroups.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Repeats the Germany_IA test with 6 different outgroup combinations
- Confirms the artefact conclusion is stable, not outgroup-dependent

### Step 06 — Ancestry models (Polish reference populations)

```bash
Rscript aadr/06_ancestry_models.R \
    ./NN_FF/merged_pw_aadr_extended \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Tests 16 models: single-source and two/three-source combinations
  using Polish medieval and Iron Age populations
- Output TSV: `./NN_FF/06_ancestry_models_NN.tsv`

### Step 07 — Extended models (non-Polish references)

```bash
Rscript aadr/07_extended_models.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Tests 23 single-source models: Slavic, Baltic, Lusatian culture,
  Pannonian Avars, Huns, Sarmatians, Steppe Bulgars, Xiongnu
- Identifies the best single-source ancient match for your ancestry
- Established result: **Ukraine_IA_Lusatian** is the best fit for
  Polish/Slavic ancestry (p ≈ 0.86–0.96)
- Output TSV: `./NN_FF/07_extended_models_NN.tsv`

### Step 08 — Two-source rotating scan (hidden components)

```bash
Rscript aadr/08_two_source.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Takes the best single-source fit and tests 23 diverse second sources
  (Uralic, Siberian, Steppe, Hunnic, Middle Eastern, Caucasian)
- Detects whether a hidden ancestry component exists beyond the Slavic baseline
- Verdict `SIGNIFICANT-2src` = real additional component found
- Verdict `UNIDENTIFIABLE` = model numerically unstable (not a finding)
- Established result: **0 significant two-source models** for Polish ancestry
  → no hidden Avar/steppe/Asian admixture
- Output TSV: `./NN_FF/08_two_source_NN.tsv`

### Step 09 — Przeworsk vs Wielbark channel test *(Polish Iron Age, optional)*

```bash
Rscript aadr/09_przeworsk_channel.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

- Tests whether your ancestry passed through the **Przeworsk** or the
  **Wielbark** Iron Age culture channel (relevant for a Greater-Poland /
  Wielkopolska genealogical background)
- Combines three probes: 1-source p-value per channel, an **f4 symmetry test**
  (B1 — are the two channels distinguishable at all? B2 — which channel does
  your genome favour?), and a 2-source model on a Lusatian base
- Admits `Poland_IA_Przeworsk` (n=2, skipped by steps 6–7's n≥3 gate)
- Established result: **UNRESOLVED but consistent** — the two channels are
  statistically indistinguishable at ~40K SNP (f4 B1 non-significant), both
  targets show symmetric affinity (f4 B2 non-significant); Przeworsk has the
  highest single-source p of the Polish IA sources. Needs WGS to resolve.
- Output TSV: `./NN_FF/przeworsk_channel_results_NN.tsv`,
  `./NN_FF/przeworsk_channel_f4_NN.tsv`

---

## 4. Full command sequence at a glance

```bash
cd ~/Claude

# -- replace NN with your label, and update the CSV path --
SAMPLE_ID=NN ./aadr/01_prepare_files.sh /path/to/rawdata.csv NN_FF

Rscript aadr/02_build_aadr_subset.R \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB ./NN_FF

./aadr/03_merge_and_test.sh ./NN_FF NN

Rscript aadr/04_run_qpadm_test.R \
    ./NN_FF/merged_pw_aadr_final_ready NN

Rscript aadr/05_robustness_outgroups.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN

Rscript aadr/06_ancestry_models.R \
    ./NN_FF/merged_pw_aadr_extended \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN

Rscript aadr/07_extended_models.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN

Rscript aadr/08_two_source.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN

# Step 09 is optional (Polish Iron Age channel test):
Rscript aadr/09_przeworsk_channel.R \
    ./NN_FF/merged_pw_aadr_final_ready \
    /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB NN
```

---

## 5. What to look at in the results

| File | What it shows |
|------|---------------|
| `NN_FF/06_ancestry_models_NN.tsv` | Best Polish-population models, p-values |
| `NN_FF/07_extended_models_NN.tsv` | Best single-source match across 23 populations; sort by `p` descending |
| `NN_FF/08_two_source_NN.tsv` | Whether a second hidden ancestry component exists |
| `NN_FF/przeworsk_channel_results_NN.tsv` | Przeworsk vs Wielbark Iron Age channel (step 09, optional) |
| `NN_FF/przeworsk_channel_f4_NN.tsv` | f4 symmetry stats for the channel test (step 09) |

**Reading the p-value:** A model is compatible with your data when `p > 0.05`.
Higher p = better fit (the model cannot be statistically rejected). The best
single-source fit typically has `p > 0.5` for Polish/Slavic ancestry.

**Reading the verdict column (steps 6–8):**
- `PASS` — model fits, p > 0.05
- `REJECT` — model rejected, p < 0.05
- `SIGNIFICANT-2src` — real second component detected (step 08 only)
- `UNIDENTIFIABLE` — weights numerically unstable; not a meaningful finding

---

## 6. Example results (PW and JW reference runs)

See the canonical analysis reports for two existing runs:

- `aadr/PW_FF_slavic_aadr.md` — FTDNA Family Finder, western Mazovia background
- `aadr/JW_MH_slavic_aadr.md` — MyHeritage, Podlasie/Uniate-Ruthenian background

---

## 7. Expected runtimes and SNP counts

| Step | Runtime | Notes |
|------|---------|-------|
| 01 | < 1 min | CSV conversion |
| 02 | 5–10 min | AADR subset extraction |
| 03 | 5–15 min | Merge; yields ~41–50K SNPs (normal for consumer arrays) |
| 04–05 | 5–10 min each | f2 + qpAdm tests |
| 06–07 | 10–20 min each | f2cache built on first run, reused on re-run |
| 08 | 20–40 min | 46 two-source models |
| 09 | 10–20 min | optional; 1-source + f4 symmetry + 2-source channel test |

**~41–50K SNPs after merge is correct and expected.** Consumer arrays cover
~520–560K positions; AADR HO has 584K; physical overlap ~57K; after
`maxmiss=0` across all ancient samples → 41–50K. Not a bug.

---

## 8. Developer reference

For the full command history, pitfall documentation, and platform-merge
workflow (FTDNA + MyHeritage union), see `aadr/README_commands.md`.

## 9. Acknowledgements

Parts of this project were developed with the assistance of Claude Code
(Anthropic) as an AI programming assistant. All methodological decisions,
implementation, testing, and validation were performed by the repository
author.
