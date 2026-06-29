# AADR / qpAdm Ancestry Analysis Pipeline (EU variant) — Claude Code Project Context

> **`aadr_eu` = West/Central-European-expanded fork of `../aadr/`.** Same code,
> same pitfalls, same invocation. The only difference: steps **07** and **08**
> append a batch of West/Central European reference populations (British /
> Anglo-Saxon, continental Germanic, Gaulish/French, Italic/Roman, Iberian +
> modern English/French/Italian/Spanish/Orcadian/Icelandic anchors). All
> Polish/Slavic/Baltic/Avar/Hun/steppe pops from `aadr/` are **retained**, not
> replaced. Intended for **admixed / non-Polish targets** (e.g. US public
> PGP-HMS genomes `AM`/`N01`/`NW`), assumed ≈ Western European + minor
> Asian/exotic admixture.
>
> **Expansion policy:** grow the pool *gradually* — every added pop costs RAM +
> disk in each `extract_f2` pass. Add a labelled batch, run 07/08, inspect, then
> add the next. The n≥3 gate in 07 / `check_pops()` in 08 silently drops any
> label missing from the installed AADR version, so a wrong/absent label never
> crashes the run.
>
> The closed-findings table below is inherited from `aadr/` and pertains to the
> **Polish targets PW/JW** — it is reference, not a constraint on EU targets.

Personal genome ancestry analysis using formal f-statistics (qpAdm) against
AADR v66.p1 HO. Reference targets carried over from `aadr/`: **PW** (FTDNA
Family Finder, western Mazovia) and **JW** (MyHeritage, Podlasie/Warsaw).
EU-target runs (e.g. `AM`/`N01`/`NW`) live in their own `*_FF/` dirs.
Companion to the G25/nMonte pipeline (`~/Claude/g25/`).

**Canonical Polish-target reports (reference):** `../aadr/PW_FF_slavic_aadr.md`
and `../aadr/JW_MH_slavic_aadr.md` (they live with the Polish pipeline, not here).
Read before re-deriving any Polish closed finding.

## aadr_eu notes — admixed / US targets (2026-06-29)

Worked the first US public-dataset genome (`AM`, FamilyTreeDNA from PGP-HMS).
Findings + fixes that future runs must not re-litigate:

### Continental check first — is the target even European?
For US/admixed samples, run an **admixture-f3** screen before the EU scan:
`qp3pop(f2dir, pop1=TARGET, pop2=European, pop3=X)` for X in
{Yoruba(African), Han(EAsian), Karitiana/Mixe/Mayan(Native American), Papuan(Oceanian)}.
**Negative est & z < −3 ⇒ real admixture with X.** All positive ⇒ European-only.
- **AM result: ALL strongly positive (est ≈ +0.07, z ≈ +50) on every axis** ⇒
  AM is European, **no detectable continental admixture**. So the EU pop scan is
  the right tool; no need for continental left-sources. (Script: `/tmp/cont_extend.R`
  pattern — extract Karitiana/Papuan, merge, f3. ~49K SNPs, moderns genotype well.)
- If a future target shows a negative continental f3, the Europe-only left panel
  **cannot** model it — add continental sources to the LEFT + deep outgroups RIGHT.

### Memory fixes in steps 07/08 (do not regress)
`extract_f2` OOM-killed (`Unicestwiony`) the first AM run. Three changes, all in
07 and 08, all required:
1. **SNP pre-filter** before `extract_f2`: `plink2 --keep-fam TARGET --geno 0
   --write-snplist` → `--extract --chr 1-22 --make-bed` → `*_dense`. The ext7/2src
   union is ~1.1M SNPs but the consumer array genotyped only ~599K; extract_f2
   loaded the full matrix (24–45 GB) before dropping to ~35.7K usable.
2. **`--chr 1-22`** on that `--make-bed` is mandatory: plink2 writes STRING chr
   codes (X/Y/MT) that admixtools `discard_snps(auto_only=TRUE)` cannot parse
   ("Could not parse chromosome numbers" — repo pitfall #11). p-link 1.9 wrote
   numeric codes, which is why the un-filtered path didn't trip it.
3. **`maxmem = 1000`** in `extract_f2` — MUST be below its own "~2515 MB without
   splitting" estimate or it won't split; that no-split f2 pass is what OOM'd.
   (The estimate badly undercounts real peak.) Lower to 500 if a run is killed.

### `extra_needed` must include `right_pops`
Original 07/08 pulled only candidate/left pops from AADR
(`setdiff(confirmed/all_left_pops, existing_fids)`) — assuming outgroups were
always Mbuti/Yoruba/Han already in the base file. Enriching `right_pops` (added
Papuan + Karitiana) silently failed: the new outgroups were never merged →
`extract_f2` "Populations missing in indfile: Papuan, Karitiana". Fixed to
`setdiff(c(..., right_pops), existing_fids)` in both. **Any new right-pop now
forces an ext7/2src rebuild** (delete `AM_FF/merged_pw_aadr_{ext7,2src}*`).

### Step 08 SIGNIFICANT-2src rule hardened (2026-06-29) — weak-base artifact
First AM 08 run flagged **25 SIGNIFICANT-2src**, ALL on base `Ukraine_IA_Lusatian`
(n=3), zero on the 9 well-powered bases — the tell of a weak-anchor artifact, not a
real hidden component (a near-empty base floats the whole model onto the second
source: w_second 0.6–0.93). Three cheap tightenings, all in `08_two_source.R`:
1. **`MIN_BASE_N <- 15L`** — a base must have n≥15. Low-n bases (Ukraine_IA_Lusatian
   n=3, Sweden_IA n=8) are **demoted to the second-source pool**, not dropped.
2. **Improvement + non-vestigial** — `SIGNIFICANT-2src` now also requires `p > p_base`
   AND **both** weights in (0.1, 0.9) (was: only w_second in (0.05, 0.95)).
3. **`|z_second| > 3`** (was >2) **+ Benjamini-Hochberg `q_second < 0.05`** across all
   ~400 models (new `q_second` column). At ~400 tests, |z|>2 ⇒ ~24 chance hits.

Re-run result: **0 SIGNIFICANT-2src** (177 OK-2src-nonsig / 74 UNIDENTIFIABLE / 66
OK-base-suff / 51 OK-second-dom / 33 REJECTED). The clincher: **max |z_second| over
the entire scan on strong bases = 1.381** — nothing even approaches |z|>2, vs 2.0–2.5
on the n=3 base. The 25 were pure weak-anchor inflation. Worth carrying back to
`aadr/08` (PW/JW) — the same lax rule lives there. **Don't re-add low-n pops as bases.**

### Consumer-array ceiling — intra-European is UNRESOLVABLE here
Enriched right set = `Mbuti, Yoruba, Han, Papuan, Karitiana` (moderns, costs ~1
SNP; 35,695 remain). Result for AM: 1-source rejects only **French (0.008)** and
**Orcadian (3e-4)** — the most *specific/drifted* modern pops — while the
steppe/Asian-shifted sources (Hungary_Conqueror 0.57, Russia_Tatarstan_Bulgar
0.80, Avars ~0.65) still pass. Continental outgroups catch only GROSS non-European
shift; a minor Magyar/Avar steppe fraction is invisible at ~35K SNP. **Verdict:
AM = broadly European, not continentally admixed, not attributable to any single
source.** Fine "which corner of Europe" needs WGS or G25/PCA — qpAdm on a consumer
array won't deliver it regardless of pop-list tuning. Same ceiling as the Polish
targets' autochthon/allochthon question. Do not re-tune pop lists hoping to break it.

## Closed findings — do not re-derive, cite instead

All confirmed by qpAdm on both PW and JW independently:

| Question | Answer | Where confirmed |
|---|---|---|
| Germany_IA real ancestry? | **NO** — geometric artefact, 6/6 outgroup combos | steps 4–5 |
| Best single-source fit (38 tested) | **Ukraine_IA_Lusatian p≈0.88–1.00** | step 7 |
| Best Poland_* single-source fit | Poland_EarlySlav p≈0.05–0.83 | step 6 |
| Avar/Pannonian ancestry real? | **NO** — Hungary/Austria_Avar rejected 1-source; 0 significant 2-source | steps 7, 8 |
| Hidden second source (steppe/Asia/ME)? | **NO** — 0 `SIGNIFICANT-2src` models out of 70 (2 bases × 35) | step 8 |
| Ashkenazi Jewish admixture? | **NO** — `Jew_Ashkenazi`/`Germany_Medieval_Jewish` w≤0 or REJECTED, both targets | step 8 (2026-06-29) |
| Magyar / Viking / Wielbark / Saxon hidden? | **NO** — all w≤0 (PW) or non-specific +bias (JW) | step 8 (2026-06-29) |
| PW vs JW ranking qpAdm | Identical (Lusatian #1) — non-individualising at ~41K SNP | steps 6–8 |
| PW vs JW *direction* of difference | JW drifts **SOUTH** (South Slavic/Pannonian; Serbia +0.25, Croatia +0.21), NOT east (Russian/Finnish gain least). PW = tighter Lusatian point | steps 7, 8 (2026-06-29) |
| PW autochthon vs allochthon? | **UNRESOLVED at ~40K SNP** — local Polish chain broken (Roman_Wielbark, Trzciniec REJECTED) but eastern Ukraine_IA_Lusatian best; all sources collinear. Needs WGS | steps 6–8 (2026-06-29) |
| Przeworsk vs Wielbark channel (Wielkopolska background)? | **UNRESOLVED — consistent, not contradicted.** f4(Mbuti,Sweden_IA;Przeworsk,Wielbark) n.s. (channels indistinguishable, both targets); f4(Mbuti,TGT;Wielbark,Przeworsk) n.s. (symmetric). Przeworsk = highest 1-src p of Polish IA (PW 0.351 / JW 0.651). Caveat: Przeworsk n=2 (Gąski, both female). Needs WGS | step 9 (2026-06-29) |
| G25/nMonte weighted artefacts | Hungary_Mig_Himod~49% + Latvia_Neo~28% in PW = bootstrap-stable attractor, falsified by step 8 | step 8 cross-check |

**Germany_IA mechanism:** shared Corded Ware/Bronze Age substrate, not Iron Age
Germanic ancestry. Closed — don't re-run without new AADR label or new target.

**Avar mechanism:** structural collinearity (sources redundant vs Mbuti/Yoruba/Han
outgroups), same class as Germany_IA. More SNPs amplify degeneracy, don't resolve it.

## Repository layout

This fork keeps ONLY the prep steps + the two EU-expanded analysis steps. The
Polish-target probes (04 Germany_IA artefact, 05 robustness, 06 Poland_* models,
09 Przeworsk channel) were byte-identical copies of `../aadr/` with no EU content,
so they were removed (2026-06-29) — run them from `../aadr/` if needed.

```
00_merge_arrays.sh              # union-merge two arrays of same person (MH+FTDNA) before step 01
01_prepare_files.sh             # consumer CSV (FTDNA/MyHeritage) -> PLINK
02_build_aadr_subset.R          # AADR population subset -> PLINK (memory-safe, pops= filter)
03_merge_and_test.sh            # merge personal genome + AADR subset -> merged_pw_aadr_final_ready
07_extended_models.R            # 1-source models outside Poland_*: Slavic/Baltic/Lusatian +
                                #   Wielbark/ScandinavianIA/Nordic + modern neighbours +
                                #   Ashkenazi/Magyar + Avar/Hun/Sarmatian/Bulgar/Xiongnu
                                #   + [aadr_eu] West/Central-EU batch: England_IA/Saxon/Medieval,
                                #     Ireland/Netherlands/Germany_AltInden/Hungary_Lombards,
                                #     France IA, Italy Roman/Piceni, Spain IA/EBA + modern
                                #     English/French/Italian_North/Spanish/Orcadian/Icelandic
08_two_source.R                 # rotating 2-source scan: well-powered bases × diverse second sources
                                #   (Uralic/Siberian/Steppe/Hunnic/ME/Caucasian +
                                #    Viking/Wielbark/Saxon/Ashkenazi/Magyar); SIGNIFICANT-2src detector
                                #   hardened: MIN_BASE_N>=15 + p>p_base + both w in (0.1,0.9) +
                                #   |z|>3 + BH q_second<0.05 (see hardening note above)
                                #   + [aadr_eu] same West/Central-EU batch as second sources
                                #     (England/Ireland/Netherlands/Germany/France/Italy/Spain +
                                #      modern English/French/Spanish/Orcadian)

examples/site.txt               # pointer to PGP-HMS public genomes (future EU-target examples land here)
README.md                       # beginner's guide (prep 00–03 + steps 07–08)
```

Canonical Polish reference reports (PW/JW) and the original command history live
in `../aadr/`, not here — this fork carries only the EU-target machinery.

## Pipeline invocation

Output directory is always derived from the **input CSV filename**, created
relative to the **working directory at invocation** (not the AADR data dir).
Run everything from `~/Claude` — f2caches and TSVs land in the run dir (e.g. FF_PW/, JW_MH/), AADR stays read-only.

```bash
cd /home/pwas/Claude

# Steps 01-03: prepare personal genome + merge with AADR
./aadr_eu/01_prepare_files.sh /path/to/FTDNA_rawdata.csv        # -> ./FTDNA_PW/
Rscript aadr_eu/02_build_aadr_subset.R \
  /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB ./FTDNA_PW
./aadr_eu/03_merge_and_test.sh ./FTDNA_PW PW                    # arg 2 = SAMPLE_ID!

# Steps 07-08: EU-expanded analysis (both take <merged_prefix> <aadr_prefix> <target>)
# (04-06, 09 are Polish-target probes — run from ../aadr/ if needed)
Rscript aadr_eu/07_extended_models.R  ./FTDNA_PW/merged_pw_aadr_final_ready \
  /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB PW
Rscript aadr_eu/08_two_source.R       ./FTDNA_PW/merged_pw_aadr_final_ready \
  /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB PW
```

**Idempotency:**
- Step 3: skips if `merged_pw_aadr_final_ready.*` already exists.
- Steps 4–8: f2cache in `<basename>_f2cache/` in CWD; re-run is safe (cache reused).
- Step 7 builds its own ext7 fileset; step 8 builds its own 2src fileset — both
  idempotent on the ext/2src BED files.
- **f2cache is target-specific in steps 7–8** (`..._<TARGET>_f2cache/`) — running
  PW then JW from the same CWD does NOT collide.

**For non-PW targets** (e.g. JW):
```bash
SAMPLE_ID=JW ./aadr_eu/01_prepare_files.sh /path/to/JW_rawdata.csv
./aadr_eu/03_merge_and_test.sh ./JW_MH JW      # arg 2 = JW, not PW
```

**Missing `pop_lookup.tsv`** (manual earlier run):
```bash
AADR_PREFIX=/usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
  ./aadr_eu/03_merge_and_test.sh ./FTDNA_PW PW
```

## Hard-won pitfalls — all pre-fixed in the scripts, do not regress

Every fix has a matching comment block in the relevant script.

1. **plink2 needs `--pedmap`**, not bare `--file` (plink1.9 syntax). (Steps 1, 3a.)

2. **plink2 cannot read AADR's native packedancestrymap format.**
   `--pgen/--pvar/--psam` do NOT accept `.geno/.snp/.ind`. (Step 2.)

3. **`admixtools::packedancestrymap_to_plink()` loads the ENTIRE panel (~129 GB
   RAM) — OOM.** Always use `read_packedancestrymap(pops=...)` instead. (Step 2.)

4. **`admixtools` has no exported `write_plink()`** — PED/MAP must be written
   manually from the genotype matrix. (Step 2.)

5. **`read_packedancestrymap()` `$ind` tibble columns are `X1/X2/X3`**, not
   `ID/Sex/Pop`. Access positionally. (Step 2.)

6. **SNPs 100% missing in a small subset → `A1=A2="0"` in `.bim`** — plink
   rejects as "Identical A1 and A2 alleles". Exclude before merge. (Step 3b.)

7. **plink2 `--pmerge`/`--pmerge-list` cannot do non-concatenating merges** —
   "under development." Use plink1.9 (`p-link --bmerge`). (Step 3e.)

8. **plink1.9 `--sort-vars` cannot combine with `--make-bed`** — sort via
   `--make-pgen` first, then convert. (Step 3c.)

9. **Cross-platform allele "conflicts" (real allele vs `.`) mis-flagged as
   "3+ alleles."** `--flip` does NOT fix this. Exclude conflicting SNP IDs from
   both filesets. (Step 3d.)

10. **[RESOLVED] plink1.9 v1.90b2p (2014) segfaulted at ~99% of `--bmerge`.**
    Current build `v1.90b7.11.d` (overlay pwr, `plink-1.90_p20260606`) is clean.
    Binary name is `p-link`, not `plink`. (Step 3e.)

11. **`extract_f2()` `auto_only=TRUE` fails to parse mixed numeric/non-numeric
    chromosome column.** Fix at plink level: `plink2 --chr 1-22` before
    `extract_f2()`. (Step 3f.)

12. **After merge, `.fam` FID defaults to sample ID, not population** —
    `extract_f2(pops=...)` and `qpadm()` match on FID. Step 2 writes
    `aadr_subset_pop_lookup.tsv`; step 3g rewrites FID post-merge. (Step 3g.)

13. **f2cache ląduje w katalogu runu** (`dirname(prefix)`) — obok plików BED,
    nie w CWD. Wszystkie skrypty 04–08 używają `file.path(outdir/dirname(prefix), ...)`.
    Warunek: katalog runu musi być zapisywalny (zawsze prawda dla `MH_PW/`, `JW_MH/`
    itp.; AADR source w `/usr/local/share/aadr/` to tylko wejście, nie prefix merga).
    **Nie uruchamiaj skryptów z katalogu root-owned** — i tak nie ma powodu. (Steps 4–8.)

14. **f2cache target-specyficzny w krokach 7–8** (`..._<TARGET>_f2cache/`) —
    krok 7 bez target-suffix miał kolizję cache PW↔JW przy tym samym basename
    pliku BED. Krok 8 analogicznie. Nie współdziel cache między targetami. (Steps 7–8.)

15. **`options(width=...)` must be set at the top of each script** — otherwise
    R wraps wide summary table columns (z_second, verdict) to next line.
    Currently `options(width=200)` in scripts 06–08. (Steps 6–8.)

## SNP-count expectations

- **Aggregate-population test** (all pops from AADR only): ~500K SNPs after
  `extract_f2(maxmiss=0)`. Normal — same panel throughout.
- **Personal-genome test** (consumer array merged with AADR): **~41–50K SNPs**.
  This is the realistic expectation, not ~16K (that was an early run with fewer
  reference populations). Root cause: FTDNA/MyHeritage chips (GWAS-designed)
  genotype ~520–560K positions; AADR HO panel has 584K positions; physical
  overlap is ~57K, then `maxmiss=0` across all ancient populations (with their
  own missingness) brings it to ~41–50K. **This is expected, not a bug.**
  - sub-30K: check self-genotyping-rate from step 1 (F_MISS should be 0).
  - MH+FTDNA union merge: gives ~49K — same as individual arrays, no improvement
    (overlap with AADR HO is the bottleneck, not array size).
- **WGS 30× subsetted to AADR panel positions**: ~200–400K SNPs expected (5–10×
  improvement); new bottleneck = ancient missingness, not the personal genome.
  Use `bcftools mpileup -R <1240K.sites>` — do NOT use all WGS positions
  (breaks f-statistic ascertainment calibration). Prefer 1240K panel over HO
  for WGS (more ancient coverage); needs hg19 alignment or liftover.

## Known AADR label gaps

- No `Germany_IA` matching G25/Davidski's Saxony-Anhalt IA cluster. Closest
  proxy: `Germany_Esperstedt_CordedWare` (n=13) — tests the CWC substrate
  hypothesis, not IA Germanic ancestry specifically. `Germany_Singen_IA` is
  Baden-Württemberg/Hallstatt (n=1, unusable).
- Polish labels coarser than G25 clusters: `Poland_EarlyMedieval_Slav` (n=46),
  `Poland_EarlySlav` (n=25), `Poland_IA_Wielbark` (n=59), `Poland_IA_Przeworsk`
  (n=2 — `PCA0011.SG`/`PCA0012.SG`, both female, Gąski/Kuyavia 100–300 CE, high-cov
  shotgun 160K/100K HO SNPs; tested only in step 9, skipped by steps 6–7's n≥3 gate),
  `Poland_Medieval-1/-2` (n=1 each). No 1:1 equivalents for Davidski's
  per-site clusters (Santok, Markowice, Płońsk, Gródek).
- `Ukraine_IA_Lusatian` (n=3, Rovantsi, Wołyń, 1000–700 BC) — best single-source
  fit for both PW and JW; n=3 means Y-haplogroup coverage is poor (1 usable male:
  „R1a1'5", too coarse to confirm CTS1211).
- Y-haplogroup data (`.anno` col 35, YFull 12.03 auto-called): YP343 (on PW's
  paternal path) found in `Poland_EarlyMedieval_Slav` n=1 only. Each ancient male
  typically gets a unique terminal label — n=1 per exact branch is normal.

## Tool versions and binary names

- R: `admixtools` (uqrmaie1/admixtools via `remotes::install_github`)
- `plink2`: `sci-biology/plink2-bin`, binary `plink2`
- `plink1.9`: `sci-biology/plink` (overlay pwr, `plink-1.90_p20260606`), binary
  **`p-link`** (not `plink`). Build `v1.90b7.11.d`, 2026-06-06. Segfault resolved.
- AADR: v66.p1 HO (`v66.p1_HO.aadr.patch.PUB.{geno,snp,ind}` +
  `v66.p1_HO.aadr.PUB.anno`), `/usr/local/share/aadr/` (root-owned, read-only).
  HO chosen over 1240K: includes modern outgroups (Mbuti, Yoruba, Han) needed
  for qpAdm right-pops. Downloaded from Harvard Dataverse `doi:10.7910/DVN/FFIDCW`
  — no `download_aadr()` R function, use Dataverse file API directly.

## What not to do

- Don't `packedancestrymap_to_plink()` on the full panel — OOM (~129 GB).
- Don't use plink2 `--pmerge` for AADR-vs-personal-genome merge — not supported.
- Don't retry `--flip` for "3+ alleles" from one-sided missingness — exclude instead.
- Don't interpret ~41–50K SNPs as a bug — it's the correct expectation for consumer arrays.
- Don't run steps 4–8 expecting f2cache in CWD — caches now go to `dirname(prefix)` (the run dir). Old caches in `~/Claude/` are stale and can be deleted.
- Don't re-litigate Germany_IA — closed, confirmed on two independent targets.
- Don't re-litigate Avar/Pannonian ancestry — closed by steps 7+8, same mechanism
  as Germany_IA (structural collinearity, not data insufficiency).
- Don't report G25/nMonte **weighted or mixed** model components (Hungary_Mig_Himod,
  Latvia_Neo_Kivutkalns, Latvia_Viking) as real ancestry — these are bootstrap-stable
  attractors from noisy tail PCs, falsified by step 8 (0 significant 2-source models).
  Use unweighted historical nMonte (L2) for reporting.
- Don't share f2cache between PW and JW runs in steps 7–8 — target-specific dirs required.
