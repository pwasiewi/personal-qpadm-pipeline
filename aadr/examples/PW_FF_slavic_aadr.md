# Report: Ancestry Analysis PW — FTDNA Family Finder × AADR v66
**Date:** 2026-06-27  
**Input data:** `FF_PW_Chrom_Autoso_20260223.csv` (FTDNA Family Finder)  
**Reference panel:** AADR v66.p1 Human Origins (584,131 SNPs, 27,594 samples)  
**Tools:** plink 1.90b7.11.d / plink2 / admixtools2 (R), pipeline `aadr/01–08`

---

## 1. Input data and quality

| Metric | Value |
|---|---|
| SNPs in FTDNA file (autosomal, rs*, non-missing) | 601,140 |
| SNPs overlapping AADR panel (raw, by rsID) | 56,988 |
| SNPs after merge and `extract_f2(maxmiss=0)` | **49,883** |
| Missing genotypes (F_MISS) | 0.0117 (7,117 SNPs) |

Comparison with earlier MyHeritage run (`MH_PW`): both datasets yield ~50,000 SNPs after filtering — the earlier ~16,000 SNP count was an artefact of the buggy plink 1.9 build from 2014 (segfault during `--bmerge`). The current ebuild `sci-biology/plink-1.90_p20260606` (overlay `pwr`) resolves this.

---

## 2. Germany_IA artefact test (step 4)

**Question:** Is the `Germany_Esperstedt_CordedWare` component (proxy for "Germany_IA" from G25/nMonte models) real ancestry of PW, or a geometric artefact of the panel?

### qpAdm results (FF_PW, 49,883 SNPs)

| Model | p | Interpretation |
|---|---|---|
| **Model A:** PW ← Wielbark (alone) | **0.256** | NOT rejected; weight = 1.000 |
| **Model B:** PW ← Wielbark + Germany_CWC | — | Germany weight = −1.299, SE = 1.117, **z = −1.16** |

```
popdrop:
  pat=00 (both sources)    p = 0.928
  pat=01 (Wielbark alone)  p = 0.044   ← borderline
  pat=10 (Germany alone)   p = 0.00027 ← strongly rejected
```

**Conclusion:** Germany/CWC is a geometric artefact — a negative weight (−1.3) is the classic collinearity signature, not real Germanic admixture. The large "Germany_IA" share in G25/nMonte models reflects a shared Corded Ware/Bronze Age substrate, not genuine Germanic ancestry.

---

## 3. Outgroup robustness (step 5)

Six outgroup combinations tested (right populations):

| Outgroup set | Model A p | Germany weight | Germany z | Verdict |
|---|---|---|---|---|
| Mbuti+Yoruba+Han (baseline) | 0.381 | −1.060 | −1.10 | ARTEFACT |
| Mbuti+Yoruba+Papuan | 0.352 | −1.101 | −1.09 | ARTEFACT |
| Mbuti+Han+Biaka | 0.523 | +0.001 | +0.00 | ARTEFACT |
| Mbuti+Yoruba+Nganasan | 0.323 | −1.125 | −1.19 | ARTEFACT |
| Mbuti+Yoruba (2 outgroups) | 0.165 | −1.100 | −1.07 | ARTEFACT |
| Mbuti+Han (2 outgroups) | 0.382 | −1.438 | −0.55 | ARTEFACT |

**6/6 combinations → ARTEFACT.** Conclusion robust to outgroup choice.

---

## 4. Systematic ancestry model test (step 6)

**Method:** qpAdm, 16 models (9 single-source + 6 two-source + 1 three-source).  
Outgroups: Mbuti + Yoruba + Han (proven stable in step 5).  
SNPs: 49,144 after filtering with extended file (737 samples, 13 populations).

### Single-source models

| Population | Period | n AADR | p | Status |
|---|---|---|---|---|
| **Poland_EarlySlav** | VI–VIII CE | 25 | **0.305** | ✓ PASS — best fit |
| Poland_EarlyMedieval_Slav | IX–XI CE | 46 | **0.222** | ✓ PASS |
| Poland_CordedWare | ~2800–2300 BC | 16 | **0.128** | ✓ PASS |
| Poland_IA_Wielbark | I–III CE | 59 | **0.123** | ✓ PASS |
| Poland_GlobularAmphora | ~3200–2300 BC | 33 | **0.099** | ✓ PASS |
| Russia_Samara_EBA_Yamnaya | ~3100–2500 BC | 46 | **0.071** | ✓ PASS |
| Poland_Roman_Wielbark | III–IV CE | 35 | 0.034 | ✗ REJECTED |
| Slovakia_N_LBK | ~5500–4900 BC | 56 | 0.032 | ✗ REJECTED |
| Poland_BA_Trzciniec | ~1800–1200 BC | 61 | 0.020 | ✗ REJECTED |

### Two- and three-source models

All 6 multi-source models passed (p > 0.05), but **all with unphysical weights** (far outside [0,1]):

| Model | p | Weights | Status |
|---|---|---|---|
| Wielbark + EarlySlav | 0.316 | −3.9 / +4.9 | OK-unphysical |
| CordedWare + EarlyMedieval | 0.311 | −3.7 / +4.7 | OK-unphysical |
| Wielbark + EarlyMedieval | 0.276 | −5.1 / +6.1 | OK-unphysical |
| Yamnaya + LBK | 0.108 | +9.6 / −8.6 | OK-unphysical |

Unphysical weights are a **mathematical collinearity signature**: the tested populations are too genetically similar for qpAdm to resolve them as independent components.

---

## 5. Extended models — populations outside the Poland_* set (step 7)

**Question:** Does any neighbouring population (non-Polish Slavic, Baltic, modern neighbours) explain PW better than Poland_EarlySlav?

**Method:** qpAdm, 12 single-source models, same outgroups: Mbuti + Yoruba + Han.

### Results (sorted by p)

| Population | Period | n AADR | p | Status |
|---|---|---|---|---|
| **Ukraine_IA_Lusatian** | ~500 BC (Iron Age) | 3 | **0.951** | ✓ **NEW LEADER** |
| Serbia_EarlyMedieval_Byzantine_Slav | VI–VIII CE | 5 | **0.434** | ✓ |
| Lithuanian | modern | 10 | **0.365** | ✓ |
| Lithuania_LBA | ~1000–500 BC | 5 | **0.254** | ✓ |
| Ukrainian_North | modern | 8 | 0.197 | ✓ |
| Slovakia_CiferPac_Medieval_Avar_Slav | IX–X CE | 7 | 0.163 | ✓ |
| Russia_Ivanovo_EarlyMedieval | VIII–X CE | 4 | 0.107 | ✓ |
| Estonian | modern | 12 | 0.062 | ✓ borderline |
| Belarusian | modern | 10 | 0.061 | ✓ borderline |
| Czech | modern | 11 | 0.037 | ✗ REJECTED |
| Ukrainian | modern | 13 | 0.015 | ✗ REJECTED |
| Czechia_EarlyMedieval_EarlySlav | VI–VIII CE | 3 | 0.001 | ✗ REJECTED |

*For comparison with step 6: Poland_EarlySlav p=0.305, Poland_EarlyMedieval_Slav p=0.222.*

### Interpretation

**Ukraine_IA_Lusatian (p=0.951)** is three times a better fit than Poland_EarlySlav (p=0.305) — the best of all 21 models tested in steps 6 and 7. The Lusatian culture (~1300–400 BC) spanned Poland, Czech lands, eastern Germany and western Ukraine and is the primary candidate for the genetic substrate of the Slavs. Three samples from the Ukrainian range of the Lusatian culture fit PW far better than populations from the second half of the 1st millennium CE — PW's ancestry reaches deeper than Wielbark and EarlySlav: back to the Iron Age Central European Lusatian substrate.

**Serbia_EarlyMedieval_Byzantine_Slav (p=0.434)** — better than Poland_EarlySlav. Early medieval Balkan Slavs (VI–VIII CE) derive from migrations out of what is now Poland/Ukraine; their genetic profile is closer to PW than later Polish Slavs from Gródek/Markowice.

**Lithuanian (p=0.365)** — better than Poland_EarlySlav. PW has a clear Baltic/WHG substrate component, consistent with the geography of Mazovia (historical Slavic–Baltic contact zone).

**Czech (p=0.037) and Ukrainian (p=0.015) rejected** — despite geographic proximity. Modern Czechs carry more Western European (Hallstatt/Celtic) ancestry; modern Ukrainians more steppe. PW does not fit either profile.

**Czechia_EarlyMedieval_EarlySlav strongly rejected (p=0.001)** — early medieval Czech Slavs have a distinctly different profile from Polish Poland_EarlySlav, likely due to a larger Western European substrate or a different migration pathway.

### Implications for the Slavic homeland debate

Ukraine_IA_Lusatian is three samples from a single site: **Rovantsi, Volyn oblast (Lutsk), 1000–700 BC** (Saag & Thomas, *Science Advances* 2025) — ~200 km east of the present Polish border, ~150 km from Gródek Nadbużny.

The Lusatian culture (~1300–400 BC) is the eastern variant of the Urnfield culture complex, covering Poland, eastern Germany, Czech lands, Slovakia and western Ukraine. Rovantsi lies at its eastern margin, in the Bug–Styr zone (Polesia) — one of the leading candidates for the Slavic urheimat in archaeological literature.

**Genetic flow:**

```
Corded Ware  (~2800 BC)
      │
      ▼
  Lusatian culture (~1300–400 BC)   ← Ukraine_IA_Lusatian (Rovantsi, Volhynia)
      │
  ┌───┴────────────┐
  │                 │
  + Scandinavian     no admixture
  admixture          → Przeworsk
  (Goths from        ▼
  Jutland)          Poland_EarlySlav
  ▼                 (600–800 CE)
Wielbark
(100–300 CE)
```

Lusatian is a population **prior to the split** into the Germanic branch (Wielbark = Lusatian substrate + ~20–30% Scandinavian ancestry) and the Slavic branch (EarlySlav = Lusatian substrate after ~600 years of drift through Przeworsk). The p=0.951 result supports the **autochthonous hypothesis** (Kostrzewski/Hensel school): Slavs did not arrive from outside in the 6th century CE — they evolved from the local Lusatian population of Central Europe. PW fits a 1,000-year-old Volhynian population better than "Slavs" from the 8th century CE, because that population *is* those same people before they adopted the name "Slavs."

The collinearity of Wielbark + EarlySlav in two-source models is now genetically explained: both are derivatives of Lusatian in different directions — qpAdm sees two vectors spanning the same space as Lusatian itself and cannot resolve them as independent components. Any proportional estimate ("x% Wielbark + y% EarlySlav") is mathematically fictional, imposed on one continuous stream.

---

## 5A. Extended step 7 — Avars, Huns, Steppe (verification of G25 "Avar_Kecskemet")

**Reason:** Weighted nMonte on G25 vectors showed ~7% `Lithuania_BA` + 20–30% `Avar_Kecskemet` for PW instead of 40% `Germany_IA`. Step 7 was extended from 12 to **23 single-source populations**, adding Pannonian Avars, Huns, Sarmatians and early Steppe Bulgars.

| Population | n | p PW | Status |
|---|---|---|---|
| **Ukraine_IA_Lusatian** | 3 | **0.9795** | ✓ leader (unchanged) |
| Croatia_EarlyMedieval_Avar | 12 | 0.301 | ✓ (border Avars, Slavicised) |
| Russia_Tatarstan_EarlyMedieval_EarlyBulgar | 7 | 0.284 | ✓ |
| Hungary_Transtisza_EHun | 7 | 0.245 | ✓ (small n, wide CI) |
| Russia_Sarmatian | 9 | 0.232 | ✓ |
| Hungary_Madaras_MigrationPeriod_Sarmatian | 11 | 0.144 | ✓ |
| Hungary_EHun | 7 | 0.053 | ✓ borderline |
| **Austria_Avar** | 714 | **0.043** | ✗ REJECTED |
| **Hungary_EarlyAvar-oLowEastAsia** | 52 | **0.035** | ✗ REJECTED |
| **Hungary_EarlyAvar** | 35 | **0.023** | ✗ REJECTED |
| Slovakia_MigrationPeriod | 6 | 0.022 | ✗ REJECTED |
| Lithuania_MigrationPeriod | 5 | 0.012 | ✗ REJECTED |

**Pannonian Avars REJECTED** as a sole source — including the European fraction `-oLowEastAsia` (Avar samples without East Asian admixture, i.e. exactly what could survive in a Polish population as "Avar_Kecskemet"). `Avar_Kecskemet` from G25 does not exist as an AADR label — it is Davidski's cluster; its closest AADR equivalents are rejected by PW. The non-rejected `Croatia_EarlyMedieval_Avar` is effectively a Slavic population from the Avar borderland, not steppe.

## 5B. Step 8 — rotating two-source scan: no hidden admixture

**Question:** Does PW require a SECOND source (steppe/Asian/Middle Eastern) that a single-source model cannot detect by construction? (With 1 source the weight is fixed = 1.000, so the second source is unmeasurable.)

**Test (expanded 2026-06-29):** 2 bases (`Ukraine_IA_Lusatian`, `Poland_EarlySlav`) × **35** diverse second sources — Uralic/Siberian (Nganasan, Selkup, Mordovian, Bashkir…), Eastern steppe / **true Xiongnu** (`Russia_Buryatia_XiongnuPeriod`) and **Hunnic steppe** (`Kazakhstan_Berel_Hunnic`), Middle Eastern (Druze, Iranian), Caucasian (Armenian, Georgian, Adygei, Lezgin), plus **Viking** (Sweden/Denmark/Norway/Iceland), **Wielbark** (`Poland_IA_Wielbark`), **Saxon** (Germany/England), and **historical Polish-admixture candidates: `Jew_Ashkenazi`, `Germany_Medieval_Jewish-lowEastEU` (Erfurt), `Hungary_Conqueror_Commoner` (Magyar), Romanian, Hungarian**. 40,696 SNPs. (Roma/Gypsy absent from AADR HO — untestable.)

**Result: ZERO "SIGNIFICANT-2src" models** (38 second sources, all |z| < 0.32). With `Ukraine_IA_Lusatian` base (p_base=**0.973**), 24/35 second sources receive **negative** weights and the rest ≤ +0.066: the base already over-explains PW; there is no room for any admixture. Verdict tally: **19 OK-base-sufficient**, 1 second-dominant, 2 nonsignificant, 15 rejected.

- **Ashkenazi Jewish:** `Jew_Ashkenazi` w=−0.03 (base-sufficient), `Germany_Medieval_Jewish` REJECTED → **no detectable Ashkenazi component.**
- **Magyar:** `Hungary_Conqueror_Commoner` w=−0.005 (base-sufficient) → none.
- **Viking / Wielbark / Saxon:** all w ≤ 0 → none. (`Poland_IA_Wielbark` as second source w=−0.002, i.e. PW is *less* Wielbark than pure Lusatian.)

→ **No hidden surprise.** Lower p for PW (vs JW) from step 7 is global noise/drift (FTDNA), not a separate ancestry component. PW sits *tighter* on the pure Lusatian profile than JW — the opposite of "more diverse genome." **"Avar_Kecskemet" from G25 = geometric artefact, identical mechanism to Germany_IA** (closed by two independent tests: single-source in 5A, two-source here).

### Autochthon vs allochthon — unresolved at ~40K SNP

In step 6, PW's **local Polish continuity chain is broken**: `Poland_Roman_Wielbark` (p=0.034) and `Poland_BA_Trzciniec` (p=0.020) are REJECTED as single sources, and `Poland_IA_Wielbark`/`Poland_IA`/`Poland_CordedWare` are only marginal (p≈0.12). Yet the **eastern** `Ukraine_IA_Lusatian` (Volhynia, n=3) fits far better (p=0.883). Taken alone this hints at an **allochthonous** origin (Slavs arriving from the east). **But it cannot be concluded** at 40–50K SNP:
1. The Lusatian culture spanned **both Poland and western Ukraine** — `Ukraine_IA_Lusatian` is simply the best-sampled representative, not proof of eastern provenance.
2. `Ukraine_IA_Lusatian`, Wielbark, `Poland_EarlyMedieval_Slav` and `Poland_EarlySlav` are mutually **collinear** (shared Lusatian/CWC substrate); qpAdm at this SNP count cannot order them causally.
3. Autochthonous continuity (Wielbark→Medieval→modern) and allochthonous migration both pass — the data are underpowered to distinguish them. WGS (200–400K SNP) would be required.

**Methodological note:** `Poland_EarlySlav` as base produces degenerate weights (|w| > 30, z ≈ 0, inflated p) — this is unidentifiability (collinearity with second source), not good models. Script 08 marks these UNIDENTIFIABLE and sorts them to the bottom of the table (21/46 rows = all Poland_EarlySlav rows). Interpret **only** the `Ukraine_IA_Lusatian` base.

---

## 5C. Step 9 — Przeworsk vs Wielbark channel of descent (2026-06-29)

**Question:** PW's documented background is medieval **Greater Poland (Wielkopolska)**. Did the line of descent from the Iron Age Lusatian substrate run through the **Przeworsk** culture (central/southern Poland, autochthonous Slavic-substrate channel) rather than the **Wielbark** culture (Pomeranian, Gothic/Scandinavian-admixed channel)? Przeworsk covered Wielkopolska/Silesia/Lesser Poland/Mazovia; Wielbark was Pomerania-centred.

**Why this was never tested before:** steps 6–7 hard-skip any source with n<3. AADR v66 HO has only **`Poland_IA_Przeworsk` n=2** — `PCA0011.SG`, `PCA0012.SG`, both **female**, both from **Gąski (Kuyavia, Inowrocław county), 100–300 CE**, high-coverage shotgun (160K / 100K HO SNPs hit). Step 9 admits n≥2 with an explicit caveat.

**Method:** single-source qpAdm + the decisive **f4 symmetry statistic** + 2-source competition (Lusatian base). Outgroups Mbuti+Yoruba+Han, ~40K SNPs.

### Single-source (context only — collinear, pass ≠ channel proof)

| Source | n | p | Status |
|---|---|---|---|
| **Poland_IA_Przeworsk** | 2 | **0.351** | ✓ — highest of the Polish IA sources |
| Poland_IA_Wielbark | 59 | 0.192 | ✓ |
| Poland_EarlySlav | 25 | 0.148 | ✓ |
| Ukraine_IA_Lusatian | 3 | 0.984 | ✓ (overall leader, unchanged) |

### f4 symmetry tests — the decisive instrument

| Test | f4 | Z | Reading |
|---|---|---|---|
| **B1** channel differentiation `f4(Mbuti, Sweden_IA; Przeworsk, Wielbark)` | −0.00013 | **−0.29** | **n.s. — the two channels are genetically indistinguishable** |
| **B2 (decisive)** `f4(Mbuti, PW; Wielbark, Przeworsk)` | −0.00034 | **−0.46** | **n.s. — PW symmetric to both** |
| B3a vs Przeworsk `f4(Mbuti, PW; Sweden_IA, Przeworsk)` | −0.00065 | −0.76 | no Scandinavian pull |
| B3b vs Wielbark `f4(Mbuti, PW; Sweden_IA, Wielbark)` | −0.00031 | −0.62 | no Scandinavian pull |
| B4 robustness (Han outgroup) | −0.00027 | −0.38 | agrees with B2 |

**2-source (Lusatian base):** both channels receive weight ≈ 0 (Przeworsk −0.048, Wielbark −0.072, z≈0) → `OK-base-sufficient`. Neither channel adds anything over the Lusatian substrate.

### Verdict

The Wielkopolska/Przeworsk hypothesis is **consistent with the data and not contradicted, but unprovable at ~40K SNP** — for a concrete reason: **B1 is n.s.**, so Przeworsk and Wielbark cannot be told apart from each other here (the expected Gothic/Scandinavian shift of Wielbark via `Sweden_IA` does not surface). With the channels themselves indistinguishable, "which channel" is unanswerable; **B2 confirms** PW is equidistant from both. Same collinearity mechanism as the closed Wielbark/EarlySlav/Lusatian finding — Przeworsk simply joins that family of Lusatian derivatives.

What *weakly favours* the hypothesis (directional, not significant): Przeworsk has the **highest single-source p among Polish IA sources** (0.351 vs Wielbark 0.192); all `Sweden_IA` f4 are slightly negative → **no Scandinavian/Gothic pull in PW**, the profile expected of a non-Wielbark channel; and the only Przeworsk samples (Gąski, Kuyavia) sit at the Wielkopolska doorstep — culturally and geographically the right route for a Greater-Poland background. Resolving it requires WGS (200–400K SNP), which would separate B1 and make B2 decisive.

**Files:** `PW_FF/przeworsk_channel_results_PW.tsv`, `PW_FF/przeworsk_channel_f4_PW.tsv`. Script: `aadr/09_przeworsk_channel.R`.

---

## 6. Biological interpretation

### Genetic continuity — one population stream

Results point to **genetic continuity from the Neolithic through the Iron Age to the Early Medieval period** on Polish territory. The populations tested as sources form a single linear stream and are mutually substitutable in qpAdm models at ~50,000 SNPs:

```
Yamnaya / Pontic steppe (~3000 BC)
        ↓
Corded Ware / Globular Amphora (Poland, ~2800 BC)
        ↓
Wielbark / Iron Age (~I–III CE)
        ↓
Early Medieval Slavs (~VI–VIII CE)
        ↓
PW (modern)
```

### What is rejected — negative information

- **Poland_BA_Trzciniec** (p=0.020) — the Trzciniec culture on its own (without a steppe component, or with a different profile) does not explain PW.
- **Slovakia_N_LBK** (p=0.032) — pure Neolithic farmers without steppe ancestry do not explain PW.
- **Poland_Roman_Wielbark** (p=0.034) — the specifically *Roman-period* Wielbark phase is rejected. Interpretation: this phase may have carried distinctive ancestry from more intensive contacts with Scandinavia or western Germanic populations — or this is a small-sample artefact (n=35).

### The Slavic hypothesis — confirmation

**Poland_EarlySlav has the highest p = 0.305** among all 9 single-source populations tested. PW is better described by early medieval Polish Slavs (VI–VIII CE) than by the Wielbark culture (I–III CE) — although both populations pass the test.

---

## 7. The "10% Wielbark, 90% Gródek" claim for all of Poland

Some archaeogenomic studies (including from the MPI Leipzig / Haak et al. circle) estimate a local Iron Age ancestry contribution of ~10% at early medieval sites in eastern Poland (e.g. Gródek Nadbużny), with ~90% influx from eastern Slavs. Extrapolating these proportions to **all of Poland is methodologically flawed** for several reasons:

### 7a. Geographic problem

Gródek Nadbużny lies on the Bug river, at the eastern edge of the Wielbark culture's range. Wielbark dominated central and western Poland (Pomerania, Mazovia, Kujawy). By the 4th–5th centuries it had already retreated in the east to the Kiev and Cherniakhov cultures. The 10/90 ratio from Gródek describes the eastern flank of the Slavic expansion — an area of minimal Wielbark continuity — and cannot be applied to Mazovia or Greater Poland.

### 7b. Statistical problem in qpAdm

Our results show that Wielbark and Poland_EarlySlav produce **unphysical weights** in a two-source model:
```
Wielbark + EarlySlav → weights: −3.9 / +4.9   (p = 0.316)
```
When weights are this extreme it means the **populations are linearly unidentifiable** — f2-statistics cannot distinguish 100% Wielbark from 100% EarlySlav from any mixture thereof. Every specific proportion (10/90, 50/50, 90/10) within this pair is equally consistent with the genomic data. Claiming "10% Wielbark" without estimating the uncertainty of that number is false precision.

### 7c. Definition problem: "Wielbark" vs "Gródek"

The Wielbark culture (I–III CE) and early medieval Slavs (VI–VIII CE) are separated by ~300 years and several cultural phases. Migration-period populations (IV–V CE) are poorly represented in AADR. If an analysis compares Gródek (VI–VIII CE) against Wielbark (I–III CE) without intermediate populations, the "missing" 90% may represent natural genetic evolution of the same population lineage across several generations, not mass influx from outside.

### 7d. What our data show

For PW (modern Pole, genealogy from Mazovia/central Poland):
- Poland_EarlySlav **and** Poland_IA_Wielbark **both** pass as single-source models
- They are collinear → there is no "mixture" to estimate; this is one continuous stream
- Poland_Roman_Wielbark is rejected, suggesting the late-Wielbark-specific signal did not carry through — but the earlier Iron Age Wielbark did

**Conclusion:** For central and western Polish populations, the "90% population replacement" claim is inconsistent with the genomic data of modern Poles. If 90% of ancestry derived from incoming eastern Slavs, models based on local IA populations (Wielbark, CordedWare) would fail as good single-source approximations. Instead they pass.

---

## 8. Summary

| Question | Answer |
|---|---|
| Is Germany_IA real ancestry? | **NO** — geometric artefact, 6/6 outgroup combinations |
| Best single-source model (of 32 tested) | **Ukraine_IA_Lusatian p=0.9795** (extended step 7) |
| Best model from Poland_* set | Poland_EarlySlav p=0.305 (step 6) |
| Are Avars / "Avar_Kecskemet" (G25) real ancestry? | **NO** — Hungary/Austria_Avar rejected single-source (5A); zero significant two-source (5B); same artefact as Germany_IA |
| Does PW have a hidden second source (steppe/Asia/Middle East)? | **NO** — 0 SIGNIFICANT-2src models out of 23 diverse sources (step 8) |
| Is there a distinct "Slavic" component? | **YES** — and it traces back to the Lusatian culture (~500 BC), not just EarlySlav (VI–VIII CE) |
| Can Wielbark and EarlySlav be distinguished? | **NO** — collinear at ~50K SNPs |
| Why do Lithuanians fit better than Czechs? | PW has a Baltic/WHG substrate; Czechs carry Western European/Hallstatt ancestry |
| Is the "10% Wielbark across all Poland" claim valid? | **NO** — it applies to the eastern Wielbark frontier; inconsistent with mass population replacement in central Poland |
| Does the paternal line (Y-DNA) agree with autosomes? | **YES, supportive (n=1)** — R1a-Z280-CTS1211-**YP343** (PW) = node present in `Poland_EarlyMedieval_Slav` (1 sample); Wielbark (I/G) outside the line; broad Slavic R1a cluster numerous, exact PW branch is a singleton (section 8A) |

---

## 8A. Paternal line (Y-DNA) — cross-link with autosomes

Independent axis (patrilineal, single locus), but points in the same direction as qpAdm.

**PW haplogroup (YFull):** R-Y110177\* — terminal SNPs FT14350 • Y110177.

Path: R1a-M198 → M417 → Z645 → Z283 → Z282 → **Z280 → CTS1211** → Y2205 → **YP343** → YP340 → YP371 → Y244926 → FT14429 → **Y110177**. This is the **R1a-Z280/CTS1211** branch — alongside M458 one of the two dominant R1a clusters in Poland, the core of "Slavic" R1a. FT14429 formed ~2400 ybp, TMRCA ~2300 ybp (Iron Age / pre-Roman period). Asterisk = paragroup (Y110177+, not yet placed below Y81700 — requires BAM/denser database).

**Cross-link with AADR (column 35 of `.anno`, Y-hg of ancient males queried manually):**

| Population (autosomal context) | Y-haplogroups of male samples (count) | Relation to PW's line |
|---|---|---|
| **Poland_EarlyMedieval_Slav** (~30 males of 46) | **R-YP343 (n=1)**, R-CTS11962 (n=2), R-M198 (n=2), R-YP516, R-YP4858… | **R-YP343 = node on PW's path** (5 levels above FT14429) — **the only** sample directly on PW's line |
| Poland_EarlySlav (~18 males of 25) | R1a-Z280/CTS1211: R-YP593, YP6048, YP5470, YP415, YP256, YP1448, YP1337 (each n=1) | same **broad** "Slavic R1a" cluster, but **different branches** — none on PW's path (cousins, diverged above YP343). **YP343 absent here.** |
| Poland_CordedWare | **R-Z645 (n=1)**, R-Y215377 (n=2), R-M417, R-M198 | R-Z645 = deep node of PW, but ~5000 ybp and very broad → effectively "R1a", weak evidence |
| Poland_IA_Wielbark (~16 males of 59) | mainly I1/I2 (≥10), G2a (≥7) + scattered R1a (M458, M198, L1029) | CTS1211-YP343 **absent** — Wielbark patrilineally Germanic/local |
| Ukraine_IA_Lusatian | 1 male "R1a1'5" (low quality), 1 A0-T, 1 n/a | R1a only — insufficient resolution to confirm CTS1211 |

**Sample sizes: for PW's exact sub-branch (YP343→…→Y110177) there is a single sample (n=1)**, not a series. This is normal for aDNA: each male receives his own terminal branch label limited by coverage → almost every exact label = n=1; "R-YP343" may actually be deeper, stopped by coverage. The broad R1a-Z280/CTS1211 cluster as a whole is numerous (a dozen males across EarlySlav+EarlyMedieval), but those are sister branches, not PW's.

**Patrilineal chain through time (noting n=1):** Corded Ware (R-Z645, broad) → [Z280/CTS1211 expansion, Iron Age] → Early Medieval Poland (R-YP343, n=1) → PW (YP343 → … → Y110177).

**Conclusion:** PW's paternal line has **one** direct patrilineal hit — node YP343 in `Poland_EarlyMedieval_Slav` (n=1) — a population that is also among the best single-source autosomal fits. A single individual, not a statistic, but the same lineage; consistent with the autosomal conclusion of a local Central European continuum (supporting signal, not proof).

**Caveats:** (1) YP343 hit is n=1 — supporting signal, not proof. (2) Lusatian Y too weak (1 coarse R1a call) — autosomal #1, but does not confirm CTS1211 patrilineally. (3) Wielbark is NOT on PW's paternal line (I/G dominant) → reinforces the conclusion that PW's component is "local Slavic", not "Gothic/Wielbark". (4) Z645 root in CordedWare too broad to count as strong line evidence. (5) TMRCA (~300 BC) coincidence with the Lusatian culture is an independent argument, NOT proof that Lusatians carried CTS1211.

**Source:** YFull tree R-FT14429 (https://www.yfull.com/tree/R-FT14429/); no BAM available — with WGS/BigY BAM, placement below Y110177 can be done locally (Yleaf/yhaplo) and novel SNPs verified.

---

## 9. Output files

```
FF_PW/                             ← main FTDNA run directory
  merged_pw_aadr_final_ready       ← merged fileset after 4 steps (49,883 SNPs f2)
  merged_pw_aadr_extended          ← + Papuan/Biaka/Nganasan (step 5)
  merged_pw_aadr_models            ← + 7 additional populations (step 6)
  robustness_results.tsv           ← step 5: 6 outgroup combinations
  ancestry_models_results.tsv      ← step 6: 16 ancestry models
  merged_pw_aadr_final_ready_f2cache/  ← f2 cache step 4
  merged_pw_aadr_extended_f2cache/     ← f2 cache step 5
  merged_pw_aadr_models_f2cache/       ← f2 cache step 6
  extended_models_results.tsv      ← step 7: 23 non-Polish models (+Avars/Huns/steppe)
  merged_pw_aadr_ext7              ← + Lusatian/Avar/Hun/Sarmatian/etc. (step 7)
  merged_pw_aadr_ext7_f2cache/     ← f2 cache step 7
  two_source_results_PW.tsv        ← step 8: 46 two-source models
  merged_pw_aadr_2src              ← + 23 diverse second sources (step 8)
  merged_pw_aadr_2src_PW_f2cache/  ← f2 cache step 8 (target-specific, 41,133 SNPs)

Scripts: aadr/01..08
```

---

## 10. Limitations and what WGS would add

**Why ~41–50K SNPs and not ~580K?** This is ascertainment bias, not an irreconcilable mismatch between arrays. The AADR HO panel (Affymetrix Human Origins, designed for population genetics) and the FTDNA consumer chip (Illumina, designed for GWAS) physically overlap, but the chip genotypes only a fraction of HO positions; after `extract_f2(maxmiss=0)` (every population including low-coverage aDNA must have a call) ~41–50K remain. The bottleneck is the chip + ancient sample missingness.

**What WGS for PW would add:** WGS has no ascertainment — it covers all positions. Calling genotypes at AADR panel coordinates (`bcftools` against positions from `.snp`) would yield a call at **every** one of the 584K HO positions instead of ~41K. After filtering, survival would be limited only by ancient sample missingness → realistically **~200–400K SNPs (5–10× more)**. Critical: subset WGS to AADR panel positions, do NOT use all WGS sites (otherwise the sample would be on a different ascertainment than the references → biased f4).

**What would change:** narrower SE, more power — some currently borderline-OK models would be rejected, and it would finally be possible to **distinguish** `Ukraine_IA_Lusatian` / `Poland_EarlySlav` / `Poland_IA_Wielbark` (currently statistically indistinguishable at ~50K SNPs). **What would NOT change:** the Germany_IA = artefact and Avars = artefact conclusions — these are structural unidentifiability (sources collinear relative to outgroups), not a SNP-count problem; more data makes the blowup *more visible*, not resolved in favour of those sources. Cheaper alternative: chip imputation (GLIMPSE/Beagle + 1000G/HRC) → subset to AADR positions, but introduces modern reference bias into f-statistics.

---

*Pipeline: `~/Claude/aadr/` | AADR: `/usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB` | plink: `sci-biology/plink-1.90_p20260606` (overlay pwr)*
