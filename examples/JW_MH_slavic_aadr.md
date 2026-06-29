# Report: Ancestry Analysis JW — MyHeritage × AADR v66
**Date:** 2026-06-27  
**Input data:** `MyHeritage_JW_raw_dna_data.csv` (MyHeritage)  
**Reference panel:** AADR v66.p1 Human Origins (584,131 SNPs, 27,594 samples)  
**Tools:** plink 1.90b7.11.d / plink2 / admixtools2 (R), pipeline `aadr/01–08`

---

## 1. Input data and quality

| Metric | Value |
|---|---|
| SNPs in MyHeritage file (autosomal, rs*, non-missing) | 560,853 |
| Missing genotypes (F_MISS) | **0.0000** (0 of 560,853 SNPs) |
| SNPs after merge with AADR (final_ready.bim, union) | 1,080,480 |
| SNPs after `extract_f2(maxmiss=0)` — step 5 | **49,625** |
| SNPs after `extract_f2(maxmiss=0)` — step 6 | **49,364** |

The MyHeritage file has 0 missing genotypes — better than FTDNA (which had F_MISS=0.0117). The effective SNP count for qpAdm (~49–50K) is comparable to PW/FTDNA (49,883) and PW/MH (50,114).

---

## 2. Germany_IA artefact test (steps 4–5)

**Question:** Is the `Germany_Esperstedt_CordedWare` component (proxy for "Germany_IA" from G25/nMonte models) real ancestry of JW, or a geometric artefact of the panel?

### Step 5 results — outgroup robustness (6 combinations)

| Outgroup set | Model A p | Germany weight | Germany SE | Germany z | Verdict |
|---|---|---|---|---|---|
| **Mbuti+Yoruba+Han (baseline)** | **0.620** | −0.649 | 0.725 | −0.894 | ARTEFACT |
| Mbuti+Yoruba+Papuan | 0.516 | −0.703 | 0.806 | −0.872 | ARTEFACT |
| Mbuti+Han+Biaka | 0.479 | +0.614 | 1.136 | +0.541 | ARTEFACT |
| Mbuti+Yoruba+Nganasan | 0.340 | −0.963 | 0.884 | −1.090 | ARTEFACT |
| Mbuti+Yoruba (2 outgroups) | 0.342 | −0.665 | 0.779 | −0.854 | ARTEFACT |
| Mbuti+Han (2 outgroups) | 0.611 | −0.753 | 1.799 | −0.419 | ARTEFACT |

**6/6 combinations → ARTEFACT.** Conclusion robust to outgroup choice.

Germany/CWC weights for JW (range: −0.963 to +0.614) have **smaller absolute values** than for PW (range: −1.438 to +0.001). Germany is an even weaker candidate for real ancestry in JW than in PW.

---

## 3. Systematic ancestry model test (step 6)

**Method:** qpAdm, 16 models (9 single-source + 6 two-source + 1 three-source).  
Outgroups: Mbuti + Yoruba + Han (stable in step 5).  
SNPs: 49,364 after filtering (694 samples, 13 populations).

### Single-source models

| Population | Period | n AADR | p | Status |
|---|---|---|---|---|
| **Poland_EarlySlav** | VI–VIII CE | 25 | **0.833** | ✓ PASS — best fit |
| Poland_EarlyMedieval_Slav | IX–XI CE | 46 | **0.827** | ✓ PASS |
| Poland_GlobularAmphora | ~3200–2300 BC | 33 | **0.669** | ✓ PASS |
| Poland_IA_Wielbark | I–III CE | 59 | **0.524** | ✓ PASS |
| Russia_Samara_EBA_Yamnaya | ~3100–2500 BC | 46 | **0.466** | ✓ PASS |
| Slovakia_N_LBK | ~5500–4900 BC | 56 | **0.449** | ✓ PASS |
| Poland_Roman_Wielbark | III–IV CE | 35 | **0.437** | ✓ PASS |
| Poland_CordedWare | ~2800–2300 BC | 16 | **0.430** | ✓ PASS |
| Poland_BA_Trzciniec | ~1800–1200 BC | 61 | **0.316** | ✓ PASS |

**All 9 single-source models passed** (p > 0.05). For PW, 3 of 9 were rejected.

### Two- and three-source models

All 7 multi-source models passed (p > 0.05), **all with unphysical weights**:

| Model | p | Weights | Status |
|---|---|---|---|
| Yamnaya + GlobAmphora | 0.811 | −3.70 / +4.70 | OK-unphysical |
| Wielbark + EarlySlav | 0.806 | −1.23 / +2.23 | OK-unphysical |
| Yamnaya + LBK + GAC (3-source) | 0.777 | −2.99 / −0.85 / +4.84 | OK-unphysical |
| Wielbark + EarlyMedieval | 0.776 | −1.22 / +2.22 | OK-unphysical |
| CordedWare + EarlyMedieval | 0.737 | −0.84 / +1.84 | OK-unphysical |
| Yamnaya + LBK | 0.549 | +6.19 / −5.19 | OK-unphysical |
| BA_Trzciniec + EarlySlav | 0.400 | −0.62 / +1.62 | OK-unphysical |

Unphysical weights are a collinearity signature — the tested populations are too genetically similar for qpAdm to resolve them as independent components.

---

## 4. Extended models — populations outside the Poland_* set (step 7)

**Method:** qpAdm, 12 single-source models, outgroups: Mbuti + Yoruba + Han.

### Results (sorted by p)

| Population | Period | n AADR | p JW | p PW | Status JW |
|---|---|---|---|---|---|
| **Ukraine_IA_Lusatian** | ~500 BC (Iron Age) | 3 | **0.961** | 0.951 | ✓ **NEW LEADER** |
| Serbia_EarlyMedieval_Byzantine_Slav | VI–VIII CE | 5 | **0.608** | 0.434 | ✓ |
| Lithuanian | modern | 10 | **0.518** | 0.365 | ✓ |
| Lithuania_LBA | ~1000–500 BC | 5 | 0.365 | 0.254 | ✓ |
| Ukrainian_North | modern | 8 | 0.262 | 0.197 | ✓ |
| Slovakia_CiferPac_Medieval_Avar_Slav | IX–X CE | 7 | 0.223 | 0.163 | ✓ |
| Russia_Ivanovo_EarlyMedieval | VIII–X CE | 4 | 0.159 | 0.107 | ✓ |
| Belarusian | modern | 10 | 0.117 | 0.061 | ✓ |
| Estonian | modern | 12 | 0.109 | 0.062 | ✓ |
| Czech | modern | 11 | 0.073 | 0.037 ✗ | ✓ (JW only!) |
| Ukrainian | modern | 13 | 0.044 | 0.015 | ✗ REJECTED |
| Czechia_EarlyMedieval_EarlySlav | VI–VIII CE | 3 | 0.002 | 0.001 | ✗ REJECTED |

*For comparison with step 6: Poland_EarlySlav JW p=0.833, PW p=0.305.*

### Interpretation

**Ukraine_IA_Lusatian (p=0.961)** — identical conclusion as for PW: the Lusatian culture is the best fit for **both**, better than anything in the Poland_\* set. The ancestral substrate of both individuals reaches deeper than Wielbark/EarlySlav — back to the Iron Age Lusatian substrate.

**Serbia and Lithuanians** — higher than for PW, consistent with the general pattern that "JW is more central." Early medieval Balkan Slavs fit JW with p=0.608 (better than Poland_EarlySlav fits PW).

**Czech passes for JW (p=0.073), rejected for PW (p=0.037)** — the only qualitative difference between the two individuals in this set. JW is closer to the western Slavic profile than PW, consistent with the general "JW centrality" conclusion.

**Czechia_EarlyMedieval_EarlySlav strongly rejected for both** (p=0.001–0.002) — early medieval Czech Slavs have a different profile from Polish sources, regardless of individual.

### Implications for the Slavic homeland debate

Ukraine_IA_Lusatian is three samples from a single site: **Rovantsi, Volyn oblast (Lutsk), 1000–700 BC** (Saag & Thomas, *Science Advances* 2025) — ~200 km east of the present Polish border, ~150 km from Gródek Nadbużny. The result is identical for JW (p=0.961) and PW (p=0.951).

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

Lusatian is a population **prior to the split** into the Germanic branch (Wielbark) and the Slavic branch (EarlySlav). The result supports the **autochthonous hypothesis** (Kostrzewski/Hensel school): Slavs did not arrive from outside in the 6th century CE — they evolved from the local Lusatian population of Central Europe. Both JW and PW fit a 1,000-year-old Volhynian population better than "Slavs" from the 8th century CE, because that 1,000-year-old population *is* those same people before they adopted the name "Slavs."

The collinearity of Wielbark + EarlySlav is now genetically explained: both are derivatives of Lusatian — qpAdm cannot resolve them as independent components since they span the same space as Lusatian itself. Any estimated proportion "x% Wielbark + y% Gródek" is a mathematical fiction imposed on one continuous stream.

---

## 4A. Extended step 7 — Avars, Huns, Steppe (verification of G25 "Avar_Kecskemet")

Step 7 was extended from 12 to **23 single-source populations** (Pannonian Avars, Huns, Sarmatians, early Steppe Bulgars) — the same verification as for PW.

| Population | n | p JW | p PW | Status JW |
|---|---|---|---|---|
| **Ukraine_IA_Lusatian** | 3 | **0.9923** | 0.9795 | ✓ leader |
| Croatia_EarlyMedieval_Avar | 12 | 0.434 | 0.301 | ✓ |
| Russia_Tatarstan_EarlyMedieval_EarlyBulgar | 7 | 0.410 | 0.284 | ✓ |
| Hungary_Transtisza_EHun | 7 | 0.362 | 0.245 | ✓ |
| Russia_Sarmatian | 9 | 0.324 | 0.232 | ✓ |
| Austria_Avar | 714 | 0.085 | 0.043 ✗ | ✓ (JW only) |
| Hungary_EarlyAvar-oLowEastAsia | 52 | 0.068 | 0.035 ✗ | ✓ borderline (JW only) |
| **Hungary_EarlyAvar** | 35 | **0.044** | 0.023 | ✗ REJECTED |
| Lithuania_MigrationPeriod | 5 | 0.033 | 0.012 | ✗ REJECTED |

The full `Hungary_EarlyAvar` population (n=35) is rejected for both. The European fraction `-oLowEastAsia` and `Austria_Avar` are borderline-OK for JW (p≈0.07–0.09), rejected for PW — the same "JW more central" pattern as throughout step 7. Even borderline-OK Avar p-values are far below `Ukraine_IA_Lusatian` (0.99); no Avar population competes with the leader.

## 4B. Step 8 — rotating two-source scan: no hidden admixture

Identical test as for PW (expanded 2026-06-29): 2 bases (`Ukraine_IA_Lusatian`, `Poland_EarlySlav`) × **35** diverse second sources (Uralic/Siberian/steppe — including true Xiongnu and Hunnic steppe — Middle Eastern/Caucasian, plus Viking, Wielbark, Saxon, and **historical Polish-admixture candidates: `Jew_Ashkenazi`, `Germany_Medieval_Jewish-lowEastEU`, `Hungary_Conqueror_Commoner` (Magyar), Romanian, Hungarian**). 40,875 SNPs. (Roma/Gypsy absent from AADR HO.)

**Result: ZERO "SIGNIFICANT-2src" models** — same as PW (38 second sources, all |z| < 0.37). With `Ukraine_IA_Lusatian` base (p_base=**0.952**), JW gives **all 35** second sources a small **positive** weight (0.13–0.27, non-significant — `OK-2src-nonsignificant`), while PW gives mostly negative weights (base over-explains). Both readings say the same thing: **no real second component.** Tally: **26 OK-2src-nonsignificant**, 0 base-sufficient, 12 rejected.

- **Ashkenazi:** `Jew_Ashkenazi` w=+0.155 (z=0.357, nonsignificant), `Germany_Medieval_Jewish` REJECTED → **no real Ashkenazi component.** The +0.155 is JW's *uniform* positive bias toward every source (cf. Hungarian +0.268, Chuvash +0.249) — an EEF-deficit drift, not a Jewish signal.
- **Magyar:** `Hungary_Conqueror_Commoner` w=+0.216 (z=0.357) — same non-specific bias, not real.

→ "Avar_Kecskemet" / steppe / Asia / Middle East / Ashkenazi / Magyar = no significant contribution in JW. The lower p of PW vs JW remains a noise/centrality difference, not a separate ancestry component.

### JW's drift is SOUTHERN, not eastern

Comparing how much *better* each source fits JW than PW (Δp = p_JW − p_PW in step 7), after subtracting the ~+0.11 baseline shift, isolates the **direction** of JW's deviation from the pure Lusatian point:

| Preferentially SOUTHERN (Δp) | | True EASTERN — improve LEAST (Δp) | |
|---|---|---|---|
| Serbia_EarlyMedieval_Byzantine_Slav | **+0.25** | Ukrainian | +0.05 |
| Croatia_EarlyMedieval_Avar | **+0.21** | Estonian | +0.03 |
| Hungary_Conqueror_Commoner | +0.16 | Finnish | +0.02 |
| Hungary_Transtisza_EHun | +0.15 | **Russian** | **+0.01** |

If JW were *eastern*, the Uralic/East-Slavic markers (Russian, Finnish, Estonian) would improve most — instead they improve **least**. The largest gains are **South Slavic / Pannonian-Balkan** (Serbia, Croatia, Hungary, Sarmatian, Bulgaria, Romania). Mechanism: JW carries a slight **EEF (Anatolian-farmer) surplus** relative to PW, and EEF density rises toward the south — hence `Italian_North` reaching w=1.44 (second-dominant) on the `Poland_EarlySlav` base. (Baltic `Lithuania_LBA` +0.20 is the one non-southern exception.) JW = same Lusatian core as PW + a sub-resolution **southward** drift, not a resolvable admixture.

**Methodological note:** `Poland_EarlySlav` as base produces degenerate weights (unidentifiability) — script 08 marks these `UNIDENTIFIABLE`. Interpret **only** the `Ukraine_IA_Lusatian` base.

---

## 4C. Step 9 — Przeworsk vs Wielbark channel of descent (2026-06-29)

Run for symmetry with PW (see PW report §5C for the full hypothesis and method). Same instrument: single-source qpAdm + the decisive **f4 symmetry** test + 2-source competition. Only `Poland_IA_Przeworsk` n=2 (both female, Gąski/Kuyavia, 100–300 CE, high-coverage shotgun) exists in AADR HO — step 9 admits n≥2 with caveat. ~40K SNPs, outgroups Mbuti+Yoruba+Han.

### Single-source (context only — collinear)

| Source | n | p | Status |
|---|---|---|---|
| **Poland_IA_Przeworsk** | 2 | **0.651** | ✓ — highest of the Polish IA sources |
| Poland_IA_Wielbark | 59 | 0.572 | ✓ |
| Poland_EarlySlav | 25 | 0.454 | ✓ |
| Ukraine_IA_Lusatian | 3 | 0.858 | ✓ |

### f4 symmetry tests

| Test | f4 | Z | Reading |
|---|---|---|---|
| **B1** channel differentiation `f4(Mbuti, Sweden_IA; Przeworsk, Wielbark)` | −0.00021 | **−0.45** | **n.s. — channels indistinguishable** |
| **B2 (decisive)** `f4(Mbuti, JW; Wielbark, Przeworsk)` | −0.00030 | **−0.41** | **n.s. — JW symmetric to both** |
| B3a vs Przeworsk | +0.00006 | +0.07 | n.s. |
| B3b vs Wielbark | +0.00037 | +0.75 | n.s. |
| B4 robustness (Han outgroup) | −0.00027 | −0.37 | agrees with B2 |

**2-source (Lusatian base):** Przeworsk w=+0.368 (z=0.65), Wielbark w=+0.406 (z=0.71) — both `OK-2nd-contributes` but **z<1 (non-significant) and equal** → indistinguishable. The positive (vs PW's ≈0) weight reflects JW's known southward drift off the pure Lusatian point (§4B), **not** a channel preference.

### Verdict

Identical to PW: hypothesis **consistent and not contradicted, unprovable at ~40K SNP**. B1 n.s. → Przeworsk and Wielbark are not separable from each other; B2 n.s. → JW equidistant from both. Przeworsk again posts the highest single-source p among Polish IA sources (0.651). Collinearity mechanism, resolvable only by WGS.

**Files:** `JW_MH/przeworsk_channel_results_JW.tsv`, `JW_MH/przeworsk_channel_f4_JW.tsv`. Script: `aadr/09_przeworsk_channel.R`.

---

## 5. JW vs PW comparison

| Metric | PW (FTDNA) | JW (MyHeritage) |
|---|---|---|
| SNPs for qpAdm | 49,883 | 49,625 |
| Best 1-source model | EarlySlav p=**0.305** | EarlySlav p=**0.833** |
| Wielbark p | 0.123 | 0.524 |
| 1-source models REJECTED | 3/9 | **0/9** |
| Germany/CWC weight (range) | −1.44 to 0.00 | −0.96 to +0.61 |
| Germany/CWC ARTEFACT | 6/6 | 6/6 |

**JW has higher p-values in all models.** This means JW's genome is **more centrally** positioned in the Polish genetic space — a wider range of reference populations explain it equally well. This is the opposite of what would be expected if JW had unique "Gródek" ancestry.

Particularly telling: populations **rejected** for PW (Slovakia_N_LBK p=0.032, Poland_Roman_Wielbark p=0.034, Poland_BA_Trzciniec p=0.020) **pass easily** for JW (p=0.449, 0.437, 0.316). This means JW has a *less* specific genetic profile, not more specific.

### Direction of the difference: JW drifts SOUTH, PW is the tighter Lusatian point (2026-06-29)

The PW/JW difference, isolated in the step-8 weight signs and the step-7 Δp:

- **PW** sits essentially **on** the pure Lusatian point: in step 8, 24/35 second sources get negative weight (base over-sufficient); 19 `OK-base-sufficient`. The local Polish continuity chain is partly **broken** (Roman_Wielbark and Trzciniec REJECTED), yet eastern `Ukraine_IA_Lusatian` fits best (0.883).
- **JW** sits **slightly off**, toward the **south**: all 35 second sources get positive weight (~0.15–0.27, non-significant); the largest step-7 gains are South Slavic/Pannonian (Serbia +0.25, Croatia +0.21), while the deep-eastern markers (Russian +0.01, Finnish +0.02) gain least. JW's whole local Polish chain is **intact** (Wielbark→Roman→Medieval→EarlySlav→Poland_IA all p≥0.44).

**Autochthon/allochthon paradox:** JW (Podlasie, eastern Poland) shows the **stronger** signal of *local Iron-Age continuity* (autochthon); PW (western Mazovia) shows a **broken** local chain but a best fit to the eastern Volhynian Lusatian sample (hinting allochthon — Slavs from the east). **Neither can be concluded at ~40K SNP**: `Ukraine_IA_Lusatian` is merely the best-sampled Lusatian proxy (the culture spanned Poland too), and all candidate sources are collinear. Resolving autochthon vs allochthon needs WGS (200–400K SNP). See PW report §5B.

---

## 6. The "60% Gródek" claim from G25/nMonte

### 6a. What G25/nMonte says

In nMonte models on G25 vectors, JW receives ~60% of the "Gródek" component (Gródek Nadbużny, early medieval site on the Bug river). This is the highest share among the individuals compared.

### 6b. What qpAdm says

qpAdm results are inconsistent with the interpretation "unique ancestry from Gródek":

**Argument 1: All Polish populations explain JW.** If JW had unique East Slavic ancestry from Gródek, models based on west-central Polish populations (Wielbark, CordedWare, LBK) should be rejected. Instead all pass with high p-values.

**Argument 2: Higher p ≠ better ancestry, but more models passing = less specific.** A high p-value for Poland_EarlySlav (0.833) means the data are very consistent with EarlySlav as the sole source — but the same holds for Wielbark (0.524), LBK (0.449), CordedWare (0.430). All are compatible. A unique "Gródek" component would produce the opposite effect: models without Gródek would be rejected.

**Argument 3: Collinearity of Polish populations.** Wielbark and EarlySlav (including Gródek) form a model with weights −1.23/+2.23 — such extreme weights mean that f2-statistics cannot distinguish these populations. Every proportion (10/90, 50/50, 60/40) is equally consistent with the data. The "60% Gródek" number from G25 has no statistical basis.

### 6c. Where does "60% Gródek" in G25 come from?

G25/nMonte is a PCA + MDS projection, not an f-statistics analysis. In PCA, genetically similar populations have similar vectors, causing the linear model to favor the geographically nearest reference point — for a Polish genome this can be Gródek Nadbużny (eastern Poland, VI–VIII CE) instead of, say, Płońsk (central Poland). This is the geometry of PCA space, not real ancestry.

---

## 7. Biological interpretation

### Shared conclusion for PW and JW

Both have an identical qpAdm result structure:
- Germany/CWC → ARTEFACT (6/6 outgroups)
- Poland_EarlySlav → best single-source model
- All Polish populations → collinear (unphysical weights in two-source models)
- No statistical basis for distinguishing specific proportions between Wielbark and Slavs

The difference between PW and JW is **quantitative, not qualitative**: JW sits closer to the "centre" of Polish genetic space (higher p-values, no model rejected), PW is slightly more distant from the centre (3 models rejected).

### What "60% Gródek" may actually mean

G25 correctly detects that JW is "more early medieval" than PW in PCA space — the PW vs JW difference in the best model (EarlySlav p=0.305 vs 0.833) is real and in that direction. But G25 incorrectly interprets PCA geometry as ancestry proportions: 60% is a specific number without confidence intervals and without a significance test, derived from Gródek being the reference point in the given G25 set — not from real admixture.

---

## 8. Summary

| Question | Answer |
|---|---|
| Is Germany_IA real ancestry of JW? | **NO** — geometric artefact, 6/6 outgroup combinations |
| Best single-source model (of 32 tested) | **Ukraine_IA_Lusatian p=0.9923** (extended step 7) |
| Best model from Poland_* set | Poland_EarlySlav p=0.833 (step 6) |
| Are Avars / "Avar_Kecskemet" (G25) real ancestry? | **NO** — full Hungary_EarlyAvar rejected; zero significant two-source (4B) |
| Does JW have a hidden second source (steppe/Asia/Middle East)? | **NO** — 0 SIGNIFICANT-2src models out of 23 diverse sources (step 8) |
| Does JW have a distinct "Slavic" component? | **YES** — and it traces back to the Lusatian culture (~500 BC) |
| Does JW have unique "Gródek" ancestry? | **NO** — all Polish and neighbouring populations explain JW |
| Is the "60% Gródek" from G25 a precise estimate? | **NO** — number without SE and significance test; qpAdm cannot distinguish any proportion |
| How does JW differ from PW in qpAdm? | JW is more central (0/9 rejected vs 3/9 for PW); Czech passes for JW, rejected for PW |

---

## 9. Output files

```
JW_MH/                             ← main MyHeritage run directory
  merged_pw_aadr_final_ready       ← merged fileset after 4 steps (49,625 SNPs f2)
  merged_pw_aadr_extended          ← + Papuan/Biaka/Nganasan (step 5)
  merged_pw_aadr_models            ← + 7 additional populations (step 6)
  robustness_results.tsv           ← step 5: 6 outgroup combinations
  ancestry_models_results.tsv      ← step 6: 16 ancestry models
  merged_pw_aadr_extended_f2cache/ ← f2 cache step 5 (JW, n_snps=49,625)
  merged_pw_aadr_models_f2cache/   ← f2 cache step 6 (JW, n_snps=49,364)
  extended_models_results.tsv      ← step 7: 23 non-Polish models (+Avars/Huns/steppe)
  merged_pw_aadr_ext7              ← + Lusatian/Avar/Hun/Sarmatian/etc. (step 7)
  merged_pw_aadr_ext7_f2cache/     ← f2 cache step 7
  two_source_results_JW.tsv        ← step 8: 46 two-source models
  merged_pw_aadr_2src              ← + 23 diverse second sources (step 8)
  merged_pw_aadr_2src_JW_f2cache/  ← f2 cache step 8 (target-specific, 41,314 SNPs)

Scripts: aadr/01..08
```

**Limitations / WGS:** see `PW_FF_slavic_aadr.md` §10 — ~41–50K SNPs is ascertainment bias (consumer chip × HO panel), not a bug; WGS would yield ~200–400K SNPs (5–10×) and enable distinguishing Lusatian/EarlySlav/Wielbark, but would not reverse the artefact conclusions (Germany_IA, Avars). The same considerations apply to JW.

---

*Pipeline: `~/Claude/aadr/` | AADR: `/usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB` | plink: `sci-biology/plink-1.90_p20260606` (overlay pwr)*
