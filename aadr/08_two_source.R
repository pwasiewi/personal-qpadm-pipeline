#!/usr/bin/env Rscript
# ==============================================================================
# 08_two_source.R -- rotating 2-source scan: searching for a hidden "surprise"
#
# Step 7 (1-source) showed that PW and JW have an IDENTICAL source ranking,
# differing only in a uniform p-level shift (PW lower). A uniform p drop with
# preserved order = global noise/drift, NOT a different component. But this
# does NOT rule out a hidden ADDITIONAL source that a 1-source model cannot
# detect (1 source => weight always 1.000, z=garbage because se=0).
#
# This script fixes a base source and rotates the SECOND source across
# diverse populations:
#   - Uralic / Siberian forest zone (Nganasan, Selkup, Mansi, Khanty,
#     Mordovian, Chuvash, Bashkir, Tatar_Astrakhan)
#   - Eastern steppe / Huns / Turko-Mongols (Mongol, Yakut, Altaian,
#     Kazakh, Russia_Buryatia_XiongnuPeriod, Kazakhstan_Berel_Hunnic)
#   - Middle East / Caucasus / south (Druze, Armenian, Iranian, Georgian,
#     Adygei, Lezgin, Turkish, Greek, Italian_North)
#   - Viking Age / Scandinavian (Sweden_Viking, Denmark_Viking, Norway_Viking,
#     Iceland_Viking) -- hidden Norse component
#   - Wielbark / Roman Germanic (Poland_IA_Wielbark) -- hidden Germanic IA
#   - Migration Period Saxon (Germany_Anderten_Medieval_Saxon, England_Saxon)
#   - Historical admixture candidates (Jew_Ashkenazi, Germany_Medieval_Jewish,
#     Hungary_Conqueror_Commoner=Magyar, Romanian, Hungarian) -- hidden minority
#
# SURPRISE DETECTOR: a model is "SIGNIFICANT-2src" only when
#   p > 0.05  AND  second-source weight in (0.05, 0.95)  AND  |z_second| > 2.
# Then the second source is a real additional component (not noise). If such
# a source appears for PW but not JW (or vice versa) -- that IS the difference.
#
# Verdicts:
#   SIGNIFICANT-2src      -- p>0.05, w_second in (0.05,0.95), |z|>2: real component
#   OK-base-sufficient    -- p>0.05, w_second<=0.05: second source absent, base alone works
#   OK-second-dominant    -- p>0.05, w_second>=0.95: base weight collapses, second dominates
#                            (not "base sufficient" -- base is actually superfluous here)
#   OK-2src-nonsignificant-- p>0.05, weights physical, z<2: model fits but second not sig
#   REJECTED              -- p<=0.05: model fails
#   UNIDENTIFIABLE        -- weight outside [-0.5,1.5]: collinear blowup, p inflated,
#                            model meaningless -- do NOT read as "good fit"
#
# NOTE: z=garbage from step 7 does NOT apply to weights here: there the weight
# was fixed=1 (1 source), so se=0. Here weights are properly estimated => z is meaningful.
#
# Usage (from a writable directory, e.g. ~/Claude):
#   Rscript aadr/08_two_source.R \
#     ./FF_PW/merged_pw_aadr_final_ready \
#     /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
#     PW
#   Rscript aadr/08_two_source.R ./JW_MH/merged_pw_aadr_final_ready <aadr> JW
#
# Output (in outdir):
#   two_source_results_<target>.tsv
#   <outdir>/merged_pw_aadr_2src.{bed,bim,fam}  (idempotent, stale detection)
#   target-specific f2cache: merged_pw_aadr_2src_<target>_f2cache/
#     (separate for PW and JW -- no collision on sequential runs)
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

# Wide output -- print.data.frame otherwise wraps columns down at width=80.
# Honour $COLUMNS if the terminal exports it, otherwise use a wide default.
options(width = max(200L, as.integer(Sys.getenv("COLUMNS", "0")), na.rm = TRUE))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 08_two_source.R <merged_prefix> <aadr_prefix> [target]")

merged_prefix <- args[1]
aadr_prefix   <- args[2]
target_id     <- if (length(args) >= 3) args[3] else "PW"
outdir        <- dirname(merged_prefix)

cat("================================================================================\n")
cat("  STEP 8: rotating 2-source scan (search for hidden additional source)\n")
cat("================================================================================\n\n")
cat(sprintf("  Merged prefix : %s\n", merged_prefix))
cat(sprintf("  AADR prefix   : %s\n", aadr_prefix))
cat(sprintf("  Target        : %s\n\n", target_id))

right_pops <- c("Mbuti", "Yoruba", "Han")

# ---- Base sources (anchors) -- rotated as first LEFT ----
base_anchors <- c(
  "Ukraine_IA_Lusatian",   # best 1-source from step 7 (p~0.98)
  "Poland_EarlySlav"       # canonical proximal Slavic anchor
)

# ---- Diverse second sources (n>=3, deliberately NOT tested in step 7) ----
# NOTE: Han/Mbuti/Yoruba are outgroups (RIGHT) -- cannot be a second source.
diverse_seconds <- c(
  # Uralic / Siberian forest zone
  "Nganasan",                      # n=33, Uralic-Siberian archetype
  "Selkup",                        # n=24, Samoyedic
  "Mansi",                         # n=10, Ob-Ugric
  "Khanty",                        # n=10, Ob-Ugric
  "Mordovian",                     # n=32, Uralic Volga
  "Chuvash",                       # n=14, Turkic Volga
  "Bashkir",                       # n=55, Turko-Uralic Ural-Volga
  "Tatar_Astrakhan",               # n=10, Turkic lower Volga
  # Eastern steppe / Huns / Turko-Mongols
  "Mongol",                        # n=68, modern Mongols
  "Yakut",                         # n=42, Siberian Turkic
  "Altaian",                       # n=26, Altaic
  "Kazakh",                        # n=67, Central Asian steppe
  "Russia_Buryatia_XiongnuPeriod", # n=7, TRUE Xiongnu (steppe ancestor of Huns)
  "Kazakhstan_Berel_Hunnic",       # n=11, Hunnic Altaic steppe
  # Middle East / Caucasus / south
  "Druze",                         # n=83, Levantine isolate (Middle East proxy)
  "Armenian",                      # n=15, Caucasian-Anatolian
  "Iranian",                       # n=39, Iranian
  "Georgian",                      # n=24, Caucasian
  "Adygei",                        # n=48, NW Caucasus
  "Lezgin",                        # n=11, NE Caucasus
  "Turkish",                       # n=48, Anatolian (steppe + Middle East)
  "Greek",                         # n=19, Southern European
  "Italian_North",                 # n=61, Southern European (farming gradient)
  # Viking Age / Scandinavian (VIII-XI CE) -- hidden Norse component
  "Sweden_Viking",                 # n=145, Swedish Vikings (largest Scandinavian sample)
  "Denmark_Viking",                # n=74, Danish Vikings
  "Norway_Viking",                 # n=28, Norwegian Vikings
  "Iceland_Viking",                # n=18, isolated Norse (less Slavic gene flow)
  # Wielbark / Roman Germanic (I-IV CE Poland) -- hidden Germanic IA component
  "Poland_IA_Wielbark",           # n=59, Wielbark culture (also tested 1-source in step 7)
  # Migration Period Germanic / Saxon (V-VIII CE)
  "Germany_Anderten_Medieval_Saxon", # n=14, Early Medieval Saxon
  "England_Saxon",                 # n=84, Anglo-Saxon (Germanic Migration Period)
  # Historical admixture candidates for Polish populations -- hidden minority component
  "Jew_Ashkenazi",                 # n=8, modern Ashkenazi Jews (Levantine+S-European)
  "Germany_Medieval_Jewish-lowEastEU", # n=19, medieval Erfurt Ashkenazi
  "Hungary_Conqueror_Commoner",    # n=25, Magyar conquerors (hidden steppe component)
  "Romanian",                      # n=10, modern Romanians (Balkan/Vlach)
  "Hungarian"                      # n=22, modern Hungarians (Pannonian)
)

# ---- Verify populations in AADR + sample sizes (n>=3) ----
cat("--- Verifying populations in AADR .ind ---\n")
ind_all <- read.table(paste0(aadr_prefix, ".ind"),
                      col.names = c("ID", "Sex", "Pop"),
                      stringsAsFactors = FALSE)
n_by_pop <- table(ind_all$Pop)

check_pops <- function(pops, kind) {
  ok <- c()
  for (p in pops) {
    n <- if (p %in% names(n_by_pop)) n_by_pop[[p]] else 0L
    if (n < 3L) {
      cat(sprintf("  SKIP %-32s  n=%d (< 3)\n", p, n))
    } else {
      cat(sprintf("  OK   %-32s  n=%d  [%s]\n", p, n, kind))
      ok <- c(ok, p)
    }
  }
  ok
}
base_ok   <- check_pops(base_anchors,    "base")
second_ok <- check_pops(diverse_seconds, "second")
cat("\n")

if (length(base_ok) == 0L)   stop("No available base sources (n>=3).")
if (length(second_ok) == 0L) stop("No available second sources (n>=3).")

all_left_pops <- unique(c(base_ok, second_ok))

# ---- Which populations need to be read from AADR ----
fam_existing  <- read.table(paste0(merged_prefix, ".fam"),
                            col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                            stringsAsFactors = FALSE)
existing_fids <- unique(fam_existing$FID)
# Han/Mbuti/Yoruba usually already in base file; do not re-read outgroups
extra_needed  <- setdiff(all_left_pops, existing_fids)

cat(sprintf("  Already in file : %s\n",
            paste(sort(intersect(c(all_left_pops, right_pops), existing_fids)),
                  collapse = ", ")))
cat(sprintf("  Need from AADR  : %s\n\n",
            if (length(extra_needed) == 0L) "none"
            else paste(extra_needed, collapse = ", ")))

# ---- Build extended file (idempotent, stale detection) ----
ext_prefix <- file.path(outdir, "merged_pw_aadr_2src")

# Stale detection: file exists but is missing newly added populations
if (file.exists(paste0(ext_prefix, ".bed")) && length(extra_needed) > 0L) {
  ext_fids <- unique(read.table(paste0(ext_prefix, ".fam"),
                                col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                                stringsAsFactors = FALSE)$FID)
  stale_pops <- setdiff(extra_needed, ext_fids)
  if (length(stale_pops) > 0L) {
    cat(sprintf("  2src exists but is missing %d populations -- removing and rebuilding.\n",
                length(stale_pops)))
    cat(sprintf("  Missing: %s\n\n", paste(stale_pops, collapse = ", ")))
    invisible(file.remove(paste0(ext_prefix, c(".bed", ".bim", ".fam"))))
  }
}

if (length(extra_needed) == 0L) {
  cat("  All populations already in file -- using original.\n\n")
  ext_prefix <- merged_prefix

} else if (file.exists(paste0(ext_prefix, ".bed"))) {
  cat(sprintf("  2src file already exists and contains all populations: %s.bed -- skipping.\n\n",
              ext_prefix))

} else {
  cat(sprintf("--- Reading %d populations from AADR ---\n", length(extra_needed)))
  gd   <- read_packedancestrymap(pref = aadr_prefix, pops = extra_needed, verbose = TRUE)
  geno <- gd$geno; ind <- gd$ind; snp <- gd$snp
  n_snp <- nrow(geno); n_ind <- ncol(geno)
  cat(sprintf("  Loaded %d SNPs x %d individuals.\n\n", n_snp, n_ind))

  extra_stem <- file.path(outdir, "extra_2src_pops")

  a1 <- snp$A1; a1[nchar(a1) != 1] <- "0"
  a2 <- snp$A2; a2[nchar(a2) != 1] <- "0"

  cat(sprintf("  Building allele matrix (%d SNPs x %d ind)...\n", n_snp, n_ind))
  ped_alleles <- matrix("0", nrow = n_ind, ncol = 2L * n_snp)
  pb_step <- max(1L, n_snp %/% 10L)
  for (i in seq_len(n_snp)) {
    g <- geno[i, ]
    o1 <- character(length(g)); o2 <- character(length(g))
    is2 <- !is.na(g) & g == 2L; is1 <- !is.na(g) & g == 1L
    is0 <- !is.na(g) & g == 0L; isna <- is.na(g)
    o1[is2] <- a1[i]; o2[is2] <- a1[i]
    o1[is1] <- a1[i]; o2[is1] <- a2[i]
    o1[is0] <- a2[i]; o2[is0] <- a2[i]
    o1[isna] <- "0";  o2[isna] <- "0"
    ped_alleles[, 2L*i-1L] <- o1
    ped_alleles[, 2L*i]    <- o2
    if (i %% pb_step == 0L) cat(sprintf("  %d%%\n", round(100L * i / n_snp)))
  }
  iid <- ind[[1]]
  ped_meta <- data.frame(FID = iid, IID = iid, PID = 0, MID = 0, SEX = 0, PHEN = -9)

  write.table(data.frame(snp$CHR, snp$SNP, snp$cm, snp$POS),
              file = paste0(extra_stem, ".map"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(cbind(ped_meta, as.data.frame(ped_alleles, stringsAsFactors = FALSE)),
              file = paste0(extra_stem, ".ped"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  ret <- system2("plink2", c("--pedmap", extra_stem, "--make-bed", "--out", extra_stem))
  if (ret != 0L) stop("plink2 --pedmap failed for extra_2src_pops")

  cat("\n--- Merging into existing file ---\n")
  ret <- system2("p-link", c(
    "--bfile",   merged_prefix,
    "--bmerge",  paste0(extra_stem, ".bed"),
                 paste0(extra_stem, ".bim"),
                 paste0(extra_stem, ".fam"),
    "--make-bed", "--out", ext_prefix
  ))
  if (ret != 0L) stop("p-link --bmerge failed for 2src")

  # Fix FID: individual IDs -> population labels
  pop_lookup_extra <- setNames(ind[[3]], ind[[1]])
  ext_fam <- read.table(paste0(ext_prefix, ".fam"),
                        col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                        stringsAsFactors = FALSE)
  for (sid in names(pop_lookup_extra)) {
    idx <- ext_fam$IID == sid
    if (any(idx)) ext_fam$FID[idx] <- pop_lookup_extra[[sid]]
  }
  write.table(ext_fam, paste0(ext_prefix, ".fam"),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  2src .fam: %d samples, %d population labels.\n\n",
              nrow(ext_fam), length(unique(ext_fam$FID))))
}

# ---- Validate: all required populations are in the file ----
all_pops_f2 <- unique(c(target_id, base_ok, second_ok, right_pops))
fam_chk <- read.table(paste0(ext_prefix, ".fam"),
                      col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                      stringsAsFactors = FALSE)
missing_in_file <- setdiff(all_pops_f2, unique(fam_chk$FID))
if (length(missing_in_file) > 0L) {
  stop("Populations missing from 2src file: ", paste(missing_in_file, collapse = ", "),
       "\n  -> Delete ", ext_prefix, ".bed and re-run.")
}

# ---- f2 cache (target-specific -> no PW vs JW collision) ----
f2_dir <- file.path(outdir, paste0(basename(ext_prefix), "_", target_id, "_f2cache"))
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)

cat("--- Computing f2 (single pass for all populations) ---\n")
cat(sprintf("  Populations (%d): %s\n", length(all_pops_f2),
            paste(all_pops_f2, collapse = ", ")))
extract_f2(ext_prefix, pops = all_pops_f2, outdir = f2_dir,
           overwrite = TRUE, format = "plink")
cat("\n")

# ---- Helper: single qpAdm call (1 or 2 sources) ----
run_qpadm <- function(left, label) {
  f2 <- tryCatch(read_f2(f2_dir, pops = c(target_id, left, right_pops)),
                 error = function(e) {
                   cat(sprintf("  ERROR f2  [%s]: %s\n", label, conditionMessage(e)))
                   NULL
                 })
  if (is.null(f2)) return(NULL)
  tryCatch(qpadm(f2, target = target_id, left = left, right = right_pops),
           error = function(e) {
             cat(sprintf("  ERROR qpAdm [%s]: %s\n", label, conditionMessage(e)))
             NULL
           })
}

# ---- Reference: p of 1-source models (base alone) ----
cat("--- Reference: base as sole source (p_base) ---\n")
base_p <- setNames(rep(NA_real_, length(base_ok)), base_ok)
for (b in base_ok) {
  m <- run_qpadm(b, sprintf("1src: %s", b))
  if (!is.null(m)) {
    pv <- m$rankdrop$p[m$rankdrop$f4rank == 0]
    base_p[[b]] <- if (length(pv)) pv else NA_real_
    cat(sprintf("  %-28s p_base=%.4f\n", b, base_p[[b]]))
  }
}
cat("\n")

# ---- 2-source scan ----
cat("--- 2-source scan: base + diverse second source ---\n\n")
results <- list()

for (b in base_ok) {
  cat(sprintf("  === BASE: %s (p_base=%.4f) ===\n", b, base_p[[b]]))
  for (s in second_ok) {
    if (s == b) next
    label <- sprintf("%s + %s", b, s)
    m <- run_qpadm(c(b, s), label)
    if (is.null(m)) {
      results[[label]] <- data.frame(
        base = b, second = s, p = NA_real_, p_base = base_p[[b]],
        w_base = NA_real_, w_second = NA_real_, se_second = NA_real_,
        z_second = NA_real_, verdict = "ERROR", stringsAsFactors = FALSE)
      next
    }
    pv     <- m$rankdrop$p[m$rankdrop$f4rank == 0]
    pv     <- if (length(pv)) pv else NA_real_
    w_b    <- m$weights$weight[m$weights$left == b]
    w_s    <- m$weights$weight[m$weights$left == s]
    se_s   <- m$weights$se[m$weights$left == s]
    z_s    <- m$weights$z[m$weights$left == s]

    # Unidentifiability detection: weights blow up outside [0,1] (e.g. 32 / -31)
    # when base and second source are too collinear relative to outgroups -> qpAdm
    # cannot separate them; p may be inflated, model is meaningless.
    # Threshold [-0.5, 1.5] allows mild over-explanation (~-0.25), catches blowup.
    degenerate <- (!is.na(w_s) && (w_s > 1.5 || w_s < -0.5)) ||
                  (!is.na(w_b) && (w_b > 1.5 || w_b < -0.5))

    # Surprise detector: model passes + second source significant and physical
    sig_second <- !is.na(pv) && pv > 0.05 && !degenerate &&
                  !is.na(w_s) && w_s > 0.05 && w_s < 0.95 &&
                  !is.na(z_s) && abs(z_s) > 2
    verdict <- if (is.na(pv)) "ERROR"
               else if (degenerate) "UNIDENTIFIABLE"
               else if (sig_second) "SIGNIFICANT-2src"
               else if (pv > 0.05 && !is.na(w_s) && w_s <= 0.05)
                 "OK-base-sufficient"
               else if (pv > 0.05 && !is.na(w_s) && w_s >= 0.95)
                 "OK-second-dominant"
               else if (pv > 0.05) "OK-2src-nonsignificant"
               else "REJECTED"

    flag <- if (verdict == "SIGNIFICANT-2src") "  <== SURPRISE" else ""
    cat(sprintf("    + %-30s p=%.4f  w(%s)=%.3f w(second)=%.3f z=%.2f  %s%s\n",
                s, ifelse(is.na(pv), -1, pv), substr(b, 1, 10),
                ifelse(is.na(w_b), NA, w_b), ifelse(is.na(w_s), NA, w_s),
                ifelse(is.na(z_s), NA, z_s), verdict, flag))

    results[[label]] <- data.frame(
      base = b, second = s, p = round(pv, 4), p_base = round(base_p[[b]], 4),
      w_base = round(w_b, 3), w_second = round(w_s, 3),
      se_second = round(se_s, 4), z_second = round(z_s, 3),
      verdict = verdict, stringsAsFactors = FALSE)
  }
  cat("\n")
}

# ---- Summary table ----
res_df <- do.call(rbind, results)
rownames(res_df) <- NULL
# Sort in 3 tiers: SIGNIFICANT-2src (0) -> meaningful (1) -> UNIDENTIFIABLE (2),
# within tier by p descending. UNIDENTIFIABLE goes to bottom despite inflated p.
tier <- ifelse(res_df$verdict == "SIGNIFICANT-2src", 0L,
        ifelse(res_df$verdict == "UNIDENTIFIABLE", 2L, 1L))
ord <- order(tier, -res_df$p)
res_df <- res_df[ord, ]

cat("================================================================================\n")
cat("  SUMMARY TABLE (SIGNIFICANT-2src on top, then by p descending)\n")
cat("================================================================================\n\n")
print(res_df[, c("base","second","p","p_base","w_second","z_second","verdict")],
      row.names = FALSE)
cat("\n")

tsv_out <- file.path(outdir, sprintf("two_source_results_%s.tsv", target_id))
write.table(res_df, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n\n", tsv_out))

# ---- Interpretation ----
n_sig <- sum(res_df$verdict == "SIGNIFICANT-2src", na.rm = TRUE)
cat("================================================================================\n")
cat("  INTERPRETATION\n")
cat("================================================================================\n\n")
if (n_sig == 0L) {
  cat(sprintf("  ZERO 'SIGNIFICANT-2src' models for %s.\n", target_id))
  cat("  None of the diverse second sources (Uralic/Siberian/steppe/\n")
  cat("  Middle Eastern/Caucasian) receives a significant weight with a Slavic base.\n")
  cat("  => No hidden 'surprise' -- lower p from step 7 is global noise/drift,\n")
  cat("     not a separate ancestry component. Profile consistent with 1-source model.\n\n")
} else {
  cat(sprintf("  FOUND %d 'SIGNIFICANT-2src' model(s) for %s:\n\n", n_sig, target_id))
  sig <- res_df[res_df$verdict == "SIGNIFICANT-2src", ]
  for (i in seq_len(nrow(sig))) {
    cat(sprintf("    %-26s + %-26s  p=%.4f (p_base=%.4f)  weight=%.3f z=%.2f\n",
                sig$base[i], sig$second[i], sig$p[i], sig$p_base[i],
                sig$w_second[i], sig$z_second[i]))
  }
  cat("\n  => This second source is a real additional component. Compare with the\n")
  cat("     other target (PW vs JW): if it appears for only one -- THAT is the difference.\n\n")
}

cat("  CONTEXT (step 7, 1-source):\n")
cat("    PW and JW: identical ranking, Ukraine_IA_Lusatian #1 (p~0.98-0.99).\n")
cat("    PW uniformly lower p than JW -> hypothesis: technical noise (FTDNA chip) or\n")
cat("    slightly higher drift, NOT a separate component. This scan verifies that.\n")
cat("================================================================================\n")
cat("  END\n")
cat("================================================================================\n")
