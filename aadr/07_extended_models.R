#!/usr/bin/env Rscript
# ==============================================================================
# 07_extended_models.R -- non-Polish 1-source qpAdm models
#
# Tests populations outside the Poland_* set that are geographically or
# culturally proximate:
#   - neighbouring Slavic (Czech, Serbian, Russian, Slovak)
#   - Baltic / Lusatian Iron Age (Lithuania_LBA, Ukraine_IA_Lusatian)
#   - Wielbark / Roman Iron Age Poland (Poland_IA_Wielbark, Poland_Roman_Wielbark)
#   - Scandinavian Iron Age (Sweden_IA, Denmark_IA, Norway_IA) -- Wielbark homeland
#   - modern populations of neighbouring countries (Czech, Lithuanian,
#     Belarusian, Ukrainian, Estonian, Norwegian, Finnish, Russian, Romanian,
#     Moldavian, Hungarian, Bulgarian)
#   - historical admixture candidates (Jew_Ashkenazi, Germany_Medieval_Jewish,
#     Hungary_Conqueror_Commoner=Magyar). NOTE: Roma/Gypsy absent from AADR HO.
#   - Huns in Pannonia (Hungary_EHun, Hungary_Transtisza_EHun)
#   - Avars in the Carpathian Basin (Hungary and Austria, early phase)
#   - Migration Period / Sarmatians (Lithuania, Slovakia, Hungary)
#   - Early Bulgars (Russia_Tatarstan, Volga steppe origin)
#
# 1-source models only -- the goal is to rank sources, not decompose.
# Populations with n<3 in AADR HO are skipped with a note.
#
# G25/nMonte context: PW showed ~7% Lithuania_BA + 20-30% Avar_Kecskemet
# (G25 label) instead of 40% Germany_AI -> test whether AADR Avar / Hun
# populations are real ancestry sources or geometric artefacts.
#
# Usage (run from ~/Claude or any writable dir):
#   Rscript aadr/07_extended_models.R \
#     ./FF_PW/merged_pw_aadr_final_ready \
#     /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
#     PW
#
# Output:
#   extended_models_results.tsv  (in CWD)
#   <outdir>/merged_pw_aadr_ext7.{bed,bim,fam}  (idempotent)
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

# Wide output -- print.data.frame otherwise wraps columns down at width=80.
# Honour $COLUMNS if the terminal exports it, otherwise use a wide default.
options(width = max(200L, as.integer(Sys.getenv("COLUMNS", "0")), na.rm = TRUE))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 07_extended_models.R <merged_prefix> <aadr_prefix> [target]")

merged_prefix <- args[1]
aadr_prefix   <- args[2]
target_id     <- if (length(args) >= 3) args[3] else "PW"
outdir        <- dirname(merged_prefix)

cat("================================================================================\n")
cat("  STEP 7: extended 1-source models (outside the Poland_* set)\n")
cat("================================================================================\n\n")
cat(sprintf("  Merged prefix : %s\n", merged_prefix))
cat(sprintf("  AADR prefix   : %s\n", aadr_prefix))
cat(sprintf("  Target        : %s\n\n", target_id))

right_pops <- c("Mbuti", "Yoruba", "Han")

# ---- Candidate source populations ----
# grouped by category; only those with n>=3 in AADR v66 HO
candidate_pops <- c(
  # Early medieval / archaeological Slavic (non-Polish)
  "Czechia_EarlyMedieval_EarlySlav",      # n=3, VI-VIII CE
  "Slovakia_CiferPac_Medieval_Avar_Slav", # n=7, post-Avar Slavs
  "Serbia_EarlyMedieval_Byzantine_Slav",  # n=5, VI-VIII CE
  "Russia_Ivanovo_EarlyMedieval",         # n=4, early medieval Russian
  # Baltic / Lusatian Iron Age
  "Lithuania_LBA",                        # n=5, LBA ~1000-500 BC
  "Ukraine_IA_Lusatian",                  # n=3, Lusatian IA (predecessor of Slavs)
  # Wielbark culture — Roman Iron Age Poland (I-IV CE, Germanic/Gothic horizon)
  "Poland_IA_Wielbark",                   # n=59, Wielbark IA (Pomerania/Masovia, I-III CE)
  "Poland_Roman_Wielbark",                # n=35, Wielbark Roman period subset
  # Scandinavian Iron Age — putative Wielbark homeland (Jastorf/Gotland connection)
  "Sweden_IA",                            # n=8, Swedish IA aggregate
  "Denmark_IA",                           # n=4, Danish IA aggregate
  "Norway_IA",                            # n=4, Norwegian IA aggregate
  # Modern neighbouring populations
  "Czech",                                # n=11, modern Czechs
  "Lithuanian",                           # n=10, modern Lithuanians
  "Belarusian",                           # n=10, modern Belarusians
  "Ukrainian",                            # n=13, modern Ukrainians
  "Ukrainian_North",                      # n=8, North Ukrainians
  "Estonian",                             # n=12, modern Estonians
  "Norwegian",                            # n=12, modern Norwegians (no Swedish aggregate in AADR v66)
  "Finnish",                              # n=17, modern Finns
  "Russian",                              # n=98, modern Russians (East Slavic)
  "Romanian",                             # n=10, modern Romanians (SE neighbour)
  "Moldavian",                            # n=10, modern Moldavians (E neighbour)
  "Hungarian",                            # n=22, modern Hungarians (Magyar + Pannonian)
  "Bulgarian",                            # n=12, modern Bulgarians (South Slavic)
  # Historical admixture candidates for Polish populations
  "Jew_Ashkenazi",                        # n=8, modern Ashkenazi Jews
  "Germany_Medieval_Jewish-lowEastEU",    # n=19, medieval Erfurt Ashkenazi (low East-EU subset)
  "Hungary_Conqueror_Commoner",           # n=25, Magyar conquerors (steppe-shifted, X CE)
  # Huns in Pannonia (V CE)
  "Hungary_EHun",                         # n=7, early Huns in Hungary ~380-430 CE
  "Hungary_Transtisza_EHun",              # n=7, Trans-Tisza early Huns
  # Avars in the Carpathian Basin (VI-IX CE)
  # oLowEastAsia = outlier subset with LOW East-Asian admixture (more European-like)
  "Hungary_EarlyAvar",                    # n=35, early Avar phase ~568-670 CE
  "Hungary_EarlyAvar-oLowEastAsia",       # n=52, European-leaning early Avars
  "Austria_Avar",                         # n=714, Pannonian Avars across all periods
  "Croatia_EarlyMedieval_Avar",           # n=12, Avar horizon at Dalmatian border
  # Migration Period / Sarmatians (IV-VI CE)
  "Lithuania_MigrationPeriod",            # n=5, Baltic Migration Period
  "Slovakia_MigrationPeriod",             # n=6, Migration Period Slovakia
  "Hungary_Madaras_MigrationPeriod_Sarmatian", # n=11, late Sarmatians in Hungary
  "Russia_Sarmatian",                     # n=9, Pontic steppe Sarmatians
  # Early Bulgars (Volga/Pontic steppe -> lower Danube, VII CE)
  "Russia_Tatarstan_EarlyMedieval_EarlyBulgar" # n=7, pre-Volga-Bulgaria phase
)

# ---- Verify populations exist in AADR + check sample sizes ----
cat("--- Verifying populations in AADR .ind ---\n")
ind_all <- read.table(paste0(aadr_prefix, ".ind"),
                      col.names = c("ID", "Sex", "Pop"),
                      stringsAsFactors = FALSE)
n_by_pop  <- table(ind_all$Pop)
confirmed <- c()
for (p in candidate_pops) {
  n <- if (p %in% names(n_by_pop)) n_by_pop[[p]] else 0L
  if (n < 3L) {
    cat(sprintf("  SKIP %-45s  n=%d (< 3)\n", p, n))
  } else {
    cat(sprintf("  OK   %-45s  n=%d\n", p, n))
    confirmed <- c(confirmed, p)
  }
}
cat("\n")

if (length(confirmed) == 0L) stop("No candidate populations with n>=3 found.")

# ---- Check which pops already in merged file ----
fam_existing <- read.table(paste0(merged_prefix, ".fam"),
                           col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                           stringsAsFactors = FALSE)
existing_fids <- unique(fam_existing$FID)
all_needed    <- unique(c(confirmed, right_pops))
extra_needed  <- setdiff(confirmed, existing_fids)

cat(sprintf("  Already in file : %s\n",
            paste(sort(intersect(all_needed, existing_fids)), collapse = ", ")))
cat(sprintf("  Need from AADR  : %s\n\n",
            if (length(extra_needed) == 0L) "none"
            else paste(extra_needed, collapse = ", ")))

# ---- Build extended merged file (idempotent) ----
ext_prefix <- file.path(outdir, "merged_pw_aadr_ext7")

# Detect stale ext7: exists but is missing newly added populations
if (file.exists(paste0(ext_prefix, ".bed")) && length(extra_needed) > 0L) {
  ext7_fids <- unique(read.table(paste0(ext_prefix, ".fam"),
                                  col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                                  stringsAsFactors = FALSE)$FID)
  stale_pops <- setdiff(extra_needed, ext7_fids)
  if (length(stale_pops) > 0L) {
    cat(sprintf("  Ext7 exists but is missing %d populations -- removing and rebuilding.\n",
                length(stale_pops)))
    cat(sprintf("  Missing: %s\n\n", paste(stale_pops, collapse = ", ")))
    invisible(file.remove(paste0(ext_prefix, c(".bed", ".bim", ".fam"))))
  }
}

if (length(extra_needed) == 0L) {
  cat("  All populations already in file -- using original.\n\n")
  ext_prefix <- merged_prefix

} else if (file.exists(paste0(ext_prefix, ".bed"))) {
  cat(sprintf("  Extended file already exists and contains all populations: %s.bed -- skipping.\n\n", ext_prefix))

} else {
  cat(sprintf("--- Reading %d populations from AADR ---\n", length(extra_needed)))
  gd   <- read_packedancestrymap(pref = aadr_prefix, pops = extra_needed, verbose = TRUE)
  geno <- gd$geno; ind <- gd$ind; snp <- gd$snp
  n_snp <- nrow(geno); n_ind <- ncol(geno)
  cat(sprintf("  Loaded %d SNPs x %d individuals.\n\n", n_snp, n_ind))

  extra_stem <- file.path(outdir, "extra_ext7_pops")

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
  if (ret != 0L) stop("plink2 --pedmap failed for extra_ext7_pops")

  cat("\n--- Merging into existing file ---\n")
  ret <- system2("p-link", c(
    "--bfile",   merged_prefix,
    "--bmerge",  paste0(extra_stem, ".bed"),
                 paste0(extra_stem, ".bim"),
                 paste0(extra_stem, ".fam"),
    "--make-bed", "--out", ext_prefix
  ))
  if (ret != 0L) stop("p-link --bmerge failed for ext7")

  # Fix FID: individual IDs -> population labels
  pop_lookup_extra <- setNames(ind[[3]], ind[[1]])
  ext_fam <- read.table(paste0(ext_prefix, ".fam"),
                        col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                        stringsAsFactors = FALSE)
  for (sid in names(pop_lookup_extra)) {
    idx <- ext_fam$IID == sid
    if (any(idx)) ext_fam$FID[idx] <- pop_lookup_extra[[sid]]
  }
  write.table(ext_fam, paste0(ext_prefix, ".fam"),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Ext7 .fam: %d samples, %d population labels.\n\n",
              nrow(ext_fam), length(unique(ext_fam$FID))))
}

# ---- Build f2 cache ----
all_pops_f2 <- unique(c(target_id, confirmed, right_pops))
f2_dir <- file.path(outdir, paste0(basename(ext_prefix), "_f2cache"))
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)

cat("--- Computing f2 (single pass for all populations) ---\n")
cat(sprintf("  Populations: %s\n", paste(all_pops_f2, collapse = ", ")))
extract_f2(ext_prefix, pops = all_pops_f2, outdir = f2_dir,
           overwrite = TRUE, format = "plink")
cat("\n")

# ---- Run 1-source qpAdm for each candidate ----
cat("--- qpAdm 1-source for candidate populations ---\n\n")
results <- list()

for (src in confirmed) {
  f2 <- tryCatch(
    read_f2(f2_dir, pops = c(target_id, src, right_pops)),
    error = function(e) {
      cat(sprintf("  ERROR f2 dla %-40s : %s\n", src, conditionMessage(e))); NULL
    }
  )
  if (is.null(f2)) {
    results[[src]] <- data.frame(source = src, p = NA_real_,
                                 weight = NA_real_, se = NA_real_, z = NA_real_,
                                 verdict = "ERROR", stringsAsFactors = FALSE)
    next
  }
  m <- tryCatch(
    qpadm(f2, target = target_id, left = src, right = right_pops),
    error = function(e) {
      cat(sprintf("  ERROR qpAdm dla %-40s : %s\n", src, conditionMessage(e))); NULL
    }
  )
  if (is.null(m)) {
    results[[src]] <- data.frame(source = src, p = NA_real_,
                                 weight = NA_real_, se = NA_real_, z = NA_real_,
                                 verdict = "ERROR", stringsAsFactors = FALSE)
    next
  }
  p_val <- m$rankdrop$p[m$rankdrop$f4rank == 0]
  w     <- m$weights$weight[m$weights$left == src]
  se    <- m$weights$se[m$weights$left == src]
  z     <- m$weights$z[m$weights$left == src]

  verdict <- if (is.na(p_val)) "ERROR" else if (p_val >= 0.05) "OK" else "REJECTED"
  cat(sprintf("  %-45s  p=%.4f  weight=%.3f  z=%.2f  => %s\n",
              src, p_val, w, z, verdict))

  results[[src]] <- data.frame(source = src, p = round(p_val, 4),
                               weight = round(w, 3), se = round(se, 3),
                               z = round(z, 3), verdict = verdict,
                               stringsAsFactors = FALSE)
}

# ---- Summary ----
res_df <- do.call(rbind, results)
rownames(res_df) <- NULL
res_df <- res_df[order(-res_df$p), ]

cat("\n================================================================================\n")
cat("  SUMMARY TABLE (sorted by p descending)\n")
cat("================================================================================\n\n")
print(res_df, row.names = FALSE)
cat("\n")

tsv_out <- file.path(outdir, "extended_models_results.tsv")
write.table(res_df, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n\n", tsv_out))

cat("================================================================================\n")
cat("  CONTEXT (Poland_* results from step 6 for comparison):\n")
cat("    Poland_EarlySlav          p ~ 0.20-0.31  (best so far)\n")
cat("    Poland_EarlyMedieval_Slav p ~ 0.14-0.22\n")
cat("    Poland_IA_Wielbark        p ~ 0.09-0.12\n")
cat("    Poland_CordedWare         p ~ 0.09-0.13\n")
cat("    Ukraine_IA_Lusatian       p >> Poland_* (best 1-source in previous run)\n")
cat("\n")
cat("  G25/nMonte CONTEXT (hypotheses to verify with qpAdm):\n")
cat("    ~7%  Lithuania_BA + ~20-30% Avar_Kecskemet (G25) visible in PW\n")
cat("    -> hypothesis: Avars and Baltic_BA are real ancestors, not an artefact\n")
cat("    -> AADR equivalents: Hungary_EarlyAvar* / Austria_Avar / Lithuania_LBA\n")
cat("================================================================================\n")
cat("  END\n")
cat("================================================================================\n")
