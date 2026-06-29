#!/usr/bin/env Rscript
# ==============================================================================
# 06_ancestry_models.R -- systematic ancestry model test
#
# Tests a set of 1-source and 2-source qpAdm models to answer:
#   - Which single populations (beyond Wielbark) also explain the target?
#   - Which 2-source models pass and with what weights?
#   - How deep in time does the explanation reach (proximal vs distal)?
#
# Populations tested as LEFT (source):
#   Proximal:  Poland_IA (n=130, Polish Iron Age = Lusatian horizon in Poland),
#              Poland_IA_Wielbark, Poland_Roman_Wielbark,
#              Poland_EarlyMedieval_Slav, Poland_EarlySlav,
#              Poland_BA_Trzciniec, Poland_CordedWare, Poland_GlobularAmphora
#   Distal:    Russia_Samara_EBA_Yamnaya, Slovakia_N_LBK
#
# Outgroups (RIGHT): Mbuti + Yoruba + Han -- proven stable in step 5.
#
# Usage (run from same dir as steps 4 and 5):
#   Rscript aadr/06_ancestry_models.R \
#     ./FF_PW/merged_pw_aadr_extended \
#     /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
#     PW
#
# Arg 1: prefix from step 5 (merged_pw_aadr_extended) -- already has Mbuti/Yoruba/Han
# Arg 2: aadr prefix (read-only; new populations will be added)
# Arg 3: target ID (default PW)
#
# Output (in outdir):
#   ancestry_models_results.tsv
#   <outdir>/merged_pw_aadr_models.{bed,bim,fam}  (idempotent)
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

# Wide output -- print.data.frame otherwise wraps columns down at width=80.
# Honour $COLUMNS if the terminal exports it, otherwise use a wide default.
options(width = max(200L, as.integer(Sys.getenv("COLUMNS", "0")), na.rm = TRUE))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 06_ancestry_models.R <ext_prefix> <aadr_prefix> [target]")
}

ext_prefix  <- args[1]   # from step 5
aadr_prefix <- args[2]
target_id   <- if (length(args) >= 3) args[3] else "PW"
outdir      <- dirname(ext_prefix)

cat("================================================================================\n")
cat("  STEP 6: Systematic ancestry model test\n")
cat("================================================================================\n\n")
cat(sprintf("  Base prefix   : %s\n", ext_prefix))
cat(sprintf("  AADR prefix   : %s\n", aadr_prefix))
cat(sprintf("  Target        : %s\n\n", target_id))

# ---- Outgroups (RIGHT) -- same set as step 5 baseline ----
right_pops <- c("Mbuti", "Yoruba", "Han")

# ---- New populations to extract from AADR ----
# Poland_IA_Wielbark, Poland_EarlyMedieval_Slav -- already in file from step 2
# Poland_GlobularAmphora -- Poland GAC (Globular Amphora Culture)
new_pops <- c(
  "Poland_IA",                   # n=130, Polish Iron Age aggregate (Lusatian horizon)
  "Poland_EarlySlav",
  "Poland_Roman_Wielbark",
  "Poland_BA_Trzciniec",
  "Poland_CordedWare",
  "Poland_GlobularAmphora",
  "Russia_Samara_EBA_Yamnaya",   # Yamnaya -- Pontic steppe representative
  "Slovakia_N_LBK"               # Early Farmer -- Neolithic farmer representative
)

# ---- Models to test ----
# Format: list(name = c(left_pop1, ...) )
models_1src <- list(
  "1src: Poland_IA"                   = c("Poland_IA"),
  "1src: Poland_IA_Wielbark"          = c("Poland_IA_Wielbark"),
  "1src: Poland_Roman_Wielbark"       = c("Poland_Roman_Wielbark"),
  "1src: Poland_EarlyMedieval_Slav"   = c("Poland_EarlyMedieval_Slav"),
  "1src: Poland_EarlySlav"            = c("Poland_EarlySlav"),
  "1src: Poland_BA_Trzciniec"         = c("Poland_BA_Trzciniec"),
  "1src: Poland_CordedWare"           = c("Poland_CordedWare"),
  "1src: Poland_GlobularAmphora"      = c("Poland_GlobularAmphora"),
  "1src: Russia_Samara_EBA_Yamnaya"   = c("Russia_Samara_EBA_Yamnaya"),
  "1src: Slovakia_N_LBK"              = c("Slovakia_N_LBK")
)

models_2src <- list(
  # Proximal mixes: Iron Age + Medieval Slavic
  "2src: Poland_IA + EarlySlav"       = c("Poland_IA",
                                           "Poland_EarlySlav"),
  "2src: Wielbark + EarlySlav"        = c("Poland_IA_Wielbark",
                                           "Poland_EarlySlav"),
  "2src: Wielbark + EarlyMedieval"    = c("Poland_IA_Wielbark",
                                           "Poland_EarlyMedieval_Slav"),
  "2src: BA_Trzciniec + EarlySlav"    = c("Poland_BA_Trzciniec",
                                           "Poland_EarlySlav"),
  "2src: CordedWare + EarlyMedieval"  = c("Poland_CordedWare",
                                           "Poland_EarlyMedieval_Slav"),
  # Distal: Steppe + Early Farmer
  "2src: Yamnaya + LBK"               = c("Russia_Samara_EBA_Yamnaya",
                                           "Slovakia_N_LBK"),
  # Distal + GAC (between LBK and Corded Ware)
  "2src: Yamnaya + GlobAmphora"       = c("Russia_Samara_EBA_Yamnaya",
                                           "Poland_GlobularAmphora")
)

models_3src <- list(
  # Classic three-component model for Central Europe
  "3src: Yamnaya + LBK + GAC"         = c("Russia_Samara_EBA_Yamnaya",
                                           "Slovakia_N_LBK",
                                           "Poland_GlobularAmphora")
)

all_models <- c(models_1src, models_2src, models_3src)
all_left_pops <- unique(unlist(all_models))

cat(sprintf("  Models to test: %d\n", length(all_models)))
cat(sprintf("  LEFT populations: %s\n\n", paste(all_left_pops, collapse = ", ")))

# ---- Check which populations are already in the base file ----
fam_ext <- read.table(paste0(ext_prefix, ".fam"),
                      col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                      stringsAsFactors = FALSE)
existing_fids  <- unique(fam_ext$FID)
extra_needed   <- setdiff(c(all_left_pops, right_pops), existing_fids)
extra_needed   <- setdiff(extra_needed, c("Mbuti", "Yoruba", "Han"))  # always present

cat(sprintf("  Already in file : %s\n",
            paste(sort(intersect(all_left_pops, existing_fids)), collapse = ", ")))
cat(sprintf("  Need from AADR  : %s\n\n",
            if (length(extra_needed) == 0) "none" else paste(extra_needed, collapse = ", ")))

# ---- Build models file (idempotent, stale detection) ----
models_prefix <- file.path(outdir, "merged_pw_aadr_models")

# Stale detection: file exists but is missing newly added LEFT populations
# (mirrors steps 7-8). Without this, adding a population to the model list and
# re-running stops at the validation below ("Delete the .bed file") instead of
# rebuilding automatically.
if (file.exists(paste0(models_prefix, ".bed")) && length(extra_needed) > 0) {
  models_fids <- unique(read.table(paste0(models_prefix, ".fam"),
                                   col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                                   stringsAsFactors = FALSE)$FID)
  stale_pops <- setdiff(extra_needed, models_fids)
  if (length(stale_pops) > 0) {
    cat(sprintf("  Models file exists but is missing %d populations -- removing and rebuilding.\n",
                length(stale_pops)))
    cat(sprintf("  Missing: %s\n\n", paste(stale_pops, collapse = ", ")))
    invisible(file.remove(paste0(models_prefix, c(".bed", ".bim", ".fam"))))
  }
}

if (file.exists(paste0(models_prefix, ".bed"))) {
  cat(sprintf("  Models file already exists: %s.bed -- skipping build.\n\n",
              models_prefix))
} else if (length(extra_needed) == 0) {
  cat("  All populations already available -- using file from step 5.\n\n")
  models_prefix <- ext_prefix
} else {
  # ---- Verify labels in AADR ----
  cat("--- Verifying populations in AADR .ind ---\n")
  ind_check <- read.table(paste0(aadr_prefix, ".ind"),
                          col.names = c("ID", "Sex", "Pop"),
                          stringsAsFactors = FALSE)
  missing_pops <- setdiff(extra_needed, unique(ind_check$Pop))
  if (length(missing_pops) > 0) {
    stop("Population(s) not found in AADR .ind: ",
         paste(missing_pops, collapse = ", "))
  }
  n_new <- table(ind_check$Pop[ind_check$Pop %in% extra_needed])
  cat("  Sample sizes:\n"); print(n_new); cat("\n")

  # ---- Read from AADR (memory-safe) ----
  cat(sprintf("--- Reading %d populations from AADR ---\n", length(extra_needed)))
  cat("  (read_packedancestrymap filters before loading -- RAM-safe)\n")
  gd   <- read_packedancestrymap(pref = aadr_prefix, pops = extra_needed, verbose = TRUE)
  geno <- gd$geno
  ind  <- gd$ind    # X1/X2/X3 -- positional access, NOT named columns
  snp  <- gd$snp
  n_snp <- nrow(geno); n_ind <- ncol(geno)
  cat(sprintf("  Loaded: %d SNPs x %d individuals.\n\n", n_snp, n_ind))

  # ---- Write MAP ----
  extra_stem <- file.path(outdir, "extra_model_pops")
  cat("  Writing extra_model_pops.map ...\n")
  write.table(data.frame(snp$CHR, snp$SNP, snp$cm, snp$POS),
              file = paste0(extra_stem, ".map"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # ---- Build allele matrix -> PED ----
  # Identical encoding as step 2 and step 5:
  #   g=2: A1/A1 | g=1: A1/A2 | g=0: A2/A2 | NA: 0/0
  a1 <- snp$A1; a1[nchar(a1) != 1] <- "0"
  a2 <- snp$A2; a2[nchar(a2) != 1] <- "0"

  cat(sprintf("  Building allele matrix (%d SNPs x %d ind, progress every 10%%)...\n",
              n_snp, n_ind))
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
  cat("  Writing extra_model_pops.ped ...\n")
  write.table(cbind(ped_meta, as.data.frame(ped_alleles, stringsAsFactors = FALSE)),
              file = paste0(extra_stem, ".ped"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # ---- PED/MAP -> binary PLINK ----
  cat("  plink2: PED/MAP -> binary ...\n")
  ret <- system2("plink2", c("--pedmap", extra_stem, "--make-bed", "--out", extra_stem))
  if (ret != 0L) stop("plink2 --pedmap failed for extra_model_pops")

  # ---- Merge with file from step 5 ----
  # PITFALL: plink2 --pmerge cannot handle non-concatenating merges -> p-link 1.9
  cat("\n--- Merging new populations into base file ---\n")
  ret <- system2("p-link", c(
    "--bfile",  ext_prefix,
    "--bmerge", paste0(extra_stem, ".bed"),
                paste0(extra_stem, ".bim"),
                paste0(extra_stem, ".fam"),
    "--make-bed",
    "--out",    models_prefix
  ))
  if (ret != 0L) stop("p-link --bmerge failed")

  # ---- Fix FID for new samples ----
  pop_lookup_new <- setNames(ind[[3]], ind[[1]])
  models_fam <- read.table(paste0(models_prefix, ".fam"),
                           col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                           stringsAsFactors = FALSE)
  for (sid in names(pop_lookup_new)) {
    idx <- models_fam$IID == sid
    if (any(idx)) models_fam$FID[idx] <- pop_lookup_new[[sid]]
  }
  write.table(models_fam, paste0(models_prefix, ".fam"),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

  n_pops <- length(unique(models_fam$FID))
  cat(sprintf("  Models file: %d samples, %d populations.\n\n",
              nrow(models_fam), n_pops))
}

# ---- Single extract_f2 pass for all required populations ----
all_pops_needed <- unique(c(target_id, all_left_pops, right_pops))
f2_dir <- file.path(outdir, paste0(basename(models_prefix), "_f2cache"))

cat("--- Computing f2 for all populations (single pass) ---\n")
cat(sprintf("  Populations: %s\n", paste(sort(all_pops_needed), collapse = ", ")))
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)

# Verify all required populations are in the file
fam_models <- read.table(paste0(models_prefix, ".fam"),
                         col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                         stringsAsFactors = FALSE)
fids_present <- unique(fam_models$FID)
missing_from_file <- setdiff(all_pops_needed, fids_present)
if (length(missing_from_file) > 0) {
  stop("Populations missing from models file: ",
       paste(missing_from_file, collapse = ", "),
       "\n  -> Delete the .bed file and re-run to force rebuild.")
}

extract_f2(models_prefix, pops = all_pops_needed, outdir = f2_dir,
           overwrite = TRUE, format = "plink")
cat("\n")

# ---- Run models ----
cat("--- qpAdm tests ---\n\n")
results <- list()

run_model <- function(f2_dir, target, left, right, label) {
  f2 <- tryCatch(
    read_f2(f2_dir, pops = c(target, left, right)),
    error = function(e) {
      cat(sprintf("    ERROR read_f2 [%s]: %s\n", label, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(f2)) {
    return(data.frame(model = label, n_src = length(left),
                      p = NA_real_, weights = NA_character_,
                      verdict = "ERROR", stringsAsFactors = FALSE))
  }

  res <- tryCatch(
    qpadm(f2, target = target, left = left, right = right),
    error = function(e) {
      cat(sprintf("    ERROR qpadm [%s]: %s\n", label, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) {
    return(data.frame(model = label, n_src = length(left),
                      p = NA_real_, weights = NA_character_,
                      verdict = "ERROR", stringsAsFactors = FALSE))
  }

  # p-value from rankdrop (f4rank=0 = model with minimum degrees of freedom)
  p <- res$rankdrop$p[res$rankdrop$f4rank == 0]
  if (length(p) == 0 || is.na(p)) p <- NA_real_

  # Weights: "Pop1=0.XX(z=N.N), Pop2=0.XX(z=N.N)"
  w_df    <- res$weights
  w_parts <- sprintf("%s=%.3f(z=%.1f)",
                     w_df$left,
                     w_df$weight,
                     w_df$z)
  w_str   <- paste(w_parts, collapse = " | ")

  feasible <- all(w_df$weight >= -0.1 & w_df$weight <= 1.1, na.rm = TRUE)
  verdict  <- if (is.na(p)) "ERROR" else if (p > 0.05 && feasible) "OK" else if (p > 0.05) "OK-unphysical" else "REJECTED"

  cat(sprintf("  %-40s p=%.4f  %s  [%s]\n", label,
              ifelse(is.na(p), -1, p), verdict, w_str))

  data.frame(model = label, n_src = length(left),
             p = round(p, 4), weights = w_str,
             verdict = verdict, stringsAsFactors = FALSE)
}

cat("  [--- 1-source models ---]\n")
for (label in names(models_1src)) {
  left <- models_1src[[label]]
  results[[label]] <- run_model(f2_dir, target_id, left, right_pops, label)
}

cat("\n  [--- 2-source models ---]\n")
for (label in names(models_2src)) {
  left <- models_2src[[label]]
  results[[label]] <- run_model(f2_dir, target_id, left, right_pops, label)
}

cat("\n  [--- 3-source models ---]\n")
for (label in names(models_3src)) {
  left <- models_3src[[label]]
  results[[label]] <- run_model(f2_dir, target_id, left, right_pops, label)
}

cat("\n")

# ---- Result table ----
res_df <- do.call(rbind, results)
rownames(res_df) <- NULL

# Sort: non-rejected first, then by p descending, then by n_src
res_ok  <- res_df[!is.na(res_df$p) & res_df$p > 0.05, ]
res_rej <- res_df[is.na(res_df$p) | res_df$p <= 0.05, ]
res_ok  <- res_ok[order(-res_ok$p), ]
res_rej <- res_rej[order(res_rej$p, na.last = TRUE), ]
res_sorted <- rbind(res_ok, res_rej)

cat("================================================================================\n")
cat("  NON-REJECTED MODELS (p > 0.05) -- possible ancestry explanations\n")
cat("================================================================================\n\n")
ok_models <- res_sorted[!is.na(res_sorted$p) & res_sorted$p > 0.05, ]
if (nrow(ok_models) == 0) {
  cat("  No models passed the test (p > 0.05).\n\n")
} else {
  print(ok_models[, c("model", "p", "weights", "verdict")], row.names = FALSE)
  cat("\n")
}

cat("================================================================================\n")
cat("  REJECTED MODELS (p <= 0.05)\n")
cat("================================================================================\n\n")
rej_models <- res_sorted[is.na(res_sorted$p) | res_sorted$p <= 0.05, ]
if (nrow(rej_models) == 0) {
  cat("  All models passed the test.\n\n")
} else {
  print(rej_models[, c("model", "p", "weights", "verdict")], row.names = FALSE)
  cat("\n")
}

tsv_out <- file.path(outdir, "ancestry_models_results.tsv")
write.table(res_sorted, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n\n", tsv_out))

# ---- Interpretation ----
cat("================================================================================\n")
cat("  INTERPRETATION\n")
cat("================================================================================\n\n")
n_ok    <- sum(!is.na(res_df$p) & res_df$p > 0.05, na.rm = TRUE)
n_total <- nrow(res_df)
cat(sprintf("  Non-rejected models: %d / %d\n\n", n_ok, n_total))

best_1src <- res_ok[res_ok$n_src == 1, ]
best_2src <- res_ok[res_ok$n_src == 2, ]
best_3src <- res_ok[res_ok$n_src == 3, ]

if (nrow(best_1src) > 0) {
  cat("  1-SOURCE MODELS THAT PASS:\n")
  for (i in seq_len(nrow(best_1src))) {
    cat(sprintf("    %-42s p=%.4f\n", best_1src$model[i], best_1src$p[i]))
  }
  cat("\n")
  cat("  -> Multiple populations explain the target equally well. Normal when populations\n")
  cat("     are genetically similar (e.g. Wielbark and Early Medieval Slavs share\n")
  cat("     common ancestry). qpAdm at ~50K SNPs cannot distinguish them precisely.\n\n")
} else {
  cat("  No valid 1-source models -- target requires at least 2 sources.\n\n")
}

if (nrow(best_2src) > 0) {
  cat("  2-SOURCE MODELS THAT PASS:\n")
  for (i in seq_len(nrow(best_2src))) {
    cat(sprintf("    %-42s p=%.4f  %s\n",
                best_2src$model[i], best_2src$p[i], best_2src$weights[i]))
  }
  cat("\n")
}

if (nrow(best_3src) > 0) {
  cat("  3-SOURCE MODELS THAT PASS:\n")
  for (i in seq_len(nrow(best_3src))) {
    cat(sprintf("    %-42s p=%.4f  %s\n",
                best_3src$model[i], best_3src$p[i], best_3src$weights[i]))
  }
  cat("\n")
}

cat("  NOTE ON ~50K SNPs:\n")
cat("  Personal genome merged with AADR yields ~50K SNPs (vs 500K+ for aggregate-\n")
cat("  population tests). At this SNP count:\n")
cat("    - qpAdm will NOT distinguish populations that are historically closely related\n")
cat("      (e.g. Wielbark vs Early Medieval Slavs vs Roman Wielbark)\n")
cat("    - A 1-source model 'passes' when one population is a linear approximation\n")
cat("      of the target -- does not imply no other components exist\n")
cat("    - Distal models (Yamnaya + LBK) test deeper ancestry history\n\n")

cat("================================================================================\n")
cat("  END\n")
cat("================================================================================\n")
