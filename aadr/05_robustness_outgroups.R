#!/usr/bin/env Rscript
# ==============================================================================
# 05_robustness_outgroups.R -- Germany_IA artefact test, robustness check
#
# Tests whether the qpAdm conclusion (Germany/CWC is a geometric artefact)
# is stable across multiple right-population (outgroup) combinations.
#
# Adds Papuan, Biaka, Nganasan from AADR to the existing merged PLINK file,
# then runs qpAdm with 6 outgroup sets. Stable conclusion = Germany |z|<2
# and weight<0 in ALL combinations.
#
# NOTE: Onge is NOT in AADR HO panel. Substitutes verified present:
#   Papuan (n=46), Biaka, Nganasan.
#
# Usage (run from same dir as step 4):
#   Rscript aadr/05_robustness_outgroups.R \
#     ./FF_PW/merged_pw_aadr_final_ready \
#     /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
#     PW
#
# Output (in CWD):
#   robustness_results.tsv
#   <outdir>/merged_pw_aadr_extended.{bed,bim,fam}  (idempotent: skip if exists)
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_robustness_outgroups.R <merged_prefix> <aadr_prefix> [target]")
}

merged_prefix <- args[1]
aadr_prefix   <- args[2]
target_id     <- if (length(args) >= 3) args[3] else "PW"
outdir        <- dirname(merged_prefix)

cat("================================================================================\n")
cat("  STEP 5: qpAdm outgroup robustness check\n")
cat("================================================================================\n\n")
cat(sprintf("  Merged prefix : %s\n", merged_prefix))
cat(sprintf("  AADR prefix   : %s\n", aadr_prefix))
cat(sprintf("  Target        : %s\n\n", target_id))

# ---- Fixed left populations ----
left_source <- "Poland_IA_Wielbark"
left_german <- "Germany_Esperstedt_CordedWare"

# Extra outgroups to pull from AADR (not in the original step-2 extraction)
extra_right_pops <- c("Papuan", "Biaka", "Nganasan")

# Right-population combinations to test
right_sets <- list(
  "Mbuti+Yoruba+Han (baseline)"  = c("Mbuti", "Yoruba", "Han"),
  "Mbuti+Yoruba+Papuan"          = c("Mbuti", "Yoruba", "Papuan"),
  "Mbuti+Han+Biaka"              = c("Mbuti", "Han",    "Biaka"),
  "Mbuti+Yoruba+Nganasan"        = c("Mbuti", "Yoruba", "Nganasan"),
  "Mbuti+Yoruba (2 outgroups)"   = c("Mbuti", "Yoruba"),
  "Mbuti+Han (2 outgroups)"      = c("Mbuti", "Han")
)

all_right_pops <- unique(unlist(right_sets))

# ---- Check which extra pops need to be extracted ----
fam_existing <- read.table(paste0(merged_prefix, ".fam"),
                           col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                           stringsAsFactors = FALSE)
existing_fids <- unique(fam_existing$FID)
extra_needed  <- setdiff(extra_right_pops, existing_fids)

cat(sprintf("  Already in merged file : %s\n",
            paste(sort(intersect(all_right_pops, existing_fids)), collapse = ", ")))
cat(sprintf("  Need from AADR         : %s\n\n",
            if (length(extra_needed) == 0) "none" else paste(extra_needed, collapse = ", ")))

# ---- Build extended merged file (idempotent) ----
ext_prefix <- file.path(outdir, "merged_pw_aadr_extended")

if (length(extra_needed) == 0) {
  cat("  All outgroups already present -- using original merged file.\n\n")
  ext_prefix <- merged_prefix

} else if (file.exists(paste0(ext_prefix, ".bed"))) {
  cat(sprintf("  Extended file already exists: %s.bed -- skipping build.\n\n", ext_prefix))

} else {
  # ---- Verify populations exist in AADR .ind before the expensive read ----
  cat("--- Verifying populations in AADR .ind ---\n")
  ind_check <- read.table(paste0(aadr_prefix, ".ind"),
                          col.names = c("ID", "Sex", "Pop"),
                          stringsAsFactors = FALSE)
  missing_extra <- setdiff(extra_needed, unique(ind_check$Pop))
  if (length(missing_extra) > 0) {
    stop("Population(s) not found in AADR .ind: ",
         paste(missing_extra, collapse = ", "),
         "\n  Check exact spelling: grep('keyword', ind$Pop, value=TRUE)")
  }
  n_extra <- table(ind_check$Pop[ind_check$Pop %in% extra_needed])
  cat("  Sample sizes: "); print(n_extra); cat("\n")

  # ---- Extract from AADR (memory-safe: pops= filters before loading) ----
  # PITFALL: packedancestrymap_to_plink() loads ALL 27594 samples (OOM).
  # read_packedancestrymap(pops=...) filters first -- always use this.
  cat(sprintf("--- Reading %s from AADR ---\n", paste(extra_needed, collapse = ", ")))
  gd   <- read_packedancestrymap(pref = aadr_prefix, pops = extra_needed, verbose = TRUE)
  geno <- gd$geno   # SNPs x individuals (EIGENSTRAT orientation)
  ind  <- gd$ind    # columns X1/X2/X3 -- PITFALL: NOT ID/Sex/Pop (positional only)
  snp  <- gd$snp
  n_snp <- nrow(geno)
  n_ind <- ncol(geno)
  cat(sprintf("  Loaded %d SNPs x %d individuals.\n\n", n_snp, n_ind))

  # ---- Write MAP ----
  extra_stem <- file.path(outdir, "extra_pops")
  cat("  Writing extra_pops.map ...\n")
  write.table(data.frame(snp$CHR, snp$SNP, snp$cm, snp$POS),
              file = paste0(extra_stem, ".map"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # ---- Build allele matrix and write PED ----
  # EIGENSTRAT dosage -> PED allele pairs (same encoding as step 2):
  #   g=2: A1/A1 (homozygous ref)
  #   g=1: A1/A2 (het)
  #   g=0: A2/A2 (homozygous alt)
  #   NA : 0/0   (missing)
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
  cat("  Writing extra_pops.ped ...\n")
  write.table(cbind(ped_meta, as.data.frame(ped_alleles, stringsAsFactors = FALSE)),
              file = paste0(extra_stem, ".ped"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # ---- PED/MAP -> binary PLINK ----
  cat("  Converting extra_pops PED/MAP -> binary PLINK ...\n")
  ret <- system2("plink2", c("--pedmap", extra_stem, "--make-bed", "--out", extra_stem))
  if (ret != 0L) stop("plink2 --pedmap failed for extra_pops")

  # ---- Merge extra_pops into existing merged file ----
  # PITFALL: plink2 --pmerge cannot handle non-concatenating merges.
  # Use p-link (plink 1.9) --bmerge -- same as step 3e.
  cat("\n--- Merging extra populations into existing merged file ---\n")
  ret <- system2("p-link", c(
    "--bfile",   merged_prefix,
    "--bmerge",  paste0(extra_stem, ".bed"),
                 paste0(extra_stem, ".bim"),
                 paste0(extra_stem, ".fam"),
    "--make-bed",
    "--out",     ext_prefix
  ))
  if (ret != 0L) stop("p-link --bmerge failed (check merge_probe log for allele conflicts)")

  # ---- Fix FID: set population labels for extra samples ----
  # After bmerge, extra samples have FID=IID (individual ID), not population.
  # extract_f2/qpAdm match by FID, so FID must equal the population name.
  pop_lookup_extra <- setNames(ind[[3]], ind[[1]])   # named vector: ID -> pop
  ext_fam <- read.table(paste0(ext_prefix, ".fam"),
                        col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                        stringsAsFactors = FALSE)
  for (sample_id in names(pop_lookup_extra)) {
    idx <- ext_fam$IID == sample_id
    if (any(idx)) ext_fam$FID[idx] <- pop_lookup_extra[[sample_id]]
  }
  write.table(ext_fam, paste0(ext_prefix, ".fam"),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

  n_ext_pops <- length(unique(ext_fam$FID))
  cat(sprintf("  Extended .fam: %d samples, %d unique population labels.\n\n",
              nrow(ext_fam), n_ext_pops))
}

# ---- Extract f2 for ALL needed populations at once ----
# One expensive extract_f2 call; multiple cheap qpAdm calls reuse the cache.
all_pops_needed <- unique(c(target_id, left_source, left_german, all_right_pops))
f2_dir <- file.path(outdir, paste0(basename(ext_prefix), "_f2cache"))

cat("--- Computing f2 statistics (all outgroup populations, single pass) ---\n")
cat(sprintf("  Populations: %s\n", paste(all_pops_needed, collapse = ", ")))
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)
extract_f2(ext_prefix, pops = all_pops_needed, outdir = f2_dir,
           overwrite = TRUE, format = "plink")
cat("\n")

# ---- Loop over right-population combinations ----
cat("--- qpAdm for each outgroup combination ---\n\n")
results <- list()

for (label in names(right_sets)) {
  right <- right_sets[[label]]
  cat(sprintf("  [%s]\n", label))

  f2 <- tryCatch(
    read_f2(f2_dir, pops = c(target_id, left_source, left_german, right)),
    error = function(e) {
      cat(sprintf("    ERROR reading f2: %s\n", conditionMessage(e))); NULL
    }
  )
  if (is.null(f2)) {
    results[[label]] <- data.frame(
      right_set = label, model_a_p = NA_real_, germany_weight = NA_real_,
      germany_se = NA_real_, germany_z = NA_real_, verdict = "ERROR",
      stringsAsFactors = FALSE
    )
    next
  }

  # Model A: Wielbark alone
  ma <- tryCatch(
    qpadm(f2, target = target_id, left = left_source, right = right),
    error = function(e) { cat(sprintf("    Model A error: %s\n", conditionMessage(e))); NULL }
  )
  a_p <- if (!is.null(ma)) {
    ma$rankdrop$p[ma$rankdrop$f4rank == 0]
  } else NA_real_

  # Model B: Wielbark + Germany
  mb <- tryCatch(
    qpadm(f2, target = target_id, left = c(left_source, left_german), right = right),
    error = function(e) { cat(sprintf("    Model B error: %s\n", conditionMessage(e))); NULL }
  )
  if (!is.null(mb)) {
    w  <- mb$weights$weight[mb$weights$left == left_german]
    se <- mb$weights$se[mb$weights$left == left_german]
    z  <- mb$weights$z[mb$weights$left == left_german]
  } else {
    w <- se <- z <- NA_real_
  }

  verdict <- if (is.na(z)) "ERROR" else if (abs(z) < 2) "ARTEFACT" else "SIGNIFICANT"

  cat(sprintf("    Model A p=%.4f | Germany: weight=%.3f se=%.3f z=%.2f => %s\n\n",
              ifelse(is.na(a_p), -1, a_p),
              ifelse(is.na(w),   0,   w),
              ifelse(is.na(se),  0,   se),
              ifelse(is.na(z),   0,   z),
              verdict))

  results[[label]] <- data.frame(
    right_set      = label,
    model_a_p      = round(a_p, 4),
    germany_weight = round(w,   3),
    germany_se     = round(se,  3),
    germany_z      = round(z,   3),
    verdict        = verdict,
    stringsAsFactors = FALSE
  )
}

# ---- Summary ----
res_df <- do.call(rbind, results)
rownames(res_df) <- NULL

cat("================================================================================\n")
cat("  SUMMARY TABLE\n")
cat("================================================================================\n\n")
print(res_df, row.names = FALSE)
cat("\n")

tsv_out <- file.path(outdir, "robustness_results.tsv")
write.table(res_df, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n\n", tsv_out))

# ---- Interpretation ----
n_ok    <- sum(res_df$verdict == "ARTEFACT", na.rm = TRUE)
n_total <- sum(res_df$verdict %in% c("ARTEFACT", "SIGNIFICANT"), na.rm = TRUE)

cat("================================================================================\n")
cat("  INTERPRETATION\n")
cat("================================================================================\n\n")
cat(sprintf("  Germany/CWC ARTEFACT (|z|<2): %d / %d outgroup combinations\n\n",
            n_ok, n_total))

if (n_total == 0L) {
  cat("  No results to interpret (all combinations returned ERROR).\n")
} else if (n_ok == n_total) {
  cat("  CONCLUSION STABLE: Germany/CWC proxy is a geometric artefact\n")
  cat("  in ALL tested outgroup combinations.\n")
  cat("  Conclusion is robust to outgroup choice -- result confirmed.\n")
} else if (n_ok >= ceiling(n_total * 0.75)) {
  cat("  CONCLUSION PARTIALLY STABLE: artefact in most combinations.\n")
  cat("  Check manually combinations with verdict SIGNIFICANT -- may indicate\n")
  cat("  gene flow between one of the outgroups and the tested populations.\n")
} else {
  cat("  CONCLUSION UNSTABLE: Germany/CWC gives inconsistent results depending\n")
  cat("  on outgroups. Requires manual analysis (see table above).\n")
}

cat("\n================================================================================\n")
cat("  END OF ROBUSTNESS CHECK\n")
cat("================================================================================\n")
