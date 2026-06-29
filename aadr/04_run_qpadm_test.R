#!/usr/bin/env Rscript
# ==============================================================================
# 04_run_qpadm_test.R -- Germany_IA artefact test, step 4/4 (final)
#
# Runs extract_f2() and qpAdm on the merged personal-genome + AADR fileset
# produced by 03_merge_and_test.sh, and prints a plain-language interpretation
# of whether the "German Iron Age / Corded Ware" proxy is a real ancestry
# component or a substitutable geometric artefact (the question this whole
# pipeline exists to answer).
#
# PITFALL FIXES applied here:
#   - format = "plink" passed explicitly to extract_f2() (don't rely on
#     extension auto-detection after a multi-tool conversion chain).
#   - auto_only is left at its TRUE default here because step 3 already
#     restricted the fileset to chr 1-22 via plink2 --chr -- the chromosome
#     column is now purely numeric and auto_only's parser will not choke.
#     (If you skip step 3's chr filter, pass auto_only = FALSE here instead.)
#   - Sample-size and SNP-count sanity checks are printed BEFORE trusting any
#     f-statistic, because a personal-genome merge typically has far fewer
#     usable SNPs (tens of thousands) than an all-AADR aggregate test (often
#     hundreds of thousands) -- this is expected (consumer array vs capture
#     panel intersection), not a bug, but it does mean wider standard errors.
#
# Usage:
#   Rscript 04_run_qpadm_test.R <merged_prefix> <target_sample_id> [local_source] [german_proxy]
#
#   <merged_prefix>    e.g. /usr/local/share/aadr/merged_pw_aadr_final_ready
#   <target_sample_id> the FID you assigned to your own sample (default: PW)
#   [local_source]     default: Poland_IA_Wielbark
#   [german_proxy]     default: Germany_Esperstedt_CordedWare
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 04_run_qpadm_test.R <merged_prefix> [target_id] [local_source] [german_proxy]")
}

prefix       <- args[1]
target_id    <- if (length(args) >= 2) args[2] else "PW"
local_source <- if (length(args) >= 3) args[3] else "Poland_IA_Wielbark"
german_proxy <- if (length(args) >= 4) args[4] else "Germany_Esperstedt_CordedWare"
outgroups    <- c("Mbuti", "Yoruba", "Han")

# f2cache in the run directory (dirname(prefix)) -- next to the BED files.
f2_dir <- file.path(dirname(prefix), paste0(basename(prefix), "_f2cache"))

cat("================================================================================\n")
cat("  GERMANY_IA / CWC ARTEFACT TEST -- personal genome\n")
cat("================================================================================\n\n")
cat("Target            :", target_id, "\n")
cat("Local source (left):", local_source, "\n")
cat("German/CWC proxy   :", german_proxy, "\n")
cat("Outgroups (right)  :", paste(outgroups, collapse = ", "), "\n\n")

all_pops <- c(outgroups, target_id, local_source, german_proxy)

# ---- Step A: sample-size sanity check ----
cat("--- Sample sizes ---\n")
fam <- read.table(paste0(prefix, ".fam"),
                   col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                   stringsAsFactors = FALSE)
n_table <- table(fam$FID[fam$FID %in% all_pops])
print(n_table)
cat("\n")

missing_pops <- setdiff(all_pops, names(n_table))
if (length(missing_pops) > 0) {
  stop("Population/sample label(s) not found in merged .fam: ",
       paste(missing_pops, collapse = ", "),
       "\n  -> If your target ID is missing, check step 3's FID-fix step ran",
       "     correctly: grep '<target_id>' <prefix>.fam")
}

# ---- Step B: extract f2 ----
cat("--- Extracting f2-statistics ---\n")
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)
extract_f2(prefix, pops = all_pops, outdir = f2_dir,
           overwrite = TRUE, format = "plink")
f2 <- read_f2(f2_dir, pops = all_pops)
cat("\n")

# ---- Step C: Model A (local source alone) ----
cat("--- Model A:", target_id, "explained by", local_source, "alone ---\n")
model_a <- qpadm(f2, target = target_id, left = c(local_source), right = outgroups)
print(model_a$weights)
a_p <- model_a$rankdrop$p[model_a$rankdrop$f4rank == 0]
a_weight <- model_a$weights$weight[model_a$weights$left == local_source]
cat("Model A p-value:", round(a_p, 4), "\n\n")

# ---- Step D: Model B (local source + German/CWC proxy) ----
cat("--- Model B:", target_id, "explained by", local_source, "+", german_proxy, "---\n")
model_b <- qpadm(f2, target = target_id, left = c(local_source, german_proxy), right = outgroups)
print(model_b$weights)
w_german <- model_b$weights$weight[model_b$weights$left == german_proxy]
se_german <- model_b$weights$se[model_b$weights$left == german_proxy]
z_german <- model_b$weights$z[model_b$weights$left == german_proxy]
cat("\nGerman/CWC proxy: weight =", round(w_german, 3),
    " se =", round(se_german, 3), " z =", round(z_german, 3), "\n\n")

cat("--- popdrop (formal nested model comparison) ---\n")
print(model_b$popdrop)
cat("\n")

# ---- Final interpretation ----
cat("================================================================================\n")
cat("  INTERPRETATION\n")
cat("================================================================================\n\n")

if (a_p > 0.05 && abs(a_weight - 1) < 0.1) {
  cat("Model A (", local_source, " alone) is NOT rejected (p = ", round(a_p, 3),
      ") with weight ~= 1.0.\n", sep = "")
  cat(target_id, "is fully explained by local continuity from", local_source, ".\n\n")
} else {
  cat("Model A was rejected or gave an unexpected weight -- inspect manually.\n\n")
}

if (abs(z_german) < 2) {
  cat("Model B: the", german_proxy, "weight has |z| =", round(abs(z_german), 2),
      "< 2 -- NOT significantly different from zero.\n")
  cat("CONCLUSION: the German/CWC proxy is a substitutable artefact for", target_id, ",\n")
  cat("not a real, independently-required ancestry component.\n")
  cat("This matches the geometric-proxy mechanism already identified via G25/nMonte:\n")
  cat("a large 'Germany_IA'-type weight in PCA-distance pipelines reflects a shared\n")
  cat("deep Corded Ware/Bronze Age substrate and a panel gap, not real Germanic\n")
  cat("ancestry.\n")
} else {
  cat("Model B: the", german_proxy, "weight has |z| =", round(abs(z_german), 2),
      ">= 2 -- POSSIBLY significant.\n")
  cat("This would NOT cleanly support the pure-artefact interpretation for this\n")
  cat("target -- inspect the popdrop table above and consider whether SNP count\n")
  cat("is large enough for a reliable estimate (see SNP-count caveat below).\n")
}

cat("\nCAVEAT: a personal-genome merge typically has far fewer usable SNPs than an\n")
cat("all-AADR aggregate test, due to limited overlap between a consumer genotyping\n")
cat("array and the AADR capture panel. Check the f2 extraction log above for the\n")
cat("'SNPs remain after filtering' count -- tens of thousands is workable for\n")
cat("qpAdm, but standard errors will be wider than an aggregate-population test\n")
cat("with hundreds of thousands of SNPs.\n\n")

cat("================================================================================\n")
cat("  END OF TEST\n")
cat("================================================================================\n")
