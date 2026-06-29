#!/usr/bin/env Rscript
# ==============================================================================
# 02_build_aadr_subset.R -- Germany_IA artefact test, step 2/3
#
# Extracts only the populations needed for the test from the (huge) AADR HO
# panel, and converts that small subset to PLINK. Must NOT try to convert the
# entire AADR panel at once -- that triggers an OOM kill (verified: 27594
# samples x 584131 SNPs needs ~129 GB RAM for packedancestrymap_to_plink()).
#
# PITFALL FIXES (all discovered the hard way in the original session):
#   1. admixtools::packedancestrymap_to_plink() loads the FULL genotype matrix
#      into memory regardless of how many populations you eventually want --
#      it has no pops= filtering. For the full 27594-sample AADR panel this
#      needs ~129 GB RAM and gets OOM-killed on typical workstations.
#      FIX: use read_packedancestrymap(pops=...) instead, which DOES filter
#      before loading (verified: 341 samples loads in ~1.6 GB, no OOM).
#   2. admixtools has NO exported write_plink() function, despite
#      extract_samples()'s documentation mentioning one. We write PED/MAP
#      manually from the genotype matrix returned by read_packedancestrymap().
#   3. read_packedancestrymap()'s $ind tibble has columns named X1/X2/X3
#      (ID/sex/population), NOT ID/Sex/Pop -- access positionally (ind[[1]],
#      ind[[3]]), not by name.
#   4. Some SNPs come back 100% missing in a small population subset (zero
#      individuals have any genotype call) -- read_packedancestrymap() then
#      writes A1=A2="0" for that SNP, which plink1.9/2 reject downstream as
#      "identical A1/A2 alleles". FIX: detect and exclude these before any
#      merge step (handled in step 3, but flagged here at write time).
#
# Usage:
#   Rscript 02_build_aadr_subset.R <aadr_prefix> <run_outdir> [pop1,pop2,...]
#
#   <aadr_prefix>  e.g. /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB
#                  (the large, shared AADR data -- read-only, used across runs)
#   <run_outdir>   the SAME directory step 1 created (e.g. ~/MyHeritage_PW),
#                  NOT the AADR directory. aadr_subset.* will be written here,
#                  alongside pw_plink.* from step 1, so step 3 finds both in
#                  one place.
#   pop list       comma-separated, optional. Default below covers the
#                  Germany_IA artefact test (Mbuti/Yoruba/Han outgroups,
#                  Poland_IA_Wielbark local source, Germany_Esperstedt_CordedWare
#                  proxy, Poland_EarlyMedieval_Slav aggregate target).
#
# Output:
#   <run_outdir>/aadr_subset.{map,ped,bed,bim,fam}
#   <run_outdir>/aadr_subset_pop_lookup.tsv
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 02_build_aadr_subset.R <aadr_prefix> <run_outdir> [pop1,pop2,...]")
}

aadr_prefix <- args[1]
outdir      <- args[2]   # run-specific dir from step 1, e.g. ~/MyHeritage_PW
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
default_pops <- c("Mbuti", "Yoruba", "Han",
                   "Poland_EarlyMedieval_Slav",
                   "Poland_IA_Wielbark",
                   "Germany_Esperstedt_CordedWare")
target_pops <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else default_pops

cat("================================================================================\n")
cat("  STEP 2: Extract AADR population subset -> PLINK\n")
cat("================================================================================\n\n")
cat("AADR prefix : ", aadr_prefix, "\n")
cat("Output dir  : ", outdir, "\n")
cat("Populations : ", paste(target_pops, collapse = ", "), "\n\n")

# ---- Sanity check: verify population labels exist before the expensive read ----
ind_check <- read.table(paste0(aadr_prefix, ".ind"), col.names = c("ID", "Sex", "Pop"),
                         stringsAsFactors = FALSE)
n_table <- table(ind_check$Pop[ind_check$Pop %in% target_pops])
missing_pops <- setdiff(target_pops, names(n_table))
if (length(missing_pops) > 0) {
  stop("Population label(s) not found in .ind file: ", paste(missing_pops, collapse = ", "),
       "\n  -> Check exact spelling with: grep('keyword', ind$Pop, value=TRUE)",
       "\n  -> AADR population names change between versions -- do not assume",
       "     G25/Davidski cluster names match AADR labels directly.")
}
cat("Sample sizes found:\n")
print(n_table)
cat("\n")
low_n <- names(n_table)[n_table < 5]
if (length(low_n) > 0) {
  cat("  WARNING: n < 5 for:", paste(low_n, collapse = ", "), "\n")
  cat("  Results involving these populations should be treated as provisional.\n\n")
}

# ---- PITFALL FIX 1: read_packedancestrymap(pops=...), NOT packedancestrymap_to_plink() ----
# This filters BEFORE loading into memory, unlike packedancestrymap_to_plink()
# which always loads all 27594 AADR samples regardless of what you need.
cat("Reading filtered genotype data (this is memory-safe; full-panel conversion is NOT)...\n")
geno_data <- read_packedancestrymap(pref = aadr_prefix, pops = target_pops)

geno <- geno_data$geno
ind  <- geno_data$ind   # PITFALL FIX 3: columns are X1/X2/X3, not ID/Sex/Pop
snp  <- geno_data$snp

n_snp <- nrow(geno)
n_ind <- ncol(geno)
cat("Loaded", n_snp, "SNPs x", n_ind, "individuals.\n\n")

# ---- PITFALL FIX 4: detect and report 100%-missing SNPs before writing ----
# These would otherwise produce A1=A2="0" in the .bim, which plink rejects
# downstream as "identical A1 and A2 alleles" (causes --make-bed to fail, or
# in plink1.9's case, a worse failure mode: a crash partway through a merge).
all_missing <- rowSums(!is.na(geno)) == 0
n_all_missing <- sum(all_missing)
if (n_all_missing > 0) {
  cat("NOTE:", n_all_missing, "of", n_snp,
      "SNPs are 100% missing across this population subset.\n")
  cat("      These will be written with A1=A2='0' and MUST be excluded before\n")
  cat("      any merge (this is expected, not a bug -- handled in step 3).\n\n")
}

# ---- PITFALL FIX 2: manual PED/MAP write (no write_plink() in admixtools) ----
outpref <- file.path(outdir, "aadr_subset")

cat("Writing", paste0(outpref, ".map"), "...\n")
write.table(
  data.frame(snp$CHR, snp$SNP, snp$cm, snp$POS),
  file = paste0(outpref, ".map"), sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)

# Defensive: non-single-character alleles (indels) -> PLINK missing code "0"
a1 <- snp$A1; a1[nchar(a1) != 1] <- "0"
a2 <- snp$A2; a2[nchar(a2) != 1] <- "0"

cat("Building allele matrix for", n_snp, "SNPs (this is the slow step; progress every 5%)...\n")
ped_alleles <- matrix("0", nrow = n_ind, ncol = 2 * n_snp)
pb_step <- max(1, n_snp %/% 20)

for (i in seq_len(n_snp)) {
  g <- geno[i, ]
  out1 <- character(length(g)); out2 <- character(length(g))
  is2 <- !is.na(g) & g == 2; is1 <- !is.na(g) & g == 1; is0 <- !is.na(g) & g == 0
  isna <- is.na(g)
  out1[is2] <- a1[i]; out2[is2] <- a1[i]
  out1[is1] <- a1[i]; out2[is1] <- a2[i]
  out1[is0] <- a2[i]; out2[is0] <- a2[i]
  out1[isna] <- "0";  out2[isna] <- "0"
  ped_alleles[, 2 * i - 1] <- out1
  ped_alleles[, 2 * i]     <- out2
  if (i %% pb_step == 0) cat("  ", round(100 * i / n_snp), "%\n", sep = "")
}

iid <- ind[[1]]   # positional access -- column is named X1, not ID
ped_meta <- data.frame(FID = iid, IID = iid, PID = 0, MID = 0, SEX = 0, PHEN = -9)
ped_out <- cbind(ped_meta, as.data.frame(ped_alleles, stringsAsFactors = FALSE))

cat("Writing", paste0(outpref, ".ped"), "...\n")
write.table(ped_out, file = paste0(outpref, ".ped"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# ---- Write a population-label lookup table for step 3 to fix the .fam FID ----
# read_packedancestrymap() gives us per-individual population in ind[[3]];
# .ped/.fam FID defaults to the individual ID (no population), so step 3
# needs this lookup to restore population labels as FID, which qpAdm/extract_f2
# require for pops= matching.
pop_lookup <- data.frame(ID = ind[[1]], Pop = ind[[3]])
write.table(pop_lookup, file = file.path(outdir, "aadr_subset_pop_lookup.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote population lookup table: aadr_subset_pop_lookup.tsv\n")

if (n_all_missing > 0) {
  bad_snp_ids <- snp$SNP[all_missing]
  writeLines(bad_snp_ids, file.path(outdir, "aadr_subset_allmissing_snps.txt"))
  cat("Wrote list of", n_all_missing, "all-missing SNP IDs:",
      "aadr_subset_allmissing_snps.txt\n")
}

cat("\n================================================================================\n")
cat("  STEP 2 DONE.\n")
cat("  Convert to binary PLINK with:\n")
cat("    plink2 --pedmap", outpref, "--make-bed --out", outpref, "\n")
cat("  Then exclude all-missing SNPs (if any) and run 03_merge_and_test.sh\n")
cat("================================================================================\n")
