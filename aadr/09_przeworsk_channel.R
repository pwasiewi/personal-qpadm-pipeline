#!/usr/bin/env Rscript
# ==============================================================================
# 09_przeworsk_channel.R -- Przeworsk vs Wielbark "channel of descent" test
#
# HYPOTHESIS (user, 2026-06-29):
#   PW has a medieval Greater-Poland (Wielkopolska) genealogical background.
#   The line of descent from the Iron Age Lusatian substrate may have run
#   through the PRZEWORSK culture (central/southern Poland, autochthonous,
#   Slavic-substrate channel) rather than the WIELBARK culture (Pomeranian,
#   Gothic/Scandinavian-admixed channel). Geographically Przeworsk covered
#   Wielkopolska/Silesia/Lesser Poland/Mazovia; Wielbark was Pomerania-centred.
#
# WHY THIS WAS NEVER TESTED BEFORE:
#   Step 7 hard-skips any source with n<3 (line ~126). AADR v66 HO has only
#   Poland_IA_Przeworsk n=2 (PCA0011.SG, PCA0012.SG -- both female, both from
#   Gaski, Kuyavia, 100-300 CE, high-coverage shotgun: 160K/100K HO SNPs hit).
#   So Przeworsk fell out of every prior run automatically. Here n>=2 is allowed
#   with an explicit caveat.
#
# METHODOLOGICAL NOTE -- read before interpreting:
#   Wielbark, EarlySlav, Lusatian are all COLLINEAR at ~50K SNP (closed finding).
#   Przeworsk (also a Lusatian derivative) will almost certainly be collinear too
#   -> single-source qpAdm will "pass" trivially (n=2 inflates SE and p further).
#   A pass/fail therefore CANNOT decide which channel. The decisive instrument is
#   the f4 symmetry statistic:
#     (B1) f4(Mbuti, Sweden_IA; Przeworsk, Wielbark)
#          -> do the channels even differ? Expect Wielbark closer to Sweden_IA
#             (Gothic/Scandinavian admixture). If |Z|<3 the two cultures are
#             genetically indistinguishable here and the question is moot at 50K.
#     (B2) f4(Mbuti, PW; Wielbark, Przeworsk)        <-- THE decisive test
#          -> Z ~ 0      : PW symmetric to both, channel UNRESOLVED.
#             Z > +3     : PW shares more drift with Przeworsk  (Przeworsk channel)
#             Z < -3     : PW shares more drift with Wielbark   (Wielbark channel)
#     (B3) f4(Mbuti, PW; Sweden_IA, Przeworsk/Wielbark) -- Scandinavian affinity ctx
#
# Usage (run from ~/Claude or any writable dir):
#   Rscript aadr/09_przeworsk_channel.R \
#     ./PW_FF/merged_pw_aadr_final_ready \
#     /usr/local/share/aadr/v66.p1_HO.aadr.patch.PUB \
#     PW
#
# Output:
#   <outdir>/przeworsk_channel_results_<TARGET>.tsv   (single-source + 2-source)
#   <outdir>/przeworsk_channel_f4_<TARGET>.tsv        (f4 symmetry tests)
#   <outdir>/merged_pw_aadr_chan9.{bed,bim,fam}       (idempotent)
# ==============================================================================

suppressPackageStartupMessages(library(admixtools))

options(width = max(200L, as.integer(Sys.getenv("COLUMNS", "0")), na.rm = TRUE))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 09_przeworsk_channel.R <merged_prefix> <aadr_prefix> [target]")

merged_prefix <- args[1]
aadr_prefix   <- args[2]
target_id     <- if (length(args) >= 3) args[3] else "PW"
outdir        <- dirname(merged_prefix)

cat("================================================================================\n")
cat("  STEP 9: Przeworsk vs Wielbark channel-of-descent test\n")
cat("================================================================================\n\n")
cat(sprintf("  Merged prefix : %s\n", merged_prefix))
cat(sprintf("  AADR prefix   : %s\n", aadr_prefix))
cat(sprintf("  Target        : %s\n\n", target_id))

right_pops <- c("Mbuti", "Yoruba", "Han")

# Channel candidates + substrate context + Scandinavian-shift probe.
# n>=2 allowed here (Przeworsk is the whole point); loud caveat for small n.
chan_pops <- c(
  "Poland_IA_Przeworsk",   # n=2  -- THE hypothesis channel (Gaski, Kuyavia, 100-300 CE)
  "Poland_IA_Wielbark",    # n=59 -- the rival channel (Pomerania/Mazovia, I-III CE)
  "Ukraine_IA_Lusatian",   # n=3  -- shared IA substrate (best 1-source so far)
  "Poland_EarlySlav",      # n=25 -- early medieval Slavic endpoint
  "Sweden_IA"              # n=8  -- Scandinavian IA, Gothic-admixture probe for Wielbark
)

MIN_N <- 2L   # deliberately below step 7's 3, to admit Przeworsk

# ---- Verify populations + sample sizes ----
cat("--- Verifying populations in AADR .ind ---\n")
ind_all <- read.table(paste0(aadr_prefix, ".ind"),
                      col.names = c("ID", "Sex", "Pop"),
                      stringsAsFactors = FALSE)
n_by_pop  <- table(ind_all$Pop)
confirmed <- c()
for (p in chan_pops) {
  n <- if (p %in% names(n_by_pop)) n_by_pop[[p]] else 0L
  flag <- if (n > 0L && n < 3L) "  <-- SMALL n, wide SE: treat result as indicative only" else ""
  if (n < MIN_N) {
    cat(sprintf("  SKIP %-28s  n=%d (< %d)\n", p, n, MIN_N))
  } else {
    cat(sprintf("  OK   %-28s  n=%d%s\n", p, n, flag))
    confirmed <- c(confirmed, p)
  }
}
cat("\n")
if (!("Poland_IA_Przeworsk" %in% confirmed))
  stop("Poland_IA_Przeworsk not available -- nothing to test.")

# ---- Which pops already in merged file ----
fam_existing <- read.table(paste0(merged_prefix, ".fam"),
                           col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                           stringsAsFactors = FALSE)
existing_fids <- unique(fam_existing$FID)
extra_needed  <- setdiff(confirmed, existing_fids)

cat(sprintf("  Already in file : %s\n",
            paste(sort(intersect(confirmed, existing_fids)), collapse = ", ")))
cat(sprintf("  Need from AADR  : %s\n\n",
            if (length(extra_needed) == 0L) "none"
            else paste(extra_needed, collapse = ", ")))

# ---- Build extended fileset (idempotent) -- pattern from step 7 ----
chan_prefix <- file.path(outdir, "merged_pw_aadr_chan9")

if (file.exists(paste0(chan_prefix, ".bed")) && length(extra_needed) > 0L) {
  chan_fids <- unique(read.table(paste0(chan_prefix, ".fam"),
                                 col.names = c("FID","IID","PID","MID","SEX","PHEN"),
                                 stringsAsFactors = FALSE)$FID)
  stale <- setdiff(extra_needed, chan_fids)
  if (length(stale) > 0L) {
    cat(sprintf("  chan9 exists but missing %d pops -- rebuilding (%s).\n\n",
                length(stale), paste(stale, collapse = ", ")))
    invisible(file.remove(paste0(chan_prefix, c(".bed", ".bim", ".fam"))))
  }
}

if (length(extra_needed) == 0L) {
  cat("  All populations already in merged file -- using it directly.\n\n")
  chan_prefix <- merged_prefix

} else if (file.exists(paste0(chan_prefix, ".bed"))) {
  cat(sprintf("  Channel fileset already complete: %s.bed -- skipping build.\n\n", chan_prefix))

} else {
  cat(sprintf("--- Reading %d populations from AADR ---\n", length(extra_needed)))
  gd   <- read_packedancestrymap(pref = aadr_prefix, pops = extra_needed, verbose = TRUE)
  geno <- gd$geno; ind <- gd$ind; snp <- gd$snp
  n_snp <- nrow(geno); n_ind <- ncol(geno)
  cat(sprintf("  Loaded %d SNPs x %d individuals.\n\n", n_snp, n_ind))

  extra_stem <- file.path(outdir, "extra_chan9_pops")
  a1 <- snp$A1; a1[nchar(a1) != 1] <- "0"
  a2 <- snp$A2; a2[nchar(a2) != 1] <- "0"

  cat(sprintf("  Building allele matrix (%d SNPs x %d ind)...\n", n_snp, n_ind))
  ped_alleles <- matrix("0", nrow = n_ind, ncol = 2L * n_snp)
  pb_step <- max(1L, n_snp %/% 10L)
  for (i in seq_len(n_snp)) {
    g <- geno[i, ]
    is2 <- !is.na(g) & g == 2L; is1 <- !is.na(g) & g == 1L; is0 <- !is.na(g) & g == 0L
    o1 <- character(length(g)); o2 <- character(length(g))
    o1[is2] <- a1[i]; o2[is2] <- a1[i]
    o1[is1] <- a1[i]; o2[is1] <- a2[i]
    o1[is0] <- a2[i]; o2[is0] <- a2[i]
    o1[is.na(g)] <- "0"; o2[is.na(g)] <- "0"
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

  if (system2("plink2", c("--pedmap", extra_stem, "--make-bed", "--out", extra_stem)) != 0L)
    stop("plink2 --pedmap failed for extra_chan9_pops")

  cat("\n--- Merging into existing file ---\n")
  if (system2("p-link", c("--bfile", merged_prefix,
                          "--bmerge", paste0(extra_stem, ".bed"),
                          paste0(extra_stem, ".bim"), paste0(extra_stem, ".fam"),
                          "--make-bed", "--out", chan_prefix)) != 0L)
    stop("p-link --bmerge failed for chan9")

  pop_lookup_extra <- setNames(ind[[3]], ind[[1]])
  chan_fam <- read.table(paste0(chan_prefix, ".fam"),
                         col.names = c("FID", "IID", "PID", "MID", "SEX", "PHEN"),
                         stringsAsFactors = FALSE)
  for (sid in names(pop_lookup_extra)) {
    idx <- chan_fam$IID == sid
    if (any(idx)) chan_fam$FID[idx] <- pop_lookup_extra[[sid]]
  }
  write.table(chan_fam, paste0(chan_prefix, ".fam"),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  chan9 .fam: %d samples, %d population labels.\n\n",
              nrow(chan_fam), length(unique(chan_fam$FID))))
}

# ---- f2 cache (target-specific, per pitfall #14) ----
all_pops_f2 <- unique(c(target_id, confirmed, right_pops))
f2_dir <- file.path(outdir, paste0(basename(chan_prefix), "_", target_id, "_f2cache"))
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)

cat("--- Computing f2 (single pass) ---\n")
cat(sprintf("  Populations: %s\n", paste(all_pops_f2, collapse = ", ")))
extract_f2(chan_prefix, pops = all_pops_f2, outdir = f2_dir,
           overwrite = TRUE, format = "plink")
f2b <- read_f2(f2_dir, pops = all_pops_f2)
cat("\n")

# ==============================================================================
# A. SINGLE-SOURCE qpAdm  (context: is each channel an adequate sole source?)
# ==============================================================================
cat("================================================================================\n")
cat("  A. Single-source qpAdm  (collinearity expected -- pass != channel proof)\n")
cat("================================================================================\n\n")

src_pops <- intersect(c("Poland_IA_Przeworsk", "Poland_IA_Wielbark",
                        "Ukraine_IA_Lusatian", "Poland_EarlySlav"), confirmed)
one_src <- list()
for (src in src_pops) {
  m <- tryCatch(qpadm(f2b, target = target_id, left = src, right = right_pops),
                error = function(e) { cat(sprintf("  ERR %s: %s\n", src, conditionMessage(e))); NULL })
  if (is.null(m)) next
  p_val <- m$rankdrop$p[m$rankdrop$f4rank == 0]
  w     <- m$weights$weight[m$weights$left == src]
  se    <- m$weights$se[m$weights$left == src]
  z     <- m$weights$z[m$weights$left == src]
  n_src <- n_by_pop[[src]]
  verdict <- if (is.na(p_val)) "ERROR" else if (p_val >= 0.05) "OK" else "REJECTED"
  cat(sprintf("  %-28s (n=%2d)  p=%.4f  weight=%.3f  z=%.2f  => %s\n",
              src, n_src, p_val, w, z, verdict))
  one_src[[src]] <- data.frame(test = "1src", source = src, n = n_src,
                               p = round(p_val, 4), weight = round(w, 3),
                               se = round(se, 3), z = round(z, 3),
                               verdict = verdict, stringsAsFactors = FALSE)
}
cat("\n")

# ==============================================================================
# B. f4 SYMMETRY TESTS  (the decisive instrument for "which channel")
# ==============================================================================
cat("================================================================================\n")
cat("  B. f4 symmetry tests\n")
cat("================================================================================\n\n")

have <- function(...) all(c(...) %in% confirmed)
f4_rows <- list()
run_f4 <- function(label, w, x, y, z, note) {
  if (!have(w, x, y, z) && !all(c(w, x, y, z) %in% c(target_id, right_pops, confirmed))) {
    cat(sprintf("  SKIP %s (missing pop)\n", label)); return(invisible())
  }
  r <- tryCatch(f4(f2b, w, x, y, z),
                error = function(e) { cat(sprintf("  ERR %s: %s\n", label, conditionMessage(e))); NULL })
  if (is.null(r)) return(invisible())
  est <- r$est[1]; se <- r$se[1]; zz <- r$z[1]; pp <- r$p[1]
  sig <- if (is.na(zz)) "?" else if (abs(zz) >= 3) "*** SIGNIFICANT" else if (abs(zz) >= 2) "* marginal" else "n.s. (symmetric)"
  cat(sprintf("  %s\n      f4(%s, %s; %s, %s) = %+.6f  SE=%.6f  Z=%+.2f  %s\n      %s\n\n",
              label, w, x, y, z, est, se, zz, sig, note))
  f4_rows[[label]] <<- data.frame(label = label,
                                  f4 = sprintf("f4(%s,%s;%s,%s)", w, x, y, z),
                                  est = round(est, 6), se = round(se, 6),
                                  z = round(zz, 3), p = signif(pp, 4),
                                  signif = sig, stringsAsFactors = FALSE)
}

# B1 -- do the two channels even differ? (Scandinavian/Gothic shift of Wielbark)
run_f4("B1 channel differentiation",
       "Mbuti", "Sweden_IA", "Poland_IA_Przeworsk", "Poland_IA_Wielbark",
       "Z<0 => Wielbark closer to Sweden_IA than Przeworsk is (Gothic admixture, channels differ).")

# B2 -- THE decisive test: which channel is the target closer to?
run_f4("B2 *** WHICH CHANNEL (decisive)",
       "Mbuti", target_id, "Poland_IA_Wielbark", "Poland_IA_Przeworsk",
       sprintf("Z>0 => %s closer to Wielbark; Z<0 => closer to Przeworsk; |Z|<3 => unresolved.", target_id))

# B3 -- Scandinavian affinity of the target via each channel (context)
run_f4("B3a Scandinavian affinity vs Przeworsk",
       "Mbuti", target_id, "Sweden_IA", "Poland_IA_Przeworsk",
       sprintf("Z<0 => %s closer to Przeworsk than to Sweden_IA (no Scandinavian pull).", target_id))
run_f4("B3b Scandinavian affinity vs Wielbark",
       "Mbuti", target_id, "Sweden_IA", "Poland_IA_Wielbark",
       "Compare with B3a: more-negative => stronger pull to that channel.")

# B4 -- robustness of B2 with a different outgroup in position 1
run_f4("B4 B2 robustness (Han outgroup)",
       "Han", target_id, "Poland_IA_Wielbark", "Poland_IA_Przeworsk",
       "Same expectation as B2; sign/|Z| should agree if B2 is real.")

# ==============================================================================
# C. TWO-SOURCE COMPETITION  (Lusatian base + each channel as 2nd source)
# ==============================================================================
cat("================================================================================\n")
cat("  C. Two-source: Ukraine_IA_Lusatian base + {Przeworsk | Wielbark}\n")
cat("     (collinearity likely -> watch for |weight|>>1 / |z|~0 = UNIDENTIFIABLE)\n")
cat("================================================================================\n\n")

two_src <- list()
base <- "Ukraine_IA_Lusatian"
if (base %in% confirmed) {
  for (second in intersect(c("Poland_IA_Przeworsk", "Poland_IA_Wielbark"), confirmed)) {
    m <- tryCatch(qpadm(f2b, target = target_id, left = c(base, second), right = right_pops),
                  error = function(e) { cat(sprintf("  ERR %s+%s: %s\n", base, second, conditionMessage(e))); NULL })
    if (is.null(m)) next
    p_val <- m$rankdrop$p[m$rankdrop$f4rank == 1]
    w2 <- m$weights$weight[m$weights$left == second]
    s2 <- m$weights$se[m$weights$left == second]
    z2 <- m$weights$z[m$weights$left == second]
    unident <- isTRUE(abs(w2) > 5 | (abs(z2) < 0.5 & abs(w2) > 1))
    verdict <- if (is.na(p_val)) "ERROR"
               else if (unident) "UNIDENTIFIABLE"
               else if (p_val < 0.05) "REJECTED"
               else if (w2 <= 0) "OK-base-sufficient"
               else "OK-2nd-contributes"
    cat(sprintf("  %s + %-22s  p=%.4f  w(2nd)=%+.3f  z=%+.2f  => %s\n",
                base, second, p_val, w2, z2, verdict))
    two_src[[second]] <- data.frame(test = "2src", source = paste0(base, "+", second),
                                    n = n_by_pop[[second]], p = round(p_val, 4),
                                    weight = round(w2, 3), se = round(s2, 3),
                                    z = round(z2, 3), verdict = verdict,
                                    stringsAsFactors = FALSE)
  }
} else {
  cat("  Base Ukraine_IA_Lusatian unavailable -- skipping 2-source.\n")
}
cat("\n")

# ==============================================================================
# Save results
# ==============================================================================
qp_df <- do.call(rbind, c(one_src, two_src))
if (!is.null(qp_df)) {
  rownames(qp_df) <- NULL
  tsv1 <- file.path(outdir, sprintf("przeworsk_channel_results_%s.tsv", target_id))
  write.table(qp_df, tsv1, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved qpAdm results: %s\n", tsv1))
}
f4_df <- do.call(rbind, f4_rows)
if (!is.null(f4_df)) {
  rownames(f4_df) <- NULL
  tsv2 <- file.path(outdir, sprintf("przeworsk_channel_f4_%s.tsv", target_id))
  write.table(f4_df, tsv2, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved f4 tests     : %s\n", tsv2))
}
cat("\n")

cat("================================================================================\n")
cat("  HOW TO READ THIS\n")
cat("================================================================================\n")
cat("  * B2 is the verdict on the user's hypothesis. If |Z|<3 the Przeworsk and\n")
cat("    Wielbark channels are statistically indistinguishable from the target's\n")
cat("    standpoint at this SNP count -- consistent with the closed collinearity\n")
cat("    finding; the Wielkopolska/Przeworsk route is then plausible but unproven.\n")
cat("  * Section A passing for Przeworsk is EXPECTED and does NOT confirm the\n")
cat("    channel (n=2, collinear). Cite B1/B2, not A.\n")
cat("  * If B1 itself is n.s., Przeworsk and Wielbark are indistinguishable from\n")
cat("    each other here -> the question cannot be answered without WGS.\n")
cat("================================================================================\n")
cat("  END STEP 9\n")
cat("================================================================================\n")
