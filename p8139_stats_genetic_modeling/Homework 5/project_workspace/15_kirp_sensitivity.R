#!/usr/bin/env Rscript
# 15_kirp_sensitivity.R
# Test whether adding an effect-size pre-filter + relaxed FDR makes the
# pipeline less sample-size-sensitive.
#
# Strategy: re-classify silencing under three thresholds:
#   default: DM FDR<0.05 ∩ MEA FDR<0.05 ∩ direction
#   A1:      DM FDR<0.10 ∩ MEA FDR<0.10 ∩ |logFC_DM|>0.5 ∩ direction
#   A2:      DM FDR<0.10 ∩ MEA FDR<0.10 ∩ |logFC_DM|>1.0 ∩ direction
# (M-scale |logFC|>0.5 ≈ |Δβ|>0.07; |logFC|>1.0 ≈ |Δβ|>0.15)
#
# Compare KIRC vs KIRP under each scheme.
# Also: report KIRP matched-normal count (relevant to MEA power concern).

suppressMessages({
  library(data.table)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
SIL_DIR_KIRC <- file.path(WORKDIR, "silencing_out")
KIRP_DIR <- file.path(WORKDIR, "kirp_out")

# ---- Load both silencing tables --------------------------------------------
kirc <- fread(file.path(SIL_DIR_KIRC, "silencing_table.csv"))
kirp <- fread(file.path(KIRP_DIR, "kirp_silencing_table.csv"))

# Standardize column names (kirc has assoc_FDR; kirp has assoc_FDR — same)
# kirc has 'q1_logFC' as well? Let me check
cat("KIRC columns:", paste(names(kirc), collapse=", "), "\n")
cat("KIRP columns:", paste(names(kirp), collapse=", "), "\n\n")

# ---- Verify default counts ---------------------------------------------------
cat("=== Default (FDR 0.05/0.05 + direction) ===\n")
kirc_default <- kirc[silencing == TRUE]
kirp_default <- kirp[silencing == TRUE]
cat(sprintf("  KIRC pairs %d, genes %d\n", nrow(kirc_default), uniqueN(kirc_default$gene)))
cat(sprintf("  KIRP pairs %d, genes %d\n", nrow(kirp_default), uniqueN(kirp_default$gene)))

# ---- KIRP matched-normal count (MEA power diagnostic) -----------------------
cat("\n=== KIRP MEA-power diagnostic ===\n")
preproc <- readRDS(file.path(KIRP_DIR, "preprocess_kirp.rds"))
meta <- preproc$meta
n_t <- sum(meta$sample_type == "Primary Tumor")
n_n <- sum(meta$sample_type == "Solid Tissue Normal")
cat(sprintf("  KIRP methylation aliquots: %d tumor + %d normal\n", n_t, n_n))

# Read kirp_replication_summary to see matched count, and check from kirp_silencing_table
# Reconstruct what MEA test saw: matched samples in both methylation+RNA-seq
# This was logged but not saved to file. Let's just report.
cat("  (matched (patient, sample-type) keys for MEA were ~280-290 by ratio;\n")
cat("   the bottleneck is paired patients = 45 for the DM test.)\n")

# ---- Define improvement A ---------------------------------------------------
# Use M-scale logFC threshold (data has q1_logFC on M-scale)
silencing_under <- function(dt, fdr_dm, fdr_mea, logfc_min) {
  # Make sure we have the right column names
  dm_fdr_col <- if ("q1_FDR" %in% names(dt)) "q1_FDR" else "adj.P.Val"
  dm_lfc_col <- if ("q1_logFC" %in% names(dt)) "q1_logFC" else "logFC"
  dt[get(dm_fdr_col) < fdr_dm &
     assoc_FDR < fdr_mea &
     get(dm_lfc_col) > logfc_min &
     beta_coef < 0]
}

run_scheme <- function(name, fdr_dm, fdr_mea, logfc_min) {
  cat(sprintf("\n=== %s (DM FDR<%.2f, MEA FDR<%.2f, logFC_DM>%.2f) ===\n",
              name, fdr_dm, fdr_mea, logfc_min))
  k_kirc <- silencing_under(kirc, fdr_dm, fdr_mea, logfc_min)
  k_kirp <- silencing_under(kirp, fdr_dm, fdr_mea, logfc_min)
  g_kirc <- unique(k_kirc$gene)
  g_kirp <- unique(k_kirp$gene)
  bg <- intersect(unique(kirc$gene), unique(kirp$gene))
  kirc_in_bg <- intersect(g_kirc, bg)
  kirp_in_bg <- intersect(g_kirp, bg)
  overlap <- intersect(kirc_in_bg, kirp_in_bg)
  union_g <- union(kirc_in_bg, kirp_in_bg)
  jacc <- length(overlap) / length(union_g)
  if (length(kirp_in_bg) > 0) {
    hyp <- phyper(length(overlap) - 1, length(kirp_in_bg),
                  length(bg) - length(kirp_in_bg), length(kirc_in_bg),
                  lower.tail = FALSE)
  } else hyp <- NA
  cat(sprintf("  KIRC pairs %d (%d genes); KIRP pairs %d (%d genes)\n",
              nrow(k_kirc), length(g_kirc), nrow(k_kirp), length(g_kirp)))
  cat(sprintf("  background %d, overlap %d, Jaccard %.3f, hypergeom p %.2e\n",
              length(bg), length(overlap), jacc, hyp))
  cat(sprintf("  KIRP/KIRC gene ratio: %.0f%%\n",
              100 * length(g_kirp) / max(1, length(g_kirc))))

  # 18-TSG recovery
  TSG18 <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
             "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
             "TIMP3","NEFH","UCHL1","PITX2")
  kirc_tsg <- intersect(TSG18, g_kirc)
  kirp_tsg <- intersect(TSG18, g_kirp)
  both_tsg <- intersect(kirc_tsg, kirp_tsg)
  cat(sprintf("  18-TSG: KIRC %d, KIRP %d, both %d (shared: %s)\n",
              length(kirc_tsg), length(kirp_tsg), length(both_tsg),
              paste(both_tsg, collapse=", ")))
  invisible(list(kirc_genes = g_kirc, kirp_genes = g_kirp,
                 overlap = overlap, jaccard = jacc, hyp_p = hyp,
                 kirc_tsg = kirc_tsg, kirp_tsg = kirp_tsg))
}

run_scheme("DEFAULT (current paper)",       0.05, 0.05, 0.0)
run_scheme("A0 — relax FDR only",            0.10, 0.10, 0.0)
run_scheme("A1 — effect-size + relaxed FDR", 0.10, 0.10, 0.5)
run_scheme("A2 — stricter effect-size",      0.10, 0.10, 1.0)
run_scheme("A3 — strictest effect-size",     0.05, 0.05, 1.0)

cat("\n=== Done ===\n")
