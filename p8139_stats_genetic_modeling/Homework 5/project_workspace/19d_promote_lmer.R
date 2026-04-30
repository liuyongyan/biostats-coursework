#!/usr/bin/env Rscript
# 19d_promote_lmer.R
# Promote the lmer MEA fit (sensitivity_lmer/) to primary by overwriting
# silencing_out/silencing_table.csv. Preserve original pooled fit as
# silencing_table_pooled.csv for provenance.
#
# Schema preserved exactly:
#   cpg, gene, promoter_group,
#   partial_r,    beta_coef,  assoc_p,
#   pct_responder, assoc_FDR,
#   q1_logFC,     q1_FDR,
#   silencing,    activation, pattern

suppressPackageStartupMessages({ library(data.table) })

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
SIL_DIR  <- file.path(WORKDIR, "silencing_out")
LMER_DIR <- file.path(WORKDIR, "sensitivity_lmer")

pooled <- fread(file.path(SIL_DIR, "silencing_table.csv"))
lmer_t <- fread(file.path(LMER_DIR, "silencing_table_lmer.csv.gz"))

cat("Pooled rows:", nrow(pooled), "\n")
cat("lmer rows:  ", nrow(lmer_t), "\n")
stopifnot(nrow(pooled) == nrow(lmer_t))

# Backup pooled
file.copy(file.path(SIL_DIR, "silencing_table.csv"),
          file.path(SIL_DIR, "silencing_table_pooled.csv"),
          overwrite = TRUE)
cat("Backed up pooled -> silencing_table_pooled.csv\n")

# pct_responder, partial_r are properties of the data (not the MEA fit) — carry over
m <- merge(
  lmer_t[, .(cpg, gene, promoter_group,
             beta_coef = beta_coef_lmer, se_lmer,
             assoc_p   = assoc_p_lmer,
             assoc_FDR = assoc_FDR_lmer,
             q1_logFC, q1_FDR,
             silencing = silencing_lmer,
             activation = activation_lmer)],
  pooled[, .(cpg, gene, partial_r, pct_responder)],
  by = c("cpg", "gene"),
  all.x = TRUE
)
stopifnot(nrow(m) == nrow(lmer_t))

# Recompute pattern from new silencing flag + same pct_responder
m[, pattern := NA_character_]
m[silencing == TRUE & pct_responder >= 0.3, pattern := "cohort_wide"]
m[silencing == TRUE & pct_responder <  0.3, pattern := "subgroup"]

# Final schema (column order matches original)
out <- m[, .(cpg, gene, promoter_group,
             partial_r, beta_coef, assoc_p,
             pct_responder, assoc_FDR,
             q1_logFC, q1_FDR,
             silencing, activation, pattern,
             se_lmer)]   # keep se_lmer as a non-canonical extra column (downstream scripts ignore unknown cols)

setorder(out, assoc_p)
fwrite(out, file.path(SIL_DIR, "silencing_table.csv"))
cat("Wrote new silencing_table.csv (lmer-primary)\n")

# Summary
n_sil <- sum(out$silencing)
n_act <- sum(out$activation)
n_cw  <- sum(out$pattern == "cohort_wide", na.rm = TRUE)
n_sg  <- sum(out$pattern == "subgroup",    na.rm = TRUE)
n_sil_genes <- uniqueN(out$gene[out$silencing])
n_sil_cpgs  <- uniqueN(out$cpg[out$silencing])

cat(sprintf("\n=== New silencing table summary (lmer-primary) ===\n"))
cat(sprintf("  Total pairs:          %d\n", nrow(out)))
cat(sprintf("  Silencing pairs:      %d\n", n_sil))
cat(sprintf("    cohort_wide:        %d (%.1f%%)\n", n_cw, 100*n_cw/n_sil))
cat(sprintf("    subgroup:           %d (%.1f%%)\n", n_sg, 100*n_sg/n_sil))
cat(sprintf("  Silencing CpGs:       %d\n", n_sil_cpgs))
cat(sprintf("  Silencing genes:      %d\n", n_sil_genes))
cat(sprintf("  Activation pairs:     %d\n", n_act))

# Pre/post comparison
old_n_sil <- sum(pooled$silencing)
old_n_act <- sum(pooled$activation)
old_genes <- uniqueN(pooled$gene[pooled$silencing])
cat(sprintf("\n=== Diff vs pooled-primary ===\n"))
cat(sprintf("  silencing pairs:  %d -> %d  (%+d)\n", old_n_sil, n_sil,  n_sil - old_n_sil))
cat(sprintf("  silencing genes:  %d -> %d  (%+d)\n", old_genes, n_sil_genes, n_sil_genes - old_genes))
cat(sprintf("  activation pairs: %d -> %d  (%+d)\n", old_n_act, n_act, n_act - old_n_act))

# VHL
vhl <- out[cpg == "cg13672843" & gene == "VHL"]
cat(sprintf("\n=== VHL/cg13672843 in new table ===\n"))
print(vhl[, .(cpg, beta_coef, assoc_p, assoc_FDR, q1_logFC, q1_FDR,
              silencing, pattern, pct_responder)])
