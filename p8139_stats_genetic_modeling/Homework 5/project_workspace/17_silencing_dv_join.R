#!/usr/bin/env Rscript
# 17_silencing_dv_join.R
#
# Join differential variability (DV) results into the silencing table and
# evaluate the relationship between the responder fraction (pi_resp,
# §2.3 of the report — a non-parametric DV proxy) and the formal varFit
# DV statistic.
#
# Key questions:
#   1. What fraction of silencing pairs are DV-FDR-significant?
#   2. Spearman rho between pct_responder and DV |t| (or log var ratio)?
#   3. Among subgroup-driven silencing pairs (pct_responder < 0.3), is the
#      DV signal stronger than among cohort-wide silencing pairs?
#   4. For the 18-gene benchmark, do silencing-recovered genes also show DV?

suppressPackageStartupMessages({
  library(data.table)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
SIL_DIR <- file.path(WORKDIR, "silencing_out")
DV_DIR  <- file.path(WORKDIR, "dv_out")
OUT_DIR <- DV_DIR

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")

# ---- Load --------------------------------------------------------------------
log_msg("[1/4] Loading silencing table and DV table...")
sil <- fread(file.path(SIL_DIR, "silencing_table.csv"))
dv  <- fread(file.path(DV_DIR, "dv_table.csv"))
log_msg("  silencing pairs: %d (silencing=%d, activation=%d)",
        nrow(sil), sum(sil$silencing), sum(sil$activation))
log_msg("  DV table CpGs:   %d", nrow(dv))

# ---- Join --------------------------------------------------------------------
log_msg("[2/4] Joining DV onto silencing table by CpG...")
sil_dv <- merge(sil, dv[, .(cpg, dv_logRatio, dv_t, dv_p, dv_FDR,
                            var_tumor, var_normal)],
                by = "cpg", all.x = TRUE)
log_msg("  Pairs with DV info: %d / %d", sum(!is.na(sil_dv$dv_p)), nrow(sil_dv))

fwrite(sil_dv, file.path(OUT_DIR, "silencing_table_with_dv.csv"))

# ---- Q1: fraction of silencing pairs that are DV-significant ----------------
log_msg("[3/4] DV statistics within the silencing list...")
sil_only <- sil_dv[silencing == TRUE]
n_sil    <- nrow(sil_only)
n_sil_dv <- sum(sil_only$dv_FDR < 0.05, na.rm = TRUE)
n_sil_dv_tumor_high <- sum(sil_only$dv_FDR < 0.05 & sil_only$dv_logRatio > 0,
                            na.rm = TRUE)
log_msg("  Silencing pairs:                       %d", n_sil)
log_msg("  Silencing AND DV (FDR<0.05):           %d (%.1f%%)",
        n_sil_dv, 100 * n_sil_dv / n_sil)
log_msg("  Silencing + DV with var(tumor)>normal: %d (%.1f%%)",
        n_sil_dv_tumor_high, 100 * n_sil_dv_tumor_high / n_sil)

# Subgroup vs cohort-wide DV signal
sg <- sil_only[pattern == "subgroup"]
cw <- sil_only[pattern == "cohort_wide"]
log_msg("  Subgroup pairs:    n=%d, median |dv_t|=%.2f, median log2 var ratio=%.2f",
        nrow(sg), median(abs(sg$dv_t), na.rm=TRUE),
        median(sg$dv_logRatio, na.rm=TRUE))
log_msg("  Cohort-wide pairs: n=%d, median |dv_t|=%.2f, median log2 var ratio=%.2f",
        nrow(cw), median(abs(cw$dv_t), na.rm=TRUE),
        median(cw$dv_logRatio, na.rm=TRUE))
sg_cw_t <- wilcox.test(abs(sg$dv_t), abs(cw$dv_t))$p.value
log_msg("  Wilcoxon |dv_t| subgroup vs cohort-wide: p=%.2e", sg_cw_t)

# ---- Q2: Spearman rho(pct_responder, DV statistic) --------------------------
log_msg("[4/4] Spearman correlation between pct_responder and DV...")
ok <- !is.na(sil_only$dv_t) & !is.na(sil_only$pct_responder)
rho_t  <- cor(sil_only$pct_responder[ok], abs(sil_only$dv_t[ok]),
              method = "spearman")
rho_lr <- cor(sil_only$pct_responder[ok], sil_only$dv_logRatio[ok],
              method = "spearman")
log_msg("  rho(pct_responder, |dv_t|):       %.3f  (n=%d)", rho_t, sum(ok))
log_msg("  rho(pct_responder, dv_logRatio):  %.3f  (n=%d)", rho_lr, sum(ok))

# ---- Per-CpG overlap on the FULL silencing universe (not just sig) ----------
# Spearman over all 137,451 tested pairs is more honest because it doesn't
# condition on silencing — but we report both.
ok_all <- !is.na(sil_dv$dv_t) & !is.na(sil_dv$pct_responder)
rho_t_all  <- cor(sil_dv$pct_responder[ok_all], abs(sil_dv$dv_t[ok_all]),
                   method = "spearman")
rho_lr_all <- cor(sil_dv$pct_responder[ok_all], sil_dv$dv_logRatio[ok_all],
                   method = "spearman")
log_msg("  Across all %d tested pairs:", sum(ok_all))
log_msg("    rho(pct_responder, |dv_t|):       %.3f", rho_t_all)
log_msg("    rho(pct_responder, dv_logRatio):  %.3f", rho_lr_all)

# ---- 18-gene benchmark: DV recovery -----------------------------------------
TSG <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
         "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
         "TIMP3","NEFH","UCHL1","PITX2")
tsg_dv <- sil_dv[gene %in% TSG, .(
  best_dv_p   = min(dv_p, na.rm = TRUE),
  best_dv_FDR = min(dv_FDR, na.rm = TRUE),
  any_dv_FDR05 = any(dv_FDR < 0.05, na.rm = TRUE),
  any_silencing = any(silencing),
  n_promoter = .N
), by = gene]
setorder(tsg_dv, best_dv_p)
log_msg("\nBiomarker DV recovery (best CpG per gene):")
print(tsg_dv)
n_dv_tsg <- sum(tsg_dv$any_dv_FDR05)
n_sil_tsg <- sum(tsg_dv$any_silencing)
log_msg("  TSGs DV-recovered (FDR<0.05):       %d / 18", n_dv_tsg)
log_msg("  TSGs silencing-recovered (any cpg): %d / 18", n_sil_tsg)

# ---- Numeric sidecar for the report ----------------------------------------
out_list <- list(
  n_dv_total            = sum(dv$dv_FDR < 0.05, na.rm = TRUE),
  n_dv_pct              = 100 * mean(dv$dv_FDR < 0.05, na.rm = TRUE),
  n_dv_tumor_high_pct   = 100 * mean(dv$dv_logRatio[dv$dv_FDR < 0.05] > 0, na.rm = TRUE),
  n_dv_only             = sum(dv$dv_FDR < 0.05) - sum(dv$dv_FDR < 0.05 & dv$cpg %in% sil$cpg[sil$q1_FDR < 0.05]),
  n_sil_pairs           = n_sil,
  n_sil_dv              = n_sil_dv,
  pct_sil_dv            = 100 * n_sil_dv / n_sil,
  rho_pi_dv_t_silonly   = rho_t,
  rho_pi_dv_t_all       = rho_t_all,
  rho_pi_dv_logR_silonly = rho_lr,
  rho_pi_dv_logR_all    = rho_lr_all,
  med_dvt_subgroup      = median(abs(sg$dv_t), na.rm=TRUE),
  med_dvt_cohort_wide   = median(abs(cw$dv_t), na.rm=TRUE),
  wilcox_sg_cw_p        = sg_cw_t,
  vhl_var_tumor         = dv[cpg == "cg13672843"]$var_tumor,
  vhl_var_normal        = dv[cpg == "cg13672843"]$var_normal,
  vhl_dv_p              = dv[cpg == "cg13672843"]$dv_p,
  vhl_dv_FDR            = dv[cpg == "cg13672843"]$dv_FDR,
  n_tsg_dv_recovered    = n_dv_tsg,
  n_tsg_sil_recovered   = n_sil_tsg
)
saveRDS(out_list, file.path(OUT_DIR, "dv_numbers_for_report.rds"))

# ---- Summary text -----------------------------------------------------------
sum_lines <- c(
  sprintf("# Silencing × DV join summary (%s)", Sys.time()),
  sprintf("Silencing pairs:                 %d", n_sil),
  sprintf("Silencing AND DV (FDR<0.05):     %d (%.1f%%)", n_sil_dv, 100*n_sil_dv/n_sil),
  sprintf("Silencing + DV tumor>normal var: %d (%.1f%%)",
          n_sil_dv_tumor_high, 100*n_sil_dv_tumor_high/n_sil),
  "",
  sprintf("rho(pct_responder, |dv_t|) silencing only: %.3f", rho_t),
  sprintf("rho(pct_responder, |dv_t|) all pairs:      %.3f", rho_t_all),
  sprintf("rho(pct_responder, dv_logRatio) all pairs: %.3f", rho_lr_all),
  "",
  sprintf("Median |dv_t|: subgroup=%.2f vs cohort-wide=%.2f (Wilcoxon p=%.2e)",
          median(abs(sg$dv_t), na.rm=TRUE),
          median(abs(cw$dv_t), na.rm=TRUE), sg_cw_t),
  "",
  sprintf("18-gene biomarker recovery: %d/18 by silencing, %d/18 by DV (FDR<0.05)",
          n_sil_tsg, n_dv_tsg)
)
writeLines(sum_lines, file.path(OUT_DIR, "silencing_dv_summary.txt"))

log_msg("Done.")
