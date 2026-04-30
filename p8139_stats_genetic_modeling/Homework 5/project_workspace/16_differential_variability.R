#!/usr/bin/env Rscript
# 16_differential_variability.R
#
# Differential variability (DV) test on TCGA-KIRC 450K methylation.
#
# Rationale: in cancer, hyper-silencing of TSGs often shows as a tumor-only
# variance increase rather than a mean shift — a subset of tumors hypermethylate
# the locus while the rest cluster with normals (Teschendorff & Widschwendter
# 2012; Teschendorff 2016 iEVORA). The responder fraction pi_resp used in
# 07_methylation_silencing.R is a non-parametric, threshold-based proxy for
# this same pattern; this script connects it to the established DV framework.
#
# Method: missMethyl::varFit (Phipson & Oshlack 2014) — Levene-style absolute
# deviation from group median, fit by limma + empirical Bayes shrinkage of the
# residual variance, returning a moderated DV statistic per CpG. Cross-checked
# against Bartlett's test on the VHL anchor CpG.
#
# Inputs:  preprocess_out/{M_matrix.rds, sample_metadata.rds}
#          q1_out/toptable_q1.csv   (for cross-tab with DM hits)
# Outputs: dv_out/dv_table.csv      (per-CpG DV stat + FDR + tumor/normal var)
#          dv_out/dv_summary.txt    (counts, VHL Bartlett, etc.)

suppressPackageStartupMessages({
  library(limma)
  library(missMethyl)
  library(data.table)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
Q1_DIR   <- file.path(WORKDIR, "q1_out")
OUT_DIR  <- file.path(WORKDIR, "dv_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")

# ---- 1. Load data ----------------------------------------------------------
log_msg("[1/5] Loading M-matrix and metadata...")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
ann  <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
stopifnot(all(colnames(M) == meta$barcode))
log_msg("  M: %d CpGs x %d samples", nrow(M), ncol(M))

# Use Primary Tumor + Solid Tissue Normal only
keep <- meta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
M    <- M[, keep]
meta <- meta[keep, , drop = FALSE]
log_msg("  After PT/STN filter: %d samples (%d tumor, %d normal)",
        ncol(M), sum(meta$sample_type == "Primary Tumor"),
        sum(meta$sample_type == "Solid Tissue Normal"))

# ---- 2. Variance per group (unmoderated, for context) ----------------------
log_msg("[2/5] Computing per-group sample variances on M-scale...")
is_tumor  <- meta$sample_type == "Primary Tumor"
is_normal <- meta$sample_type == "Solid Tissue Normal"
var_tumor  <- matrixStats::rowVars(M[, is_tumor,  drop = FALSE])
var_normal <- matrixStats::rowVars(M[, is_normal, drop = FALSE])
log_ratio  <- log2(var_tumor / var_normal)
log_msg("  log2 var(tumor)/var(normal) median: %.3f, IQR: [%.3f, %.3f]",
        median(log_ratio, na.rm = TRUE),
        quantile(log_ratio, 0.25, na.rm = TRUE),
        quantile(log_ratio, 0.75, na.rm = TRUE))
log_msg("  Probes with var(tumor) > var(normal): %d / %d (%.1f%%)",
        sum(log_ratio > 0, na.rm = TRUE), length(log_ratio),
        100 * mean(log_ratio > 0, na.rm = TRUE))

# ---- 3. missMethyl::varFit  (Phipson & Oshlack 2014) -----------------------
log_msg("[3/5] Running missMethyl::varFit (moderated DV test)...")
group  <- factor(ifelse(is_tumor, "Tumor", "Normal"), levels = c("Normal", "Tumor"))
design <- model.matrix(~ group)
t0 <- Sys.time()
fit_dv <- varFit(M, design = design, coef = c(1, 2))
log_msg("  varFit done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))

dv_top <- topVar(fit_dv, coef = 2, number = nrow(M), sort = FALSE)
dv_top$cpg <- rownames(dv_top)
log_msg("  topVar columns: %s", paste(colnames(dv_top), collapse = ", "))

# Build canonical DV table
dv_dt <- data.table(
  cpg          = dv_top$cpg,
  dv_logRatio  = dv_top$LogVarRatio,    # log fold-change of variances (tumor vs normal)
  dv_t         = dv_top$t,
  dv_p         = dv_top$P.Value,
  dv_FDR       = dv_top$Adj.P.Value,
  var_tumor    = var_tumor[match(dv_top$cpg, names(var_tumor))],
  var_normal   = var_normal[match(dv_top$cpg, names(var_normal))]
)
setorder(dv_dt, dv_p)

n_dv <- sum(dv_dt$dv_FDR < 0.05, na.rm = TRUE)
n_dv_hyper <- sum(dv_dt$dv_FDR < 0.05 & dv_dt$dv_logRatio > 0, na.rm = TRUE)
n_dv_hypo  <- sum(dv_dt$dv_FDR < 0.05 & dv_dt$dv_logRatio < 0, na.rm = TRUE)
log_msg("  DV FDR<0.05: %d (%.1f%%) — %d tumor>normal var, %d normal>tumor var",
        n_dv, 100 * n_dv / nrow(dv_dt), n_dv_hyper, n_dv_hypo)

fwrite(dv_dt, file.path(OUT_DIR, "dv_table.csv"))

# ---- 4. Cross-check: Bartlett's test on VHL anchor + 5k random probes ------
log_msg("[4/5] Cross-checking DV with Bartlett's test (VHL anchor + 5k random probes)...")

bartlett_safe <- function(y, g) {
  out <- tryCatch(bartlett.test(y, g)$p.value,
                  error = function(e) NA_real_,
                  warning = function(w) bartlett.test(y, g)$p.value)
  out
}

# VHL CpG
vhl_cpg <- "cg13672843"
if (vhl_cpg %in% rownames(M)) {
  bart_vhl <- bartlett_safe(M[vhl_cpg, ], group)
  vhl_dv <- dv_dt[cpg == vhl_cpg]
  log_msg("  VHL/cg13672843: var_tumor = %.3f, var_normal = %.3f, log2 ratio = %.2f",
          vhl_dv$var_tumor, vhl_dv$var_normal,
          log2(vhl_dv$var_tumor / vhl_dv$var_normal))
  log_msg("  VHL/cg13672843: varFit DV p = %.2e (FDR = %.2e); Bartlett p = %.2e",
          vhl_dv$dv_p, vhl_dv$dv_FDR, bart_vhl)
} else {
  log_msg("  WARNING: VHL CpG not in M-matrix")
  bart_vhl <- NA_real_
}

# Random sample of 5,000 CpGs for Bartlett vs varFit Spearman check
set.seed(8139)
sub_ix <- sample.int(nrow(M), min(5000, nrow(M)))
bart_p <- vapply(sub_ix, function(i) bartlett_safe(M[i, ], group),
                 numeric(1))
sub_cpg <- rownames(M)[sub_ix]
sub_dv  <- dv_dt[match(sub_cpg, dv_dt$cpg), dv_p]
ok_ix   <- is.finite(bart_p) & is.finite(sub_dv) & bart_p > 0 & sub_dv > 0
spear_bart_var <- cor(-log10(bart_p[ok_ix]), -log10(sub_dv[ok_ix]),
                      method = "spearman")
log_msg("  Spearman rho(-log10 Bartlett p, -log10 varFit p) on %d probes: %.3f",
        sum(ok_ix), spear_bart_var)

# ---- 5. Cross-tab against DM (Q1) and write summary ------------------------
log_msg("[5/5] Cross-tabulating DV vs DM hits...")
dm <- fread(file.path(Q1_DIR, "toptable_q1.csv"))
dm_sig <- dm$cpg[dm$adj.P.Val < 0.05]
dv_sig <- dv_dt$cpg[dv_dt$dv_FDR < 0.05]
both   <- intersect(dm_sig, dv_sig)
dv_only <- setdiff(dv_sig, dm_sig)
dm_only <- setdiff(dm_sig, dv_sig)
n_total <- nrow(dv_dt)

log_msg("  DM only (FDR<0.05):     %d", length(dm_only))
log_msg("  DV only (FDR<0.05):     %d", length(dv_only))
log_msg("  DM and DV (FDR<0.05):   %d", length(both))
log_msg("  Neither:                %d", n_total - length(union(dm_sig, dv_sig)))

# Among DV-significant, fraction with tumor variance > normal (cancer signature)
dv_sig_dt <- dv_dt[dv_FDR < 0.05]
n_tumor_high <- sum(dv_sig_dt$dv_logRatio > 0, na.rm = TRUE)
log_msg("  Among DV-sig: %d (%.1f%%) have var(tumor) > var(normal)",
        n_tumor_high, 100 * n_tumor_high / nrow(dv_sig_dt))

# Summary file
summary_lines <- c(
  sprintf("# Differential variability summary (%s)", Sys.time()),
  sprintf("Samples: %d tumor + %d normal", sum(is_tumor), sum(is_normal)),
  sprintf("Probes tested: %d", nrow(dv_dt)),
  sprintf("DV FDR < 0.05: %d (%.1f%%)", n_dv, 100 * n_dv / nrow(dv_dt)),
  sprintf("  var(tumor) > var(normal): %d", n_dv_hyper),
  sprintf("  var(normal) > var(tumor): %d", n_dv_hypo),
  sprintf("DM FDR < 0.05: %d", length(dm_sig)),
  sprintf("DM AND DV (FDR<0.05): %d", length(both)),
  sprintf("DM only:               %d", length(dm_only)),
  sprintf("DV only:               %d", length(dv_only)),
  "",
  sprintf("VHL/cg13672843: var_tumor=%.3f, var_normal=%.3f, log2 ratio=%.2f",
          dv_dt[cpg == vhl_cpg]$var_tumor, dv_dt[cpg == vhl_cpg]$var_normal,
          log2(dv_dt[cpg == vhl_cpg]$var_tumor / dv_dt[cpg == vhl_cpg]$var_normal)),
  sprintf("VHL varFit DV p=%.2e (FDR=%.2e); Bartlett p=%.2e",
          dv_dt[cpg == vhl_cpg]$dv_p, dv_dt[cpg == vhl_cpg]$dv_FDR, bart_vhl),
  "",
  sprintf("Spearman rho(-log10 Bartlett p, -log10 varFit p) on 5k probes: %.3f",
          spear_bart_var)
)
writeLines(summary_lines, file.path(OUT_DIR, "dv_summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_dv.txt"))
log_msg("Done. Outputs in %s", OUT_DIR)
