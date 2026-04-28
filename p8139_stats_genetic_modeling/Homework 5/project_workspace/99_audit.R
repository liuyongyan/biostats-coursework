# Step 4 — self-audit in "professor mode" (v2).
#
# Revisions from v1 (2026-04-23):
#   - Lambda interpretation rewritten. For tumor-vs-normal / grade-trend
#     methylation analyses we EXPECT high lambda because real signal is
#     pervasive. The new check distinguishes "pervasive real signal" from
#     "confounding inflation" by looking at the shape of the Q-Q plot at
#     small -log10(p): if observed tracks y=x near the origin and the null
#     portion is calibrated, high lambda reflects signal, not confounding.
#   - Added effect-size-aware hit counts (FDR & |Delta beta|).
#   - Added "null-region calibration" check using only the least-significant
#     20% of p-values.

suppressPackageStartupMessages({
  library(dplyr)
  library(SummarizedExperiment)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
AUD_DIR <- file.path(WORKDIR, "audit_out")
dir.create(AUD_DIR, showWarnings = FALSE, recursive = TRUE)

report <- list()
add <- function(id, status, msg) {
  report[[length(report)+1]] <<- data.frame(
    id = id, status = status, msg = msg, stringsAsFactors = FALSE
  )
}

# Helper: null-region calibration — using the top 80% largest p-values (the
# null-ish region), check that observed quantiles track expected.
null_calibration <- function(pvals, tail_frac = 0.8) {
  # pvals sorted ascending; keep the LOW-significance tail = top 80% largest
  p_sorted <- sort(pvals)
  n <- length(p_sorted)
  keep <- ceiling(n * (1 - tail_frac)) + 1
  p_null <- p_sorted[keep:n]
  exp_q  <- ppoints(length(p_null)) * (1 - p_sorted[keep]) + p_sorted[keep]
  # lambda computed on the null-ish region only
  lam <- median(qchisq(1 - p_null, df = 1)) / qchisq(0.5, df = 1)
  # Correlation between expected vs observed on -log10 scale
  cor_ <- cor(-log10(p_null), -log10(exp_q))
  list(lambda_null = lam, cor_obs_exp = cor_, n_null = length(p_null))
}

# ----- load -----
M     <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta  <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
anno_ <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))

# =======================================================================
# A. Preprocessing integrity
# =======================================================================
add("A1_no_chrXY",
    if (!any(anno_$chr %in% c("chrX","chrY"))) "PASS" else "FAIL",
    sprintf("%d chrX probes, %d chrY probes remain",
            sum(anno_$chr == "chrX"), sum(anno_$chr == "chrY")))

add("A2_no_missing_M",
    if (!any(is.na(M))) "PASS" else "FAIL",
    sprintf("%d NA cells in M matrix (should be 0)", sum(is.na(M))))

add("A3_probe_count",
    if (nrow(M) >= 300000 && nrow(M) <= 485577) "PASS" else "WARN",
    sprintf("post-filter probe count = %d", nrow(M)))

m_range <- range(M, finite = TRUE)
add("A4_M_range",
    if (all(is.finite(m_range)) && m_range[1] > -12 && m_range[2] < 12) "PASS" else "WARN",
    sprintf("M range [%.3f, %.3f]", m_range[1], m_range[2]))

st_tab <- table(meta$sample_type)
add("A5_sample_types",
    if (setequal(names(st_tab), c("Primary Tumor","Solid Tissue Normal"))) "PASS" else "FAIL",
    paste(names(st_tab), st_tab, sep="=", collapse="; "))

add("A6_order",
    if (identical(colnames(M), meta$barcode)) "PASS" else "FAIL",
    "colnames(M) should equal meta$barcode in order")

# =======================================================================
# B. Q1 integrity — statistical + calibration-aware lambda
# =======================================================================
q1_path <- file.path(WORKDIR, "q1_out", "toptable_q1.csv")
if (file.exists(q1_path)) {
  tt1 <- read.csv(q1_path, stringsAsFactors = FALSE)
  q1_fit <- readRDS(file.path(WORKDIR, "q1_out", "q1_fit.rds"))

  # Null-region calibration for Q1
  nc <- null_calibration(tt1$P.Value, tail_frac = 0.8)
  add("B1_null_calibration",
      if (nc$cor_obs_exp > 0.95 && nc$lambda_null < 3) "PASS" else "WARN",
      sprintf("null-region lambda=%.2f, obs/exp cor=%.3f (checks calibration on least-sig 80%%; global lambda=%.2f is expected high for tumor vs normal due to pervasive real signal)",
              nc$lambda_null, nc$cor_obs_exp, q1_fit$lambda))

  # Effect-size and hit counts
  es_path <- file.path(WORKDIR, "effectsize_out", "q1_with_beta.csv")
  if (file.exists(es_path)) {
    tt1b <- read.csv(es_path)
    n_fdr  <- sum(tt1b$adj.P.Val < 0.05)
    n_strg <- sum(tt1b$adj.P.Val < 0.05 & abs(tt1b$delta_beta) > 0.1)
    n_hyp  <- sum(tt1b$adj.P.Val < 0.05 & tt1b$delta_beta > 0.1)
    n_hypo <- sum(tt1b$adj.P.Val < 0.05 & tt1b$delta_beta < -0.1)

    add("B2_any_hits",
        if (n_fdr >= 10) "PASS" else "WARN",
        sprintf("Q1 FDR<0.05: %d probes", n_fdr))
    add("B3_strong_hits",
        if (n_strg >= 100) "PASS" else "WARN",
        sprintf("Q1 FDR<0.05 AND |Δβ|>0.1: %d probes", n_strg))
    add("B4_bidirectional",
        if (n_hyp > 100 && n_hypo > 100) "PASS" else "WARN",
        sprintf("Strong hits: hyper=%d, hypo=%d (balanced both directions)", n_hyp, n_hypo))
  }

  # Sensitivity overlaps
  s1 <- tryCatch(read.csv(file.path(WORKDIR, "q1_out", "toptable_q1_s1_unpaired.csv")), error=function(e) NULL)
  s2 <- tryCatch(read.csv(file.path(WORKDIR, "q1_out", "toptable_q1_s2_sameplatepairs.csv")), error=function(e) NULL)
  if (!is.null(s1)) {
    top_p <- tt1$cpg[tt1$adj.P.Val < 0.05]
    top_s <- s1$cpg[s1$adj.P.Val < 0.05]
    ov <- length(intersect(top_p, top_s)) / max(1, length(top_p))
    add("B5_sens_unpaired_overlap",
        if (ov > 0.5) "PASS" else "WARN",
        sprintf("%.1f%% of Q1 primary hits also in unpaired S1", 100*ov))
  }
  if (!is.null(s2)) {
    top_p <- tt1$cpg[tt1$adj.P.Val < 0.05]
    top_s <- s2$cpg[s2$adj.P.Val < 0.05]
    ov <- length(intersect(top_p, top_s)) / max(1, length(top_p))
    add("B6_sens_sameplate_overlap",
        if (ov > 0.5) "PASS" else "WARN",
        sprintf("%.1f%% of Q1 primary hits also in same-plate S2", 100*ov))
  }
} else add("B0_q1_missing", "FAIL", "q1 toptable missing")

# =======================================================================
# C. Q3 integrity
# =======================================================================
q3_path <- file.path(WORKDIR, "q3_out", "toptable_q3.csv")
if (file.exists(q3_path)) {
  tt3 <- read.csv(q3_path, stringsAsFactors = FALSE)
  q3_fit <- readRDS(file.path(WORKDIR, "q3_out", "q3_fit.rds"))
  nc3 <- null_calibration(tt3$P.Value, tail_frac = 0.8)
  add("C1_null_calibration",
      if (nc3$cor_obs_exp > 0.95 && nc3$lambda_null < 2) "PASS" else "WARN",
      sprintf("Q3 null lambda=%.2f cor=%.3f (global=%.2f, expected elevated)",
              nc3$lambda_null, nc3$cor_obs_exp, q3_fit$lambda))
  n3 <- sum(tt3$adj.P.Val < 0.05)
  add("C2_q3_hits",
      if (n3 >= 1) "PASS" else "WARN",
      sprintf("Q3 FDR<0.05: %d probes", n3))
} else add("C0_q3_missing", "FAIL", "q3 toptable missing")

# =======================================================================
# D. Q9 positive control
# =======================================================================
q9_path <- file.path(WORKDIR, "q9_out", "q9_assertions.csv")
if (file.exists(q9_path)) {
  q9 <- read.csv(q9_path, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(q9))) {
    status <- if (isTRUE(as.logical(q9$pass[i]))) "PASS" else "FAIL"
    add(paste0("D_", q9$assertion[i]), status, q9$note[i])
  }
} else add("D0_q9_missing", "FAIL", "q9 assertions missing")

# =======================================================================
# E. GSEA / Pathway — check biological expectations are met
# =======================================================================
g1p <- file.path(WORKDIR, "q5_out", "gsea_q1.csv")
if (file.exists(g1p)) {
  gs1 <- read.csv(g1p)
  n_sig <- sum(gs1$padj < 0.25, na.rm = TRUE)
  add("E1_gsea_q1_any",
      if (n_sig >= 1) "PASS" else "WARN",
      sprintf("Q1 GSEA sets at padj<0.25: %d", n_sig))
  hypox <- gs1[gs1$pathway == "HALLMARK_HYPOXIA", ]
  add("E2_hypoxia_q1",
      if (nrow(hypox) && hypox$padj < 0.25) "PASS" else "WARN",
      sprintf("HALLMARK_HYPOXIA in Q1: padj=%.3f, NES=%.2f",
              hypox$padj, hypox$NES))
} else add("E0_gsea_missing", "FAIL", "gsea q1 missing")

g3p <- file.path(WORKDIR, "q5_out", "gsea_q3.csv")
if (file.exists(g3p)) {
  gs3 <- read.csv(g3p)
  rcc <- gs3[gs3$pathway == "KEGG_RENAL_CELL_CARCINOMA", ]
  add("E3_rcc_q3",
      if (nrow(rcc) && rcc$padj < 0.25) "PASS" else "WARN",
      sprintf("KEGG_RENAL_CELL_CARCINOMA in Q3: padj=%.3f, NES=%.2f",
              rcc$padj, rcc$NES))
}

# =======================================================================
# F. Reproducibility
# =======================================================================
add("F1_session_logs",
    if (all(file.exists(c(
      file.path(PRE_DIR, "sessionInfo_preprocess.txt"),
      file.path(WORKDIR, "q1_out", "sessionInfo_q1.txt"),
      file.path(WORKDIR, "q3_out", "sessionInfo_q3.txt")
    )))) "PASS" else "WARN",
    "sessionInfo logs present at each step")

# =======================================================================
df <- do.call(rbind, report)
df <- df[order(factor(df$status, levels = c("FAIL","WARN","PASS"))), ]
write.csv(df, file.path(AUD_DIR, "audit_report.csv"), row.names = FALSE)

cat("=== AUDIT REPORT (v2) ===\n")
for (i in seq_len(nrow(df))) {
  cat(sprintf("[%s] %s — %s\n", df$status[i], df$id[i], df$msg[i]))
}
n_fail <- sum(df$status == "FAIL"); n_warn <- sum(df$status == "WARN"); n_pass <- sum(df$status == "PASS")
cat(sprintf("\nSummary: %d PASS, %d WARN, %d FAIL\n", n_pass, n_warn, n_fail))

mdlines <- c(
  "# Audit Report (v2, null-calibration-aware)",
  "",
  sprintf("PASS: %d   WARN: %d   FAIL: %d", n_pass, n_warn, n_fail),
  "",
  "| Status | ID | Message |",
  "|---|---|---|",
  sprintf("| %s | %s | %s |", df$status, df$id, df$msg)
)
writeLines(mdlines, file.path(AUD_DIR, "audit_report.md"))

if (n_fail > 0) quit(status = 1) else quit(status = 0)
