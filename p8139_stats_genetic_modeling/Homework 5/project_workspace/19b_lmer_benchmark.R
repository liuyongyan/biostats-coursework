#!/usr/bin/env Rscript
# 19b_lmer_benchmark.R
# Benchmark: how long would it actually take to refit all 137,451 MEA pairs
# with lme4::lmer (random intercept for patient) on the 342-sample matched set?
# Extrapolates from a 200-pair sample.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(lme4)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
SIL_DIR  <- file.path(WORKDIR, "silencing_out")
set.seed(8139)

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")

# Same data prep as 07/19 ------------------------------------------------------
log_msg("Loading data...")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
M2beta <- function(m) 2^m / (2^m + 1)
rna_se <- readRDS(file.path(DATA_DIR, "kirc_rnaseq_SE.rds"))
rna_cd <- as.data.frame(colData(rna_se))[, c("barcode", "patient", "sample_type")]
keep <- rna_cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
rna_se <- rna_se[, keep]; rna_cd <- rna_cd[keep, , drop = FALSE]
gene_names <- as.character(rowData(rna_se)$gene_name)
tpm <- assay(rna_se, "tpm_unstrand")
rna_keys <- paste(rna_cd$patient, rna_cd$sample_type, sep = "|")
unique_keys <- unique(rna_keys)
tpm_avg <- vapply(unique_keys, function(k) {
  ix <- which(rna_keys == k)
  if (length(ix) == 1) tpm[, ix] else rowMeans(tpm[, ix, drop = FALSE])
}, numeric(nrow(tpm)))
colnames(tpm_avg) <- unique_keys
keep_g <- !is.na(gene_names) & gene_names != ""
tpm_avg <- tpm_avg[keep_g, , drop = FALSE]
gene_lab <- gene_names[keep_g]
tpm_mat <- rowsum(tpm_avg, group = gene_lab, reorder = FALSE)
log2tpm <- log2(tpm_mat + 1)

meth_keys_all <- paste(meta$patient, meta$sample_type, sep = "|")
beta_full <- M2beta(M); rm(M); gc(verbose = FALSE)
unique_meth_keys <- unique(meth_keys_all)
if (length(unique_meth_keys) == ncol(beta_full)) {
  colnames(beta_full) <- meth_keys_all
  beta_meth <- beta_full
} else {
  beta_meth <- vapply(unique_meth_keys, function(k) {
    ix <- which(meth_keys_all == k)
    if (length(ix) == 1) beta_full[, ix] else rowMeans(beta_full[, ix, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(beta_full)))
  colnames(beta_meth) <- unique_meth_keys
}
rm(beta_full); gc(verbose = FALSE)
meta_unique <- meta[!duplicated(meth_keys_all), ]
rownames(meta_unique) <- paste(meta_unique$patient, meta_unique$sample_type, sep = "|")
matched_keys <- intersect(unique_meth_keys, colnames(log2tpm))
match_meta <- meta_unique[matched_keys, , drop = FALSE]
beta_mat <- beta_meth[, matched_keys]
log2tpm_mat <- log2tpm[, matched_keys]

# Covariate columns
match_meta$age_years <- match_meta$age_at_diagnosis / 365.25
match_meta$age_years[is.na(match_meta$age_years)] <- median(match_meta$age_years, na.rm = TRUE)
match_meta$age_c <- match_meta$age_years - mean(match_meta$age_years)
match_meta$gender[is.na(match_meta$gender) | match_meta$gender == ""] <- "unknown"
tss_table <- table(match_meta$tss)
small_tss <- names(tss_table)[tss_table < 5]
match_meta$tss_grp <- ifelse(match_meta$tss %in% small_tss, "OTHER", match_meta$tss)
match_meta$patient <- as.character(match_meta$patient)
log_msg("  Matched samples: %d (n_patients = %d)",
        nrow(match_meta), length(unique(match_meta$patient)))

# Pick 200 random pairs -------------------------------------------------------
pool <- fread(file.path(SIL_DIR, "silencing_table.csv"))
pool <- pool[!is.na(beta_coef)]
n_bench <- 200
samp_ix <- sample(seq_len(nrow(pool)), n_bench)
bench <- pool[samp_ix, .(cpg, gene)]
log_msg("  Benchmarking on %d randomly chosen (CpG, gene) pairs", n_bench)

# Time loop -------------------------------------------------------------------
fit_lmer <- function(cpg, gene) {
  df <- data.frame(
    y       = log2tpm_mat[gene, ],
    beta    = beta_mat[cpg, ],
    type    = match_meta$sample_type,
    age_c   = match_meta$age_c,
    gender  = match_meta$gender,
    tss_grp = match_meta$tss_grp,
    patient = match_meta$patient
  )
  tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(y ~ beta + type + age_c + gender + tss_grp + (1 | patient),
                 data = df, REML = TRUE,
                 control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
    )),
    error = function(e) NULL
  )
}

# Warm up (compilation / JIT effects)
log_msg("  Warm-up (3 fits)...")
for (i in 1:3) fit_lmer(bench$cpg[i], bench$gene[i])

log_msg("  Timing %d fits...", n_bench)
t0 <- Sys.time()
ok <- 0L
for (i in seq_len(n_bench)) {
  m <- fit_lmer(bench$cpg[i], bench$gene[i])
  if (!is.null(m)) ok <- ok + 1L
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
per_fit <- elapsed / n_bench
N_TOTAL <- 137451L
total_sec <- per_fit * N_TOTAL

cat("\n================================================================\n")
cat(" lmer benchmark (random intercept for patient, 342 samples)\n")
cat("================================================================\n")
cat(sprintf("  Pairs benchmarked:        %d\n", n_bench))
cat(sprintf("  Successful fits:          %d\n", ok))
cat(sprintf("  Total elapsed:            %.2f s\n", elapsed))
cat(sprintf("  Per-fit:                  %.4f s\n", per_fit))
cat(sprintf("  Per-fit (ms):             %.1f ms\n", per_fit * 1000))
cat("\n  Extrapolated for the full screen (137,451 pairs):\n")
cat(sprintf("    Single-thread:          %.0f s = %.1f min = %.2f hours\n",
            total_sec, total_sec / 60, total_sec / 3600))
cat(sprintf("    With 8 cores (parallel):%.1f min = %.2f hours\n",
            total_sec / 8 / 60, total_sec / 8 / 3600))
cat("\n  Compared to the existing partial-correlation pipeline:\n")
cat("    QR-based vectorized fit:  ~0.1 s for all 137,451 pairs\n")
cat(sprintf("    Slowdown factor:          %.0fx (single-thread)\n",
            total_sec / 0.1))
cat("================================================================\n")

writeLines(c(
  sprintf("n_bench=%d", n_bench),
  sprintf("ok=%d", ok),
  sprintf("elapsed_sec=%.4f", elapsed),
  sprintf("per_fit_sec=%.4f", per_fit),
  sprintf("extrapolated_total_sec=%.1f", total_sec),
  sprintf("extrapolated_total_hours=%.2f", total_sec / 3600)
), file.path(WORKDIR, "sensitivity_tumoronly", "lmer_benchmark.txt"))
