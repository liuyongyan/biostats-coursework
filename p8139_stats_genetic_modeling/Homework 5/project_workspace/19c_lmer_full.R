#!/usr/bin/env Rscript
# 19c_lmer_full.R
# Full mixed-effects refit of all 137,451 (CpG, gene) MEA pairs using
# lme4::lmer with a patient-level random intercept. Triangulates against
# the pooled (07) and tumor-only (19) fits.
#
# Model:
#   log2(TPM_gene + 1) = alpha + beta * beta_CpG + gamma' Z + u_patient + eps
#   Z = (sample_type, age_c, gender, tss_grp);  u_patient ~ N(0, sigma_p^2)
# Inference on beta: Wald z-test (asymptotic; at n=342 essentially identical
# to Satterthwaite t).

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(lme4)
  library(parallel)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
Q1_DIR   <- file.path(WORKDIR, "q1_out")
SIL_DIR  <- file.path(WORKDIR, "silencing_out")
TO_DIR   <- file.path(WORKDIR, "sensitivity_tumoronly")
OUT_DIR  <- file.path(WORKDIR, "sensitivity_lmer")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(8139)

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")

# --- Step 1: load data (mirrors 07 / 19) --------------------------------------
log_msg("[1/5] Loading data...")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
ann  <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
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

# Covariates
match_meta$age_years <- match_meta$age_at_diagnosis / 365.25
match_meta$age_years[is.na(match_meta$age_years)] <- median(match_meta$age_years, na.rm = TRUE)
match_meta$age_c <- match_meta$age_years - mean(match_meta$age_years)
match_meta$gender[is.na(match_meta$gender) | match_meta$gender == ""] <- "unknown"
tss_table <- table(match_meta$tss)
small_tss <- names(tss_table)[tss_table < 5]
match_meta$tss_grp <- ifelse(match_meta$tss %in% small_tss, "OTHER", match_meta$tss)
match_meta$patient <- as.character(match_meta$patient)

log_msg("  Samples: %d (n_patients = %d, paired = %d)",
        nrow(match_meta),
        length(unique(match_meta$patient)),
        sum(table(match_meta$patient) >= 2))

# --- Step 2: promoter pair set (mirrors 07) ------------------------------------
log_msg("[2/5] Building promoter pairs...")
PROMOTER_GROUPS <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
ann_dt <- data.table(
  cpg   = rownames(ann),
  names = as.character(ann$UCSC_RefGene_Name),
  groups = as.character(ann$UCSC_RefGene_Group)
)[names != "" & !is.na(names)]
ann_dt[, names_list  := strsplit(names,  ";", fixed = TRUE)]
ann_dt[, groups_list := strsplit(groups, ";", fixed = TRUE)]
ann_dt[, ok_len := lengths(names_list) == lengths(groups_list)]
ann_dt <- ann_dt[ok_len == TRUE]
expanded <- ann_dt[, .(gene = unlist(names_list), group = unlist(groups_list)), by = cpg]
expanded <- expanded[group %in% PROMOTER_GROUPS]
PROMOTER_RANK <- c("TSS200" = 1, "TSS1500" = 2, "5'UTR" = 3, "1stExon" = 4)
expanded[, group_rank := PROMOTER_RANK[group]]
setorder(expanded, cpg, gene, group_rank)
expanded <- expanded[!duplicated(expanded[, .(cpg, gene)])]
expanded <- expanded[gene %in% rownames(log2tpm_mat) & cpg %in% rownames(beta_mat)]
log_msg("  Pairs to fit: %d", nrow(expanded))

pairs_cpg  <- expanded$cpg
pairs_gene <- expanded$gene
N_PAIRS    <- nrow(expanded)

# --- Step 3: per-pair lmer fit ------------------------------------------------
meta_df <- data.frame(
  type    = factor(match_meta$sample_type),
  age_c   = match_meta$age_c,
  gender  = factor(match_meta$gender),
  tss_grp = factor(match_meta$tss_grp),
  patient = factor(match_meta$patient)
)

fit_one <- function(i) {
  df <- meta_df
  df$y    <- log2tpm_mat[pairs_gene[i], ]
  df$beta <- beta_mat[pairs_cpg[i], ]
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(y ~ beta + type + age_c + gender + tss_grp + (1 | patient),
                 data = df, REML = TRUE)
    )),
    error = function(e) NULL
  )
  if (is.null(m)) return(c(NA_real_, NA_real_, NA_real_))
  bh <- tryCatch(lme4::fixef(m), error = function(e) NULL)
  vm <- tryCatch(as.matrix(lme4::vcov.merMod(m)), error = function(e) NULL)
  if (is.null(bh) || is.null(vm) || !"beta" %in% names(bh)) {
    return(c(NA_real_, NA_real_, NA_real_))
  }
  est <- unname(bh["beta"])
  se  <- sqrt(vm["beta", "beta"])
  if (!is.finite(se) || se <= 0) return(c(est, NA_real_, NA_real_))
  z   <- est / se
  pv  <- 2 * pnorm(-abs(z))
  c(est, se, pv)
}

n_cores <- max(1, parallel::detectCores() - 1)
log_msg("[3/5] Fitting %d lmer models on %d cores...", N_PAIRS, n_cores)

t0 <- Sys.time()
res_list <- parallel::mclapply(seq_len(N_PAIRS), fit_one,
                               mc.cores = n_cores, mc.preschedule = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
log_msg("  Done in %.0f s = %.1f min  (%.2f ms/fit avg)",
        elapsed, elapsed/60, 1000*elapsed/N_PAIRS)

normalize <- function(x) {
  if (is.numeric(x) && length(x) == 3) return(x)
  c(NA_real_, NA_real_, NA_real_)
}
res_list <- lapply(res_list, normalize)
res_mat <- do.call(rbind, res_list)
stopifnot(ncol(res_mat) == 3, nrow(res_mat) == N_PAIRS)
n_fail <- sum(is.na(res_mat[, 1]))
log_msg("  Convergence failures / NA fits: %d (%.2f%%)",
        n_fail, 100*n_fail/N_PAIRS)

# --- Step 4: build output table + silencing flags -----------------------------
log_msg("[4/5] Computing FDR + silencing flags...")
q1 <- fread(file.path(Q1_DIR, "toptable_q1.csv"))
q1_sub <- q1[, .(cpg, q1_logFC = logFC, q1_FDR = adj.P.Val)]

out <- data.table(
  cpg            = pairs_cpg,
  gene           = pairs_gene,
  promoter_group = expanded$group,
  beta_coef_lmer = res_mat[, 1],
  se_lmer        = res_mat[, 2],
  assoc_p_lmer   = res_mat[, 3]
)
out[, assoc_FDR_lmer := p.adjust(assoc_p_lmer, method = "BH")]
out <- merge(out, q1_sub, by = "cpg", all.x = TRUE)
out[, silencing_lmer := !is.na(q1_FDR) & !is.na(assoc_FDR_lmer) &
                        q1_FDR < 0.05 & assoc_FDR_lmer < 0.05 &
                        q1_logFC > 0 & beta_coef_lmer < 0]
out[, activation_lmer := !is.na(q1_FDR) & !is.na(assoc_FDR_lmer) &
                         q1_FDR < 0.05 & assoc_FDR_lmer < 0.05 &
                         q1_logFC < 0 & beta_coef_lmer > 0]
fwrite(out, file.path(OUT_DIR, "silencing_table_lmer.csv.gz"))

# --- Step 5: triangulate against pooled and tumor-only -------------------------
log_msg("[5/5] Triangulating: lmer vs pooled vs tumor-only...")
pooled <- fread(file.path(SIL_DIR, "silencing_table.csv"))
to     <- fread(file.path(TO_DIR,  "silencing_table_tumoronly.csv"))

triple <- merge(
  pooled[, .(cpg, gene, beta_p = beta_coef, p_p = assoc_p,
             FDR_p = assoc_FDR, sil_p = silencing, act_p = activation)],
  to[, .(cpg, gene, beta_to = beta_coef_to, p_to = assoc_p_to,
         FDR_to = assoc_FDR_to, sil_to = silencing_to, act_to = activation_to)],
  by = c("cpg", "gene"))
triple <- merge(
  triple,
  out[, .(cpg, gene, beta_lm = beta_coef_lmer, p_lm = assoc_p_lmer,
          FDR_lm = assoc_FDR_lmer, sil_lm = silencing_lmer, act_lm = activation_lmer)],
  by = c("cpg", "gene"))
fwrite(triple, file.path(OUT_DIR, "triple_compare.csv.gz"))

key_str <- function(d, col) d[get(col) == TRUE, paste(cpg, gene, sep = "|")]
sil_p  <- key_str(triple, "sil_p")
sil_to <- key_str(triple, "sil_to")
sil_lm <- key_str(triple, "sil_lm")

jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))

valid_pl <- triple[!is.na(beta_p) & !is.na(beta_lm)]
valid_tl <- triple[!is.na(beta_to) & !is.na(beta_lm)]

vhl <- triple[cpg == "cg13672843" & gene == "VHL"]

TSG <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
         "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
         "TIMP3","NEFH","UCHL1","PITX2")
tsg_tab <- triple[gene %in% TSG, .(
  pooled     = any(sil_p,  na.rm=TRUE),
  tumor_only = any(sil_to, na.rm=TRUE),
  lmer       = any(sil_lm, na.rm=TRUE)
), by = gene]
tsg_tab[, all_three_agree := pooled == tumor_only & tumor_only == lmer]
fwrite(tsg_tab, file.path(OUT_DIR, "tsg_three_way_compare.csv"))

summary_lines <- c(
  "===================================================================",
  " THREE-WAY MEA SENSITIVITY: pooled  vs  tumor-only  vs  lmer",
  "===================================================================",
  sprintf(" Pairs in all three fits:                  %d", nrow(triple)),
  sprintf(" lmer fit failures (NA beta):              %d", n_fail),
  sprintf(" lmer wall-clock (parallel, %d cores):     %.0f s = %.1f min",
          n_cores, elapsed, elapsed/60),
  "",
  " (a) beta_MEA agreement (Spearman / Pearson)",
  sprintf("     pooled vs lmer:          rho_S=%.4f  rho_P=%.4f",
          cor(valid_pl$beta_p,  valid_pl$beta_lm,  method = "spearman"),
          cor(valid_pl$beta_p,  valid_pl$beta_lm,  method = "pearson")),
  sprintf("     tumor-only vs lmer:      rho_S=%.4f  rho_P=%.4f",
          cor(valid_tl$beta_to, valid_tl$beta_lm, method = "spearman"),
          cor(valid_tl$beta_to, valid_tl$beta_lm, method = "pearson")),
  "",
  " (b) Silencing list size + Jaccard pair-wise",
  sprintf("     Pooled        n_silencing = %d", length(sil_p)),
  sprintf("     Tumor-only    n_silencing = %d", length(sil_to)),
  sprintf("     lmer          n_silencing = %d", length(sil_lm)),
  sprintf("     Jaccard pooled vs lmer:        %.3f", jac(sil_p,  sil_lm)),
  sprintf("     Jaccard tumor-only vs lmer:    %.3f", jac(sil_to, sil_lm)),
  sprintf("     Jaccard pooled vs tumor-only:  %.3f", jac(sil_p,  sil_to)),
  sprintf("     Triple intersection:           %d",
          length(Reduce(intersect, list(sil_p, sil_to, sil_lm)))),
  "",
  " (c) VHL / cg13672843",
  sprintf("     pooled    : beta=%+.3f  p=%.2e  FDR=%.2e  sil=%s",
          vhl$beta_p,  vhl$p_p,  vhl$FDR_p,  vhl$sil_p),
  sprintf("     tumor-only: beta=%+.3f  p=%.2e  FDR=%.2e  sil=%s",
          vhl$beta_to, vhl$p_to, vhl$FDR_to, vhl$sil_to),
  sprintf("     lmer      : beta=%+.3f  p=%.2e  FDR=%.2e  sil=%s",
          vhl$beta_lm, vhl$p_lm, vhl$FDR_lm, vhl$sil_lm),
  "",
  " (d) Predefined TSG silencing recovery (any-CpG level)",
  sprintf("     Pooled:     %d / %d", sum(tsg_tab$pooled),     nrow(tsg_tab)),
  sprintf("     Tumor-only: %d / %d", sum(tsg_tab$tumor_only), nrow(tsg_tab)),
  sprintf("     lmer:       %d / %d", sum(tsg_tab$lmer),       nrow(tsg_tab)),
  sprintf("     Genes with three-way agreement: %d / %d",
          sum(tsg_tab$all_three_agree), nrow(tsg_tab)),
  "==================================================================="
)
writeLines(summary_lines, file.path(OUT_DIR, "summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n--- TSG three-way table ---\n")
print(tsg_tab)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
log_msg("Done. Outputs in %s", OUT_DIR)
