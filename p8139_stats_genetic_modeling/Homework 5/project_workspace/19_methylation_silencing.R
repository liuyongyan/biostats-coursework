#!/usr/bin/env Rscript
# 19_methylation_silencing.R
# Identify CpGs whose promoter methylation is associated with transcriptional
# silencing of the cognate gene in TCGA-KIRC.
#
# Statistical model (per (CpG, gene) pair):
#   log2(TPM_gene + 1) ~ beta_CpG + sample_type + age + gender + TSS
# Tested across all (patient, sample_type) keys with both modalities present.
# Implemented as partial-correlation residualization for vectorized speed.
#
# Silencing definition:
#   Q1 paired-DM FDR < 0.05  AND
#   association FDR < 0.05   AND
#   sign(Q1 logFC) opposite to sign(beta_coef)   (hyper -> mRNA down, or hypo -> mRNA up)
#
# Pattern partition (within silencing CpGs):
#   pct_responder >= 0.3 -> cohort_wide
#   pct_responder <  0.3 -> subgroup
# pct_responder = fraction of tumor samples with |beta_tumor - mean(beta_normal)| > 0.2

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
Q1_DIR   <- file.path(WORKDIR, "q1_out")
OUT_DIR  <- file.path(WORKDIR, "silencing_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(8139)

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")

# === Step 1: Load methylation =================================================
log_msg("[1/8] Loading methylation matrix and metadata...")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
ann  <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
log_msg("  M-matrix: %d CpGs x %d samples", nrow(M), ncol(M))
M2beta <- function(m) 2^m / (2^m + 1)

# === Step 2: Load RNA-seq and aggregate by (patient, sample_type) =============
log_msg("[2/8] Loading RNA-seq SummarizedExperiment...")
rna_se <- readRDS(file.path(DATA_DIR, "kirc_rnaseq_SE.rds"))
log_msg("  RNA-seq raw: %d genes x %d aliquots", nrow(rna_se), ncol(rna_se))
rna_cd <- as.data.frame(colData(rna_se))[, c("barcode", "patient", "sample_type")]
keep <- rna_cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
rna_se <- rna_se[, keep]
rna_cd <- rna_cd[keep, , drop = FALSE]
log_msg("  RNA-seq PT+STN: %d aliquots", ncol(rna_se))

gene_names <- as.character(rowData(rna_se)$gene_name)
tpm <- assay(rna_se, "tpm_unstrand")
rna_keys <- paste(rna_cd$patient, rna_cd$sample_type, sep = "|")
unique_keys <- unique(rna_keys)
# Average TPM across aliquots within (patient, sample_type)
tpm_avg <- vapply(unique_keys, function(k) {
  ix <- which(rna_keys == k)
  if (length(ix) == 1) tpm[, ix] else rowMeans(tpm[, ix, drop = FALSE])
}, numeric(nrow(tpm)))
colnames(tpm_avg) <- unique_keys
# Aggregate Ensembl IDs to HGNC gene symbols (sum TPMs across transcripts of same gene)
keep_g <- !is.na(gene_names) & gene_names != ""
tpm_avg <- tpm_avg[keep_g, , drop = FALSE]
gene_lab <- gene_names[keep_g]
tpm_mat <- rowsum(tpm_avg, group = gene_lab, reorder = FALSE)
log2tpm <- log2(tpm_mat + 1)
log_msg("  log2TPM matrix: %d HGNC genes x %d (patient, sample_type) keys",
        nrow(log2tpm), ncol(log2tpm))

# === Step 3: Match samples ====================================================
log_msg("[3/8] Matching RNA-seq to methylation by (patient, sample_type)...")
meth_keys_all <- paste(meta$patient, meta$sample_type, sep = "|")
n_dup <- sum(duplicated(meth_keys_all))
if (n_dup > 0) {
  log_msg("  Aggregating %d duplicate methylation aliquots: averaging beta within (patient, sample_type)",
          n_dup)
}
beta_full <- M2beta(M)
rm(M); gc(verbose = FALSE)

# Aggregate beta_full columns by meth_keys_all (mean across replicate aliquots)
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
# One metadata row per unique key (take first occurrence)
meta_unique <- meta[!duplicated(meth_keys_all), ]
rownames(meta_unique) <- paste(meta_unique$patient, meta_unique$sample_type, sep = "|")

matched_keys <- intersect(unique_meth_keys, colnames(log2tpm))
log_msg("  Methylation unique (patient, sample_type) keys: %d", length(unique_meth_keys))
log_msg("  Matched keys (both modalities): %d", length(matched_keys))

match_meta <- meta_unique[matched_keys, , drop = FALSE]
beta_mat <- beta_meth[, matched_keys]
log2tpm <- log2tpm[, matched_keys]
log_msg("  Sample type breakdown (matched):")
print(table(match_meta$sample_type))

# === Step 4: Build promoter (CpG, gene) pairs =================================
log_msg("[4/8] Building promoter (CpG, gene) pairs...")
PROMOTER_GROUPS <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
ann_dt <- data.table(
  cpg   = rownames(ann),
  names = as.character(ann$UCSC_RefGene_Name),
  groups = as.character(ann$UCSC_RefGene_Group)
)[names != "" & !is.na(names)]

# Split parallel arrays
ann_dt[, names_list  := strsplit(names,  ";", fixed = TRUE)]
ann_dt[, groups_list := strsplit(groups, ";", fixed = TRUE)]
ann_dt[, ok_len := lengths(names_list) == lengths(groups_list)]
ann_dt <- ann_dt[ok_len == TRUE]

# Expand
expanded <- ann_dt[, .(gene = unlist(names_list), group = unlist(groups_list)),
                   by = cpg]
# Promoter only
expanded <- expanded[group %in% PROMOTER_GROUPS]
# Dedup (cpg, gene): a CpG can be annotated to same gene under multiple groups
# (e.g., TSS200 & 5'UTR); keep the most-promoter-proximal one
PROMOTER_RANK <- c("TSS200" = 1, "TSS1500" = 2, "5'UTR" = 3, "1stExon" = 4)
expanded[, group_rank := PROMOTER_RANK[group]]
setorder(expanded, cpg, gene, group_rank)
expanded <- expanded[!duplicated(expanded[, .(cpg, gene)])]
log_msg("  Raw promoter (CpG, gene) pairs: %d", nrow(expanded))

# Filter to genes in RNA-seq and CpGs in beta_mat
expanded <- expanded[gene %in% rownames(log2tpm) & cpg %in% rownames(beta_mat)]
log_msg("  After RNA-seq + beta intersection: %d pairs (%d CpGs, %d genes)",
        nrow(expanded), uniqueN(expanded$cpg), uniqueN(expanded$gene))

# === Step 5: Covariate design =================================================
log_msg("[5/8] Building covariate design matrix...")
match_meta$age_years <- match_meta$age_at_diagnosis / 365.25
match_meta$age_years[is.na(match_meta$age_years)] <- median(match_meta$age_years, na.rm = TRUE)
match_meta$age_c <- match_meta$age_years - mean(match_meta$age_years)
match_meta$gender[is.na(match_meta$gender) | match_meta$gender == ""] <- "unknown"
# Collapse small TSS levels
tss_table <- table(match_meta$tss)
small_tss <- names(tss_table)[tss_table < 5]
match_meta$tss_grp <- ifelse(match_meta$tss %in% small_tss, "OTHER", match_meta$tss)
log_msg("  Collapsed %d small-TSS levels (<5 samples each) into 'OTHER'", length(small_tss))

X_cov <- model.matrix(~ sample_type + age_c + gender + tss_grp, data = match_meta)
qrX <- qr(X_cov)
log_msg("  Design matrix: %d x %d, rank %d", nrow(X_cov), ncol(X_cov), qrX$rank)
n  <- nrow(X_cov)
p_cov <- qrX$rank
df_resid <- n - p_cov - 1   # extra -1 for the beta coefficient
Q <- qr.Q(qrX)
log_msg("  df_resid = %d", df_resid)

# Residualizer for matrix Y where columns are samples (rows are CpGs / genes)
# y_resid_i = y_i - Q (Q' y_i)  for each row vector
residualize_rowwise <- function(Y) {
  # Y: m x n  with samples as columns
  Y - tcrossprod(Y %*% Q, Q)   # = Y - Y Q Q'
}

# === Step 6: Residualize ======================================================
log_msg("[6/8] Residualizing beta and log2TPM rows...")
unique_cpgs  <- unique(expanded$cpg)
unique_genes <- unique(expanded$gene)
beta_sub <- beta_mat[unique_cpgs, ]   # m_c x n
tpm_sub  <- log2tpm[unique_genes, ]    # m_g x n

beta_resid <- residualize_rowwise(beta_sub)
tpm_resid  <- residualize_rowwise(tpm_sub)
ss_beta <- rowSums(beta_resid^2)
ss_tpm  <- rowSums(tpm_resid^2)
log_msg("  Residualized: %d CpGs, %d genes", length(unique_cpgs), length(unique_genes))

# === Step 7: Per-pair partial correlation =====================================
log_msg("[7/8] Computing per-pair partial correlations (%d pairs)...", nrow(expanded))
cpg_ix  <- match(expanded$cpg,  unique_cpgs)
gene_ix <- match(expanded$gene, unique_genes)

t0 <- Sys.time()
chunk <- 50000
np <- nrow(expanded)
r_vec <- numeric(np)
for (start in seq(1, np, by = chunk)) {
  end <- min(start + chunk - 1, np)
  ic <- cpg_ix[start:end]; ig <- gene_ix[start:end]
  num <- rowSums(beta_resid[ic, , drop = FALSE] * tpm_resid[ig, , drop = FALSE])
  den <- sqrt(ss_beta[ic] * ss_tpm[ig])
  r_vec[start:end] <- ifelse(den > 0, num / den, NA_real_)
}
log_msg("  Done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))

# r -> t -> p; beta_coef = r * sd(y_resid)/sd(x_resid)
r_clip <- pmin(pmax(r_vec, -0.9999), 0.9999)
t_vec <- r_clip * sqrt(df_resid) / sqrt(1 - r_clip^2)
p_vec <- 2 * pt(-abs(t_vec), df = df_resid)
sd_beta <- sqrt(ss_beta / df_resid)
sd_tpm  <- sqrt(ss_tpm  / df_resid)
beta_coef_vec <- r_vec * sd_tpm[gene_ix] / sd_beta[cpg_ix]

# === Step 8: pct_responder, Q1 join, silencing flag ===========================
log_msg("[8/8] Computing pct_responder, joining Q1 DM, defining silencing flag...")
is_tumor  <- match_meta$sample_type == "Primary Tumor"
is_normal <- match_meta$sample_type == "Solid Tissue Normal"
n_tumor <- sum(is_tumor)
beta_norm_mean <- rowMeans(beta_sub[, is_normal, drop = FALSE], na.rm = TRUE)
beta_tumor_only <- beta_sub[, is_tumor, drop = FALSE]
delta <- abs(beta_tumor_only - beta_norm_mean)
pct_resp_per_cpg <- rowSums(delta > 0.2, na.rm = TRUE) / n_tumor

# Q1 toptable
q1 <- fread(file.path(Q1_DIR, "toptable_q1.csv"))
q1_sub <- q1[, .(cpg, q1_logFC = logFC, q1_FDR = adj.P.Val)]

# Output table
out <- data.table(
  cpg   = expanded$cpg,
  gene  = expanded$gene,
  promoter_group = expanded$group,
  partial_r  = r_vec,
  beta_coef  = beta_coef_vec,
  assoc_p    = p_vec,
  pct_responder = pct_resp_per_cpg[cpg_ix]
)
out[, assoc_FDR := p.adjust(assoc_p, method = "BH")]
out <- merge(out, q1_sub, by = "cpg", all.x = TRUE)

# Silencing (strict): hypermethylation in tumor + within-cohort high β -> low mRNA
# This is the classical TSG epigenetic-second-hit pattern (Knudson 1971; Jones & Baylin 2007).
out[, silencing := !is.na(q1_FDR) & !is.na(assoc_FDR) &
                    q1_FDR < 0.05 & assoc_FDR < 0.05 &
                    q1_logFC > 0 & beta_coef < 0]
# Activation (secondary): hypomethylation in tumor + high β -> high mRNA
# This is the cancer-testis-antigen / oncogene de-silencing pattern; reported but not central.
out[, activation := !is.na(q1_FDR) & !is.na(assoc_FDR) &
                    q1_FDR < 0.05 & assoc_FDR < 0.05 &
                    q1_logFC < 0 & beta_coef > 0]
out[, pattern := NA_character_]
out[silencing == TRUE & pct_responder >= 0.3, pattern := "cohort_wide"]
out[silencing == TRUE & pct_responder <  0.3, pattern := "subgroup"]

setorder(out, assoc_p)
fwrite(out, file.path(OUT_DIR, "silencing_table.csv"))

# === Summary ==================================================================
n_pairs    <- nrow(out)
n_q1       <- sum(out$q1_FDR < 0.05, na.rm = TRUE)
n_assoc    <- sum(out$assoc_FDR < 0.05, na.rm = TRUE)
n_sil      <- sum(out$silencing)
n_act      <- sum(out$activation)
n_cw       <- sum(out$pattern == "cohort_wide", na.rm = TRUE)
n_sg       <- sum(out$pattern == "subgroup",    na.rm = TRUE)
sil_cpgs   <- uniqueN(out$cpg[out$silencing])
sil_genes  <- uniqueN(out$gene[out$silencing])

cat("\n========== SUMMARY ==========\n")
cat(sprintf("(CpG, gene) pairs tested:               %d\n", n_pairs))
cat(sprintf("  Q1 paired-DM FDR<0.05:                %d (%.1f%%)\n", n_q1, 100*n_q1/n_pairs))
cat(sprintf("  beta-mRNA assoc FDR<0.05:             %d (%.1f%%)\n", n_assoc, 100*n_assoc/n_pairs))
cat(sprintf("  SILENCING (DM ∩ assoc ∩ hyper->down): %d (%.2f%%)  [PRIMARY]\n",
            n_sil, 100*n_sil/n_pairs))
cat(sprintf("    Unique CpGs:    %d\n", sil_cpgs))
cat(sprintf("    Unique genes:   %d\n", sil_genes))
cat(sprintf("    cohort_wide (pct_responder >= 0.3): %d\n", n_cw))
cat(sprintf("    subgroup    (pct_responder <  0.3): %d\n", n_sg))
cat(sprintf("  ACTIVATION (hypo->up; secondary):     %d (%.2f%%)\n",
            n_act, 100*n_act/n_pairs))

cat("\n--- VHL promoter CpGs ---\n")
vhl <- out[gene == "VHL"][order(assoc_p)]
print(vhl[, .(cpg, group = promoter_group, partial_r, beta_coef,
              assoc_p, assoc_FDR, q1_logFC, q1_FDR, pct_responder, silencing, pattern)])

cat("\n--- Top 30 silencing genes (best assoc p per gene) ---\n")
top_g <- out[silencing == TRUE,
  .(best_p = min(assoc_p), n_cpgs = .N,
    median_partial_r = round(median(partial_r), 3),
    median_pct_responder = round(median(pct_responder), 3)),
  by = gene][order(best_p)][1:30]
print(top_g)

cat("\n--- Pre-registered KIRC TSGs (sanity check) ---\n")
TSG <- c("VHL","RASSF1A","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3",
         "GATA5","BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
         "TIMP3","NEFH","UCHL1","PITX2","SPG20")
tsg_summary <- out[gene %in% TSG, .(
  best_assoc_p = min(assoc_p, na.rm = TRUE),
  n_promoter_cpgs = .N,
  n_silencing = sum(silencing),
  n_activation = sum(activation),
  any_silencing = any(silencing),
  any_cohort_wide = any(pattern == "cohort_wide", na.rm = TRUE),
  any_subgroup = any(pattern == "subgroup", na.rm = TRUE)
), by = gene]
setorder(tsg_summary, best_assoc_p)
print(tsg_summary)
n_recovered <- sum(tsg_summary$any_silencing)
cat(sprintf("\nLiterature TSG silencing recovery: %d / %d (%.0f%%)\n",
            n_recovered, nrow(tsg_summary), 100*n_recovered/nrow(tsg_summary)))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
log_msg("Done. Outputs in %s", OUT_DIR)
