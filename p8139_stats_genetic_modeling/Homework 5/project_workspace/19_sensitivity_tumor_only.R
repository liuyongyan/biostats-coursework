#!/usr/bin/env Rscript
# 19_sensitivity_tumor_only.R
# Independence-sensitivity check addressing advisor feedback (Apr 2026):
# does the pooled MEA model — which includes 24 paired adjacent-normal samples
# alongside 318 tumors and uses sample_type as a fixed covariate — give the
# same silencing list as a tumor-only regression on the 318 independent
# patients?
#
# If the two analyses agree closely on (a) beta_MEA estimates, (b) the
# silencing list, and (c) the canonical VHL/cg13672843 signal, then the
# fixed-effects pooled formulation is empirically justified despite the 7%
# paired structure (24/342) violating strict independence. If they diverge,
# the silencing screen should be rerun on tumor-only or with a mixed-effects
# patient random intercept.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
Q1_DIR   <- file.path(WORKDIR, "q1_out")
SIL_DIR  <- file.path(WORKDIR, "silencing_out")
OUT_DIR  <- file.path(WORKDIR, "sensitivity_tumoronly")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(8139)

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")

# === Step 1-3: same data loading as 07_methylation_silencing.R =================
log_msg("[1/7] Loading methylation matrix and metadata...")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
ann  <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
M2beta <- function(m) 2^m / (2^m + 1)

log_msg("[2/7] Loading RNA-seq...")
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

log_msg("[3/7] Matching by (patient, sample_type)...")
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

# === Step 4: TUMOR-ONLY subset ================================================
log_msg("[4/7] Subsetting to tumor-only (318 independent patients)...")
match_meta_full <- meta_unique[matched_keys, , drop = FALSE]
is_tumor_full <- match_meta_full$sample_type == "Primary Tumor"
tumor_keys <- matched_keys[is_tumor_full]
match_meta <- meta_unique[tumor_keys, , drop = FALSE]
beta_mat <- beta_meth[, tumor_keys]
log2tpm_mat <- log2tpm[, tumor_keys]
log_msg("  Tumor-only sample count: %d", length(tumor_keys))
stopifnot(length(unique(match_meta$patient)) == length(tumor_keys))   # one row per patient

# Promoter pair set (same as 07)
log_msg("  Building promoter pairs...")
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
log_msg("  Promoter pairs to test: %d", nrow(expanded))

# === Step 5: Covariate design (no sample_type) =================================
log_msg("[5/7] Building covariate design (age + sex + TSS, NO sample_type)...")
match_meta$age_years <- match_meta$age_at_diagnosis / 365.25
match_meta$age_years[is.na(match_meta$age_years)] <- median(match_meta$age_years, na.rm = TRUE)
match_meta$age_c <- match_meta$age_years - mean(match_meta$age_years)
match_meta$gender[is.na(match_meta$gender) | match_meta$gender == ""] <- "unknown"
tss_table <- table(match_meta$tss)
small_tss <- names(tss_table)[tss_table < 5]
match_meta$tss_grp <- ifelse(match_meta$tss %in% small_tss, "OTHER", match_meta$tss)

X_cov <- model.matrix(~ age_c + gender + tss_grp, data = match_meta)
qrX <- qr(X_cov)
n  <- nrow(X_cov)
p_cov <- qrX$rank
df_resid <- n - p_cov - 1
Q <- qr.Q(qrX)
log_msg("  Design matrix: %d x %d, rank %d; df_resid = %d",
        nrow(X_cov), ncol(X_cov), qrX$rank, df_resid)

residualize_rowwise <- function(Y) Y - tcrossprod(Y %*% Q, Q)

# === Step 6: Per-pair partial correlation =====================================
log_msg("[6/7] Residualizing and computing %d partial correlations...", nrow(expanded))
unique_cpgs  <- unique(expanded$cpg)
unique_genes <- unique(expanded$gene)
beta_sub <- beta_mat[unique_cpgs, ]
tpm_sub  <- log2tpm_mat[unique_genes, ]
beta_resid <- residualize_rowwise(beta_sub)
tpm_resid  <- residualize_rowwise(tpm_sub)
ss_beta <- rowSums(beta_resid^2)
ss_tpm  <- rowSums(tpm_resid^2)

cpg_ix  <- match(expanded$cpg,  unique_cpgs)
gene_ix <- match(expanded$gene, unique_genes)

t0 <- Sys.time()
chunk <- 50000; np <- nrow(expanded)
r_vec <- numeric(np)
for (start in seq(1, np, by = chunk)) {
  end <- min(start + chunk - 1, np)
  ic <- cpg_ix[start:end]; ig <- gene_ix[start:end]
  num <- rowSums(beta_resid[ic, , drop = FALSE] * tpm_resid[ig, , drop = FALSE])
  den <- sqrt(ss_beta[ic] * ss_tpm[ig])
  r_vec[start:end] <- ifelse(den > 0, num / den, NA_real_)
}
log_msg("  Done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))

r_clip <- pmin(pmax(r_vec, -0.9999), 0.9999)
t_vec <- r_clip * sqrt(df_resid) / sqrt(1 - r_clip^2)
p_vec <- 2 * pt(-abs(t_vec), df = df_resid)
sd_beta <- sqrt(ss_beta / df_resid)
sd_tpm  <- sqrt(ss_tpm  / df_resid)
beta_coef_vec <- r_vec * sd_tpm[gene_ix] / sd_beta[cpg_ix]

# Q1 paired-DM toptable (unchanged — paired test on 160 patients, independent of MEA)
q1 <- fread(file.path(Q1_DIR, "toptable_q1.csv"))
q1_sub <- q1[, .(cpg, q1_logFC = logFC, q1_FDR = adj.P.Val)]

out <- data.table(
  cpg   = expanded$cpg,
  gene  = expanded$gene,
  promoter_group = expanded$group,
  partial_r_to  = r_vec,
  beta_coef_to  = beta_coef_vec,
  assoc_p_to    = p_vec
)
out[, assoc_FDR_to := p.adjust(assoc_p_to, method = "BH")]
out <- merge(out, q1_sub, by = "cpg", all.x = TRUE)
out[, silencing_to := !is.na(q1_FDR) & !is.na(assoc_FDR_to) &
                      q1_FDR < 0.05 & assoc_FDR_to < 0.05 &
                      q1_logFC > 0 & beta_coef_to < 0]
out[, activation_to := !is.na(q1_FDR) & !is.na(assoc_FDR_to) &
                       q1_FDR < 0.05 & assoc_FDR_to < 0.05 &
                       q1_logFC < 0 & beta_coef_to > 0]
fwrite(out, file.path(OUT_DIR, "silencing_table_tumoronly.csv"))

# === Step 7: Compare against pooled fit =======================================
log_msg("[7/7] Comparing tumor-only vs pooled silencing screen...")
pooled <- fread(file.path(SIL_DIR, "silencing_table.csv"))
key <- c("cpg", "gene")
m <- merge(
  pooled[, .(cpg, gene, beta_coef_p = beta_coef, assoc_p_p = assoc_p,
             assoc_FDR_p = assoc_FDR, silencing_p = silencing,
             activation_p = activation, q1_logFC, q1_FDR,
             pct_responder)],
  out[, .(cpg, gene, beta_coef_to, assoc_p_to, assoc_FDR_to,
          silencing_to, activation_to)],
  by = key, all = FALSE
)

valid <- m[!is.na(beta_coef_p) & !is.na(beta_coef_to)]

# (a) beta_MEA agreement
spearman_beta <- cor(valid$beta_coef_p, valid$beta_coef_to, method = "spearman")
pearson_beta  <- cor(valid$beta_coef_p, valid$beta_coef_to, method = "pearson")
sign_concord  <- mean(sign(valid$beta_coef_p) == sign(valid$beta_coef_to), na.rm = TRUE)

# (b) silencing list overlap
sil_p  <- m[silencing_p  == TRUE,  paste(cpg, gene, sep = "|")]
sil_to <- m[silencing_to == TRUE,  paste(cpg, gene, sep = "|")]
overlap <- length(intersect(sil_p, sil_to))
union_n <- length(union(sil_p, sil_to))
jaccard_sil <- overlap / union_n

# Activation list overlap
act_p  <- m[activation_p  == TRUE, paste(cpg, gene, sep = "|")]
act_to <- m[activation_to == TRUE, paste(cpg, gene, sep = "|")]
jaccard_act <- length(intersect(act_p, act_to)) / length(union(act_p, act_to))

# (c) VHL/cg13672843
vhl <- m[cpg == "cg13672843" & gene == "VHL"]

# (d) TSG recovery
TSG <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
         "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
         "TIMP3","NEFH","UCHL1","PITX2")
tsg_recovery <- m[gene %in% TSG, .(
  any_sil_pooled = any(silencing_p),
  any_sil_tumoronly = any(silencing_to)
), by = gene]
tsg_recovery[, agree := any_sil_pooled == any_sil_tumoronly]

# Discordant pairs: silencing in one but not the other
discord_pool_only <- m[silencing_p == TRUE & silencing_to == FALSE]
discord_to_only   <- m[silencing_p == FALSE & silencing_to == TRUE]

# === Save and report ==========================================================
summary_lines <- c(
  "================================================================",
  " SENSITIVITY: tumor-only (318 indep) vs pooled (318+24, sample_type fixed effect)",
  "================================================================",
  sprintf(" Pairs evaluated (both fits non-NA):    %d",        nrow(valid)),
  "",
  " (a) beta_MEA agreement",
  sprintf("     Spearman rho:                       %.4f",     spearman_beta),
  sprintf("     Pearson  rho:                       %.4f",     pearson_beta),
  sprintf("     Sign concordance:                   %.1f%%",   100*sign_concord),
  "",
  " (b) Silencing list (q1 FDR<.05 & MEA FDR<.05 & directional)",
  sprintf("     Pooled silencing pairs:             %d",       sum(m$silencing_p,  na.rm=TRUE)),
  sprintf("     Tumor-only silencing pairs:         %d",       sum(m$silencing_to, na.rm=TRUE)),
  sprintf("     Intersection:                       %d",       overlap),
  sprintf("     Jaccard:                            %.3f",     jaccard_sil),
  sprintf("     Activation Jaccard:                 %.3f",     jaccard_act),
  sprintf("     Pooled-only (lost in tumor-only):   %d",       nrow(discord_pool_only)),
  sprintf("     Tumor-only-only (gained):           %d",       nrow(discord_to_only)),
  "",
  " (c) VHL / cg13672843 (positive control)",
  sprintf("     pooled  : beta=%+.3f, p=%.2e, FDR=%.2e, sil=%s",
          vhl$beta_coef_p, vhl$assoc_p_p, vhl$assoc_FDR_p, vhl$silencing_p),
  sprintf("     tumor-only: beta=%+.3f, p=%.2e, FDR=%.2e, sil=%s",
          vhl$beta_coef_to, vhl$assoc_p_to, vhl$assoc_FDR_to, vhl$silencing_to),
  "",
  " (d) Predefined TSG silencing recovery (any-CpG level)",
  sprintf("     Pooled recovers:    %d / %d", sum(tsg_recovery$any_sil_pooled),    nrow(tsg_recovery)),
  sprintf("     Tumor-only recovers:%d / %d", sum(tsg_recovery$any_sil_tumoronly), nrow(tsg_recovery)),
  sprintf("     Genes with same call:%d / %d (%.0f%%)",
          sum(tsg_recovery$agree), nrow(tsg_recovery),
          100*mean(tsg_recovery$agree)),
  "================================================================"
)
writeLines(summary_lines, file.path(OUT_DIR, "summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n--- TSG-level disagreements ---\n")
print(tsg_recovery[agree == FALSE])

fwrite(tsg_recovery, file.path(OUT_DIR, "tsg_recovery_compare.csv"))
fwrite(m, file.path(OUT_DIR, "pair_level_compare.csv.gz"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
log_msg("Done. Outputs in %s", OUT_DIR)
