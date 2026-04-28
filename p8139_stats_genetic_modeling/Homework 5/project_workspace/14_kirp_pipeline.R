#!/usr/bin/env Rscript
# 28_kirp_pipeline.R
# Run the same DM + MEA + silencing + pathway pipeline on TCGA-KIRP that
# we ran on TCGA-KIRC (scripts 04, 05, 19, 20). Compare results.
#
# This is the external-replication test: thresholds, filters, and gene-set
# libraries are IDENTICAL to the KIRC run. The only thing changed is the
# cohort.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(limma)
  library(DMRcate)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(fgsea)
  library(msigdbr)
})

set.seed(8139)
WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DATA_DIR <- file.path(WORKDIR, "data_kirp")
OUT_DIR <- file.path(WORKDIR, "kirp_out")
SIL_DIR_KIRC <- file.path(WORKDIR, "silencing_out")
dir.create(OUT_DIR, showWarnings = FALSE)

t_start <- Sys.time()

# ============================================================================
# 1. Load KIRP methylation, build β-matrix, apply same preprocessing as KIRC
# ============================================================================
cat("[1/8] Loading KIRP 450K methylation SE...\n")
meth_se <- readRDS(file.path(DATA_DIR, "kirp_450k_SE.rds"))
cat(sprintf("  initial: %d CpGs × %d aliquots\n", nrow(meth_se), ncol(meth_se)))

beta_mat <- assay(meth_se)
cd <- as.data.frame(colData(meth_se))
cat("  sample types:\n"); print(table(cd$sample_type))

# Annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- rownames(beta_mat) %in% rownames(ann)
beta_mat <- beta_mat[keep, , drop = FALSE]
cat(sprintf("  after annotation intersect: %d CpGs\n", nrow(beta_mat)))

# Drop chr X / chr Y
ann_sub <- ann[rownames(beta_mat), ]
beta_mat <- beta_mat[!ann_sub$chr %in% c("chrX", "chrY"), , drop = FALSE]
cat(sprintf("  after drop chrX/Y: %d CpGs\n", nrow(beta_mat)))

# Cross-reactive + SNP-overlap removal via DMRcate
beta_mat <- rmSNPandCH(beta_mat, dist = 2, mafcut = 0.05,
                       and = TRUE, rmcrosshyb = TRUE, rmXY = FALSE)
cat(sprintf("  after DMRcate filter: %d CpGs\n", nrow(beta_mat)))

# Drop probes with > 20% missingness
mr <- rowMeans(is.na(beta_mat))
beta_mat <- beta_mat[mr <= 0.2, , drop = FALSE]
cat(sprintf("  after >20%% missing drop: %d CpGs\n", nrow(beta_mat)))

# Complete-case across all samples
ok <- complete.cases(beta_mat)
beta_mat <- beta_mat[ok, , drop = FALSE]
cat(sprintf("  after complete-case: %d CpGs\n", nrow(beta_mat)))

# Clip and convert to M
beta_mat[beta_mat < 0.001] <- 0.001
beta_mat[beta_mat > 0.999] <- 0.999
M_mat <- log2(beta_mat / (1 - beta_mat))

# Build sample metadata
meta <- data.frame(
  barcode = colnames(beta_mat),
  patient = substr(colnames(beta_mat), 1, 12),
  sample_type = cd$sample_type[match(colnames(beta_mat), cd$barcode)],
  stringsAsFactors = FALSE
)
cat("  meta sample types after filter:\n"); print(table(meta$sample_type))

saveRDS(list(beta = beta_mat, M = M_mat, meta = meta),
        file.path(OUT_DIR, "preprocess_kirp.rds"))

# ============================================================================
# 2. Identify paired patients, run DM (paired moderated t)
# ============================================================================
cat("\n[2/8] DM test on paired KIRP patients...\n")
tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
cat(sprintf("  paired patients: %d\n", length(paired_pat)))

pick_bc <- function(p, type) {
  meta$barcode[meta$patient == p & meta$sample_type == type][1]
}
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")
M_diff <- M_mat[, bc_T] - M_mat[, bc_N]
cat(sprintf("  M_diff: %d CpGs × %d patients\n", nrow(M_diff), ncol(M_diff)))

fit <- lmFit(M_diff, design = matrix(1, ncol(M_diff), 1, dimnames = list(NULL, "x")))
fit <- eBayes(fit)
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
tt$cpg <- rownames(tt)
ann_kept <- ann[rownames(tt), ]
tt$chr <- ann_kept$chr; tt$pos <- ann_kept$pos
tt$gene <- ann_kept$UCSC_RefGene_Name; tt$feat <- ann_kept$UCSC_RefGene_Group

n_sig <- sum(tt$adj.P.Val < 0.05)
n_hyper <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
n_hypo  <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)
cat(sprintf("  FDR < 0.05: %d (%.1f%%); hyper %d, hypo %d\n",
            n_sig, 100 * n_sig / nrow(tt), n_hyper, n_hypo))
fwrite(tt, file.path(OUT_DIR, "toptable_kirp_dm.csv"))

# ============================================================================
# 3. Build (CpG, gene) promoter pairs
# ============================================================================
cat("\n[3/8] Building promoter (CpG, gene) pairs...\n")
prom_groups <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
expand_pairs <- function(cpg, name_str, feat_str) {
  if (is.na(name_str) || name_str == "") return(NULL)
  names <- strsplit(name_str, ";")[[1]]
  feats <- strsplit(feat_str, ";")[[1]]
  if (length(names) != length(feats)) return(NULL)
  ok <- feats %in% prom_groups
  if (!any(ok)) return(NULL)
  data.table(cpg = cpg, gene = names[ok], feat = feats[ok])
}
pairs_list <- vector("list", nrow(tt))
for (i in seq_len(nrow(tt))) {
  pairs_list[[i]] <- expand_pairs(tt$cpg[i], tt$gene[i], tt$feat[i])
}
pairs <- rbindlist(pairs_list)
# Deduplicate by (cpg, gene), keeping most-promoter-proximal feat
group_priority <- c("TSS200" = 1, "1stExon" = 2, "5'UTR" = 3, "TSS1500" = 4)
pairs[, prio := group_priority[feat]]
pairs <- pairs[order(prio)][!duplicated(pairs[, .(cpg, gene)])]
pairs[, prio := NULL]
cat(sprintf("  unique (CpG, gene) promoter pairs: %d (%d unique CpGs, %d unique genes)\n",
            nrow(pairs), uniqueN(pairs$cpg), uniqueN(pairs$gene)))

# ============================================================================
# 4. Load KIRP RNA-seq, build TPM matrix, match samples
# ============================================================================
cat("\n[4/8] Loading KIRP RNA-seq, matching samples...\n")
rna_se <- readRDS(file.path(DATA_DIR, "kirp_rnaseq_SE.rds"))
rna_cd <- as.data.frame(colData(rna_se))[, c("barcode","patient","sample_type")]
keep <- rna_cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
rna_se <- rna_se[, keep]; rna_cd <- rna_cd[keep, , drop = FALSE]
gene_names <- as.character(rowData(rna_se)$gene_name)
tpm <- assay(rna_se, "tpm_unstrand")

# Aggregate replicate aliquots
rna_keys <- paste(rna_cd$patient, rna_cd$sample_type, sep = "|")
unique_rna_keys <- unique(rna_keys)
tpm_avg <- vapply(unique_rna_keys, function(k) {
  ix <- which(rna_keys == k)
  if (length(ix) == 1) tpm[, ix] else rowMeans(tpm[, ix, drop = FALSE])
}, numeric(nrow(tpm)))
colnames(tpm_avg) <- unique_rna_keys
keep_g <- !is.na(gene_names) & gene_names != ""
tpm_avg <- tpm_avg[keep_g, , drop = FALSE]
gene_lab <- gene_names[keep_g]
log2tpm <- log2(rowsum(tpm_avg, group = gene_lab, reorder = FALSE) + 1)

# Aggregate methylation aliquots same way
meth_keys <- paste(meta$patient, meta$sample_type, sep = "|")
unique_meth_keys <- unique(meth_keys)
beta_agg <- vapply(unique_meth_keys, function(k) {
  ix <- which(meth_keys == k)
  if (length(ix) == 1) beta_mat[, ix] else rowMeans(beta_mat[, ix, drop = FALSE], na.rm = TRUE)
}, numeric(nrow(beta_mat)))
colnames(beta_agg) <- unique_meth_keys
meta_unique <- meta[!duplicated(meth_keys), ]
rownames(meta_unique) <- paste(meta_unique$patient, meta_unique$sample_type, sep = "|")

matched_keys <- intersect(unique_meth_keys, colnames(log2tpm))
beta_agg <- beta_agg[, matched_keys, drop = FALSE]
log2tpm  <- log2tpm[, matched_keys, drop = FALSE]
match_meta <- meta_unique[matched_keys, , drop = FALSE]
n_t <- sum(match_meta$sample_type == "Primary Tumor")
n_n <- sum(match_meta$sample_type == "Solid Tissue Normal")
cat(sprintf("  matched (patient, sample-type) keys: %d (%d tumor + %d normal)\n",
            length(matched_keys), n_t, n_n))

# Restrict pairs to genes available in RNA-seq
pairs <- pairs[gene %in% rownames(log2tpm) & cpg %in% rownames(beta_agg)]
cat(sprintf("  pairs after intersection with RNA-seq: %d\n", nrow(pairs)))

# ============================================================================
# 5. MEA test (partial correlation, sample_type-adjusted)
# ============================================================================
cat("\n[5/8] MEA test (vectorized partial correlation)...\n")
# Build covariate design: sample_type only (KIRP clinical may lack age/sex/TSS)
# To keep replication faithful we'll add available covariates from clinical
clinical <- readRDS(file.path(DATA_DIR, "clinical.rds"))
clin_meta <- clinical[match(match_meta$patient, clinical$submitter_id),
                      c("submitter_id", "age_at_index", "gender",
                        "tissue_or_organ_of_origin")]
match_meta$age <- as.numeric(clin_meta$age_at_index)
match_meta$sex <- ifelse(clin_meta$gender == "male", 1, 0)
match_meta$tss <- substr(match_meta$patient, 6, 7)  # TSS code

# Build design (intercept + sample_type + age + sex + TSS dummies)
sample_type_num <- as.numeric(match_meta$sample_type == "Primary Tumor")
# TSS as factor — drop reference level
tss_factor <- factor(match_meta$tss)
tss_design <- model.matrix(~ tss_factor)[, -1, drop = FALSE]
# age/sex with NA → impute median for design (only as covariate; preserves df)
match_meta$age[is.na(match_meta$age)] <- median(match_meta$age, na.rm = TRUE)
match_meta$sex[is.na(match_meta$sex)] <- 0

X <- cbind(1,
           sample_type = sample_type_num,
           age = match_meta$age,
           sex = match_meta$sex,
           tss_design)
df_resid <- ncol(beta_agg) - ncol(X) - 1  # one extra for the β term being tested

# QR residualization
qrX <- qr(X)
resid_y <- function(Y) Y - qrX$qr %*% (qr.coef(qrX, t(Y)))  # not used; use lm.fit

# Use limma::lmFit + topTable approach? No: we want one-coefficient test of β_CpG.
# Faster: residualize log2TPM and β against X; then partial correlation.
cat(sprintf("  design df = %d resid; running %d pair tests...\n",
            df_resid, nrow(pairs)))
# Vectorized residualization
resid_mat <- function(Y, X) Y - X %*% qr.solve(crossprod(X), crossprod(X, Y))
# Y = log2tpm matrix (rows: genes, cols: samples) — residualize columns
log2tpm_resid <- t(resid_mat(t(log2tpm), X))
beta_resid    <- t(resid_mat(t(beta_agg), X))

# Per-pair test: partial correlation between residualized log2tpm[gene] and residualized beta_agg[cpg]
n_samples <- ncol(beta_agg)
test_one <- function(cpg_id, gene_id) {
  bv <- beta_resid[cpg_id, ]
  yv <- log2tpm_resid[gene_id, ]
  ok <- is.finite(bv) & is.finite(yv)
  if (sum(ok) < 10) return(c(NA_real_, NA_real_, NA_real_))
  bv <- bv[ok]; yv <- yv[ok]
  r <- cor(bv, yv)
  n <- length(bv)
  beta_coef <- r * sd(yv) / sd(bv)
  t_stat <- r * sqrt((n - 2 - ncol(X)) / (1 - r^2 + 1e-12))
  p <- 2 * pt(abs(t_stat), df = n - 2 - ncol(X), lower.tail = FALSE)
  c(r, beta_coef, p)
}

# Run in batches with mapply
cat("  computing per-pair statistics...\n")
mea_res <- mapply(test_one, pairs$cpg, pairs$gene)
mea_dt <- data.table(
  cpg = pairs$cpg,
  gene = pairs$gene,
  feat = pairs$feat,
  partial_r = mea_res[1, ],
  beta_coef = mea_res[2, ],
  assoc_p   = mea_res[3, ]
)
mea_dt[, assoc_FDR := p.adjust(assoc_p, method = "BH")]
cat(sprintf("  MEA done. FDR < 0.05: %d / %d (%.1f%%)\n",
            sum(mea_dt$assoc_FDR < 0.05, na.rm = TRUE),
            nrow(mea_dt),
            100 * sum(mea_dt$assoc_FDR < 0.05, na.rm = TRUE) / nrow(mea_dt)))

# ============================================================================
# 6. Silencing definition (DM ∩ MEA ∩ direction concordance)
# ============================================================================
cat("\n[6/8] Silencing definition...\n")
mea_dt[, q1_logFC := tt$logFC[match(cpg, tt$cpg)]]
mea_dt[, q1_FDR   := tt$adj.P.Val[match(cpg, tt$cpg)]]
mea_dt[, silencing := q1_FDR < 0.05 & assoc_FDR < 0.05 &
                     q1_logFC > 0 & beta_coef < 0]
mea_dt[, activation := q1_FDR < 0.05 & assoc_FDR < 0.05 &
                      q1_logFC < 0 & beta_coef > 0]

# Per-CpG responder fraction (using tumor β only)
tumor_keys <- match_meta$sample_type == "Primary Tumor"
normal_keys <- match_meta$sample_type == "Solid Tissue Normal"
beta_tumor <- beta_agg[, tumor_keys, drop = FALSE]
beta_normal <- beta_agg[, normal_keys, drop = FALSE]
beta_normal_mean <- rowMeans(beta_normal, na.rm = TRUE)
n_tumor_keys <- sum(tumor_keys)
pct_responder <- numeric(length(unique(mea_dt$cpg)))
names(pct_responder) <- unique(mea_dt$cpg)
unique_cpgs <- intersect(names(pct_responder), rownames(beta_tumor))
for (cpg_id in unique_cpgs) {
  bt <- beta_tumor[cpg_id, ]
  bn <- beta_normal_mean[cpg_id]
  pct_responder[cpg_id] <- mean(abs(bt - bn) > 0.2, na.rm = TRUE)
}
mea_dt[, pct_responder := pct_responder[cpg]]

n_sil <- sum(mea_dt$silencing, na.rm = TRUE)
n_act <- sum(mea_dt$activation, na.rm = TRUE)
n_sil_genes <- length(unique(mea_dt$gene[mea_dt$silencing == TRUE]))
n_act_genes <- length(unique(mea_dt$gene[mea_dt$activation == TRUE]))
cat(sprintf("  silencing pairs: %d (%d unique genes)\n", n_sil, n_sil_genes))
cat(sprintf("  activation pairs: %d (%d unique genes)\n", n_act, n_act_genes))

fwrite(mea_dt, file.path(OUT_DIR, "kirp_silencing_table.csv"))

# ============================================================================
# 7. Pathway enrichment: raw-DM GSEA + silencing ORA + immune dropout test
# ============================================================================
cat("\n[7/8] Pathway enrichment (replicating immune dropout test)...\n")
# Per-gene best DM p (signed by DM logFC)
best_dm <- mea_dt[, .(best_p = min(q1_FDR), best_logFC = q1_logFC[which.min(q1_FDR)]),
                  by = gene]
best_dm <- best_dm[!is.na(best_p)]
ranks_dm <- with(best_dm, sign(best_logFC) * (-log10(pmax(best_p, 1e-300))))
names(ranks_dm) <- best_dm$gene

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
kegg_list <- split(kegg$gene_symbol, kegg$gs_name)
all_sets <- c(hallmark_list, kegg_list)

cat("  raw-DM GSEA...\n")
gsea_dm <- fgsea(pathways = all_sets, stats = ranks_dm, minSize = 15, maxSize = 500)
gsea_dm <- gsea_dm[order(padj)]
fwrite(gsea_dm[, .(pathway, pval, padj, NES, size, ES)],
       file.path(OUT_DIR, "kirp_gsea_rawdm.csv"))

cat("  silencing ORA...\n")
sil_genes <- unique(mea_dt$gene[mea_dt$silencing == TRUE])
all_tested_genes <- unique(mea_dt$gene)
ora_results <- data.frame()
for (path_name in names(all_sets)) {
  path_genes <- all_sets[[path_name]]
  k <- length(intersect(sil_genes, path_genes))
  K <- length(intersect(all_tested_genes, path_genes))
  n <- length(sil_genes)
  N <- length(all_tested_genes)
  if (K < 5) next
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  ora_results <- rbind(ora_results, data.frame(
    pathway = path_name, k = k, K = K, n = n, N = N, p = p
  ))
}
ora_results$padj <- p.adjust(ora_results$p, method = "BH")
ora_results <- ora_results[order(ora_results$padj), ]
fwrite(ora_results, file.path(OUT_DIR, "kirp_ora_silencing.csv"))

# Immune-pathway dropout
immune_paths <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_COMPLEMENT",
                  "HALLMARK_ALLOGRAFT_REJECTION",
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                  "KEGG_HEMATOPOIETIC_CELL_LINEAGE")
imm_dt <- merge(
  gsea_dm[pathway %in% immune_paths, .(pathway, dm_padj = padj, dm_NES = NES)],
  data.table(ora_results)[pathway %in% immune_paths, .(pathway, sil_padj = padj, sil_p = p, k, K)],
  by = "pathway", all = TRUE)
imm_dt <- imm_dt[order(dm_padj)]
fwrite(imm_dt, file.path(OUT_DIR, "kirp_immune_dropout.csv"))
cat("  immune dropout summary:\n"); print(imm_dt)

# ============================================================================
# 8. Compare with KIRC silencing list
# ============================================================================
cat("\n[8/8] Cross-cohort comparison KIRC vs KIRP...\n")
kirc_sil <- fread(file.path(SIL_DIR_KIRC, "silencing_table.csv"))
kirc_genes <- unique(kirc_sil$gene[kirc_sil$silencing == TRUE])
kirp_genes <- unique(mea_dt$gene[mea_dt$silencing == TRUE])

bg <- intersect(unique(mea_dt$gene), unique(kirc_sil$gene))
kirc_in_bg <- intersect(kirc_genes, bg)
kirp_in_bg <- intersect(kirp_genes, bg)
overlap <- intersect(kirc_in_bg, kirp_in_bg)
cat(sprintf("  background (genes tested in both): %d\n", length(bg)))
cat(sprintf("  KIRC silencing in bg: %d  | KIRP silencing in bg: %d  | overlap: %d\n",
            length(kirc_in_bg), length(kirp_in_bg), length(overlap)))

if (length(kirp_in_bg) > 0) {
  hyp_p <- phyper(length(overlap) - 1, length(kirp_in_bg),
                  length(bg) - length(kirp_in_bg), length(kirc_in_bg),
                  lower.tail = FALSE)
  jacc <- length(overlap) / length(union(kirc_in_bg, kirp_in_bg))
  cat(sprintf("  Jaccard: %.3f\n", jacc))
  cat(sprintf("  hypergeometric p (overlap > random): %.2e\n", hyp_p))
} else {
  hyp_p <- NA; jacc <- NA
}

# 18-TSG recovery in KIRP
TSG18 <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
           "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
           "TIMP3","NEFH","UCHL1","PITX2")
kirp_tsg_sil <- intersect(TSG18, kirp_genes)
kirc_tsg_sil <- intersect(TSG18, kirc_genes)
both_tsg <- intersect(kirp_tsg_sil, kirc_tsg_sil)
cat(sprintf("\n  18-TSG recovery:\n"))
cat(sprintf("    KIRC: %d/18  %s\n", length(kirc_tsg_sil), paste(kirc_tsg_sil, collapse = ", ")))
cat(sprintf("    KIRP: %d/18  %s\n", length(kirp_tsg_sil), paste(kirp_tsg_sil, collapse = ", ")))
cat(sprintf("    both: %d  %s\n", length(both_tsg), paste(both_tsg, collapse = ", ")))

# Summary CSV
summary_dt <- data.frame(
  metric = c("kirp_paired_patients", "kirp_dm_n_sig",
             "kirp_silencing_pairs", "kirp_silencing_genes",
             "background_in_both", "kirc_sil_in_bg", "kirp_sil_in_bg",
             "overlap", "jaccard", "hypergeometric_p",
             "kirp_TSG_recovered", "kirc_TSG_recovered", "both_TSG_recovered",
             "immune_HALLMARK_INFLAM_kirp_dm_padj",
             "immune_HALLMARK_INFLAM_kirp_sil_padj"),
  value  = c(length(paired_pat), n_sig, n_sil, n_sil_genes,
             length(bg), length(kirc_in_bg), length(kirp_in_bg),
             length(overlap),
             round(jacc %||% NA, 4),
             formatC(hyp_p %||% NA, format = "e", digits = 2),
             length(kirp_tsg_sil), length(kirc_tsg_sil), length(both_tsg),
             tryCatch(formatC(imm_dt$dm_padj[imm_dt$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"],
                              format = "e", digits = 2), error = function(e) NA),
             tryCatch(formatC(imm_dt$sil_padj[imm_dt$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"],
                              format = "g", digits = 2), error = function(e) NA))
)
fwrite(summary_dt, file.path(OUT_DIR, "kirp_replication_summary.csv"))

t_end <- Sys.time()
cat(sprintf("\nTotal: %.1f min\n", as.numeric(difftime(t_end, t_start, units = "mins"))))
cat(sprintf("Outputs in %s\n", OUT_DIR))
