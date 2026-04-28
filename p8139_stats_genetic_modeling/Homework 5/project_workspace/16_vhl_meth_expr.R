# C2 analysis: VHL methylation × expression at cg13672843
# Test whether the Layer-2 hit cg13672843 is mechanistically silencing VHL:
# Expectation: in tumor samples, high methylation at cg13672843 → low VHL mRNA
# If true, this converts the Layer-2 narrative from speculation to evidence.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
OUT_DIR  <- file.path(WORKDIR, "vhl_expr_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Step 1: Download TCGA-KIRC RNA-seq STAR Counts (if not cached) ----------
rna_rds <- file.path(DATA_DIR, "kirc_rnaseq_SE.rds")
if (!file.exists(rna_rds)) {
  cat("Querying TCGA-KIRC RNA-seq STAR Counts ...\n")
  q <- GDCquery(
    project = "TCGA-KIRC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  cat("  ", nrow(getResults(q)), "files queried\n")
  GDCdownload(q, directory = file.path(DATA_DIR, "GDCdata"))
  cat("Preparing SE...\n")
  rna_se <- GDCprepare(q, directory = file.path(DATA_DIR, "GDCdata"))
  saveRDS(rna_se, rna_rds)
  cat("Saved", rna_rds, "\n")
} else {
  cat("Loading cached RNA-seq SE...\n")
  rna_se <- readRDS(rna_rds)
}

cat("RNA-seq SE: ", nrow(rna_se), "genes ×", ncol(rna_se), "samples\n")

# --- Step 2: Get VHL expression (TPM) -----------------------------------------
# VHL Ensembl ID: ENSG00000134086
gene_ids <- rownames(rna_se)
vhl_idx <- grep("ENSG00000134086", gene_ids)
if (length(vhl_idx) == 0) {
  # Fall back: search by gene_name in rowData
  rd <- as.data.frame(rowData(rna_se))
  vhl_idx <- which(rd$gene_name == "VHL")
}
stopifnot(length(vhl_idx) >= 1)
cat("VHL row index:", vhl_idx, " gene_id:", rownames(rna_se)[vhl_idx], "\n")

# Try to get TPM from assays (TCGAbiolinks STAR-Counts SE has tpm_unstrand)
assay_names <- assayNames(rna_se)
cat("Available assays:", paste(assay_names, collapse=", "), "\n")
if ("tpm_unstrand" %in% assay_names) {
  vhl_expr <- assay(rna_se, "tpm_unstrand")[vhl_idx, ]
  expr_lab <- "VHL TPM"
} else if ("fpkm_unstrand" %in% assay_names) {
  vhl_expr <- assay(rna_se, "fpkm_unstrand")[vhl_idx, ]
  expr_lab <- "VHL FPKM"
} else {
  # raw counts; do a simple library-size normalization
  cts <- assay(rna_se, "unstranded")
  libsize <- colSums(cts)
  vhl_expr <- cts[vhl_idx, ] / libsize * 1e6  # CPM-like
  expr_lab <- "VHL CPM"
}
cat("VHL expression: n=", length(vhl_expr), " range=[",
    formatC(min(vhl_expr), digits=2, format="f"), ",",
    formatC(max(vhl_expr), digits=2, format="f"), "]\n")

# Sample metadata for RNA-seq
rna_meta <- as.data.frame(colData(rna_se))[, c("barcode","sample_type","patient")]
rna_meta$vhl_log2 <- log2(vhl_expr + 1)
cat("RNA-seq sample types:\n"); print(table(rna_meta$sample_type))

# --- Step 3: Match RNA-seq samples to methylation samples by patient ----------
M     <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta  <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
M2beta <- function(m) 2^m / (2^m + 1)

target_cpg <- "cg13672843"  # VHL 1stExon Layer-2 hit
if (!target_cpg %in% rownames(M)) stop("Target CpG not in M matrix")
beta_target <- M2beta(M[target_cpg, ])  # named by methylation barcode
cat("\ncg13672843 β stats: range [", formatC(min(beta_target),3,format="f"), ",",
    formatC(max(beta_target),3,format="f"), "], mean=",
    formatC(mean(beta_target),3,format="f"), "\n")

# Methylation per-patient β at cg13672843, by sample type
meth_df <- data.frame(
  barcode = meta$barcode,
  patient = meta$patient,
  sample_type = meta$sample_type,
  beta = beta_target[meta$barcode]
)

# RNA-seq per-patient log2(TPM+1)
rna_df <- rna_meta %>%
  mutate(simple_type = ifelse(sample_type == "Primary Tumor", "Primary Tumor",
                       ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", NA))) %>%
  filter(!is.na(simple_type))

# Merge by patient + sample_type
joined <- merge(
  meth_df,    rna_df %>% select(patient, simple_type, vhl_log2),
  by.x = c("patient","sample_type"),  by.y = c("patient","simple_type"),
  all = FALSE
)
# Drop sample-pair duplicates (when one patient has two RNA-seq aliquots, take mean)
joined <- joined %>% group_by(patient, sample_type) %>%
  summarise(beta = mean(beta), vhl_log2 = mean(vhl_log2), .groups = "drop") %>%
  as.data.frame()

cat("\nJoined methylation × RNA-seq:\n")
print(table(joined$sample_type))

# --- Step 4: Correlations ----------------------------------------------------
tumor_df  <- joined %>% filter(sample_type == "Primary Tumor")
normal_df <- joined %>% filter(sample_type == "Solid Tissue Normal")

cor_results <- data.frame(
  contrast = c("Tumor (all)", "Normal (all)"),
  n = c(nrow(tumor_df), nrow(normal_df)),
  pearson_r  = c(cor(tumor_df$beta,  tumor_df$vhl_log2,  method="pearson"),
                 cor(normal_df$beta, normal_df$vhl_log2, method="pearson")),
  pearson_p  = c(cor.test(tumor_df$beta, tumor_df$vhl_log2,  method="pearson")$p.value,
                 cor.test(normal_df$beta, normal_df$vhl_log2, method="pearson")$p.value),
  spearman_rho = c(cor(tumor_df$beta,  tumor_df$vhl_log2,  method="spearman"),
                   cor(normal_df$beta, normal_df$vhl_log2, method="spearman")),
  spearman_p   = c(cor.test(tumor_df$beta, tumor_df$vhl_log2,  method="spearman", exact=FALSE)$p.value,
                   cor.test(normal_df$beta, normal_df$vhl_log2, method="spearman", exact=FALSE)$p.value)
)
cat("\nCorrelation results:\n")
print(cor_results)
write.csv(cor_results, file.path(OUT_DIR, "vhl_corr_results.csv"), row.names = FALSE)

# --- Step 5: Define methylation-high subgroup, compare expression ------------
# Use a tumor-specific β threshold: top quartile vs bottom quartile of cg13672843 in tumor samples
q3 <- quantile(tumor_df$beta, 0.75)
q1 <- quantile(tumor_df$beta, 0.25)
tumor_df <- tumor_df %>%
  mutate(meth_group = case_when(
    beta >= q3 ~ "Methylation-high (Q4)",
    beta <= q1 ~ "Methylation-low (Q1)",
    TRUE       ~ "Middle (Q2-Q3)"
  ))

high <- tumor_df %>% filter(meth_group == "Methylation-high (Q4)")
low  <- tumor_df %>% filter(meth_group == "Methylation-low (Q1)")
wt   <- wilcox.test(high$vhl_log2, low$vhl_log2, alternative = "less")  # high < low?
cat(sprintf("\nQ4 vs Q1 VHL expression: high mean log2(TPM+1)=%.2f (n=%d), low mean=%.2f (n=%d), Wilcoxon (high<low) p=%.4g\n",
            mean(high$vhl_log2), nrow(high),
            mean(low$vhl_log2),  nrow(low), wt$p.value))

# --- Step 6: Plots -----------------------------------------------------------
p_scatter <- ggplot(joined, aes(x = beta, y = vhl_log2, color = sample_type)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
  scale_color_manual(values = c("Primary Tumor" = "#d73027",
                                 "Solid Tissue Normal" = "#4575b4"),
                     name = NULL) +
  labs(x = expression("cg13672843 " * beta),
       y = expression(log[2]("VHL TPM" + 1)),
       title = "VHL methylation × expression",
       subtitle = sprintf("Tumor: Spearman ρ = %.2f (p = %.2g, n = %d) | Normal: ρ = %.2f (p = %.2g, n = %d)",
                          cor_results$spearman_rho[1], cor_results$spearman_p[1], cor_results$n[1],
                          cor_results$spearman_rho[2], cor_results$spearman_p[2], cor_results$n[2])) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

p_box <- ggplot(tumor_df, aes(x = meth_group, y = vhl_log2, fill = meth_group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("Methylation-high (Q4)" = "#d73027",
                                "Middle (Q2-Q3)" = "#fee090",
                                "Methylation-low (Q1)" = "#4575b4"),
                     guide = "none") +
  labs(x = NULL,
       y = expression(log[2]("VHL TPM" + 1)),
       title = "VHL expression by cg13672843 methylation quartile (tumors only)",
       subtitle = sprintf("Wilcoxon (Q4 < Q1) p = %.3g", wt$p.value)) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 8))

combined <- p_scatter / p_box + plot_layout(heights = c(1, 1))
ggsave(file.path(OUT_DIR, "vhl_meth_expr.png"),
       combined, width = 7, height = 7, dpi = 150)
ggsave(file.path(OUT_DIR, "vhl_meth_expr_scatter.png"),
       p_scatter, width = 6, height = 4, dpi = 150)
ggsave(file.path(OUT_DIR, "vhl_meth_expr_box.png"),
       p_box, width = 6, height = 4, dpi = 150)

# --- Step 7: Save the joined data --------------------------------------------
write.csv(joined,    file.path(OUT_DIR, "vhl_meth_expr_joined.csv"),    row.names = FALSE)
write.csv(tumor_df,  file.path(OUT_DIR, "vhl_meth_expr_tumor_only.csv"), row.names = FALSE)

cat("\nDone. Outputs in", OUT_DIR, "\n")
sessionInfo()
