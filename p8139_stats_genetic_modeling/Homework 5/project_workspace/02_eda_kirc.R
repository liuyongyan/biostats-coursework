# TCGA-KIRC methylation + clinical EDA.
#
# Purpose (Step 1, post-download):
#   - Load downloaded 450K beta-value files into a SummarizedExperiment.
#   - Cross-check sample barcodes (tumor / normal / paired patients).
#   - Report missingness at probe and sample level (no filtering yet).
#   - Report clinical metadata completeness.
#   - Quick unsupervised look (PCA on top-variance CpGs) to see if tumor and
#     normal separate — sanity check.
#
# Outputs written to eda_outputs/.
#
# We do NOT filter probes, normalize, or run differential tests here. All of
# that belongs in Step 3 with the chosen statistical plan. This script exists
# only to confirm the dataset is what we think it is before we commit to a
# research question.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DATA_DIR <- file.path(WORKDIR, "data_kirc")
OUT_DIR <- file.path(WORKDIR, "eda_outputs_kirc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

query_meth <- readRDS(file.path(DATA_DIR, "query_meth_450k.rds"))
clinical   <- readRDS(file.path(DATA_DIR, "clinical.rds"))

cat("=== GDCprepare: assemble beta-value SummarizedExperiment ===\n")
se_file <- file.path(DATA_DIR, "kirc_450k_SE.rds")
if (file.exists(se_file)) {
  cat("Loading cached SE from", se_file, "\n")
  se <- readRDS(se_file)
} else {
  se <- GDCprepare(
    query = query_meth,
    directory = file.path(DATA_DIR, "GDCdata"),
    save = FALSE,
    summarizedExperiment = TRUE
  )
  saveRDS(se, se_file)
}
cat("SE dims (CpG x sample):", dim(se)[1], "x", dim(se)[2], "\n")

# --- Sample-type decomposition --------------------------------------------
colmeta <- as.data.frame(colData(se))
bc <- colnames(se)
patient <- substr(bc, 1, 12)
sample_type_code <- substr(bc, 14, 15)
sample_type_label <- dplyr::recode(sample_type_code,
  "01" = "Primary Tumor",
  "02" = "Recurrent Tumor",
  "05" = "Additional New Primary",
  "06" = "Metastatic",
  "11" = "Solid Tissue Normal",
  .default = paste0("Other_", sample_type_code)
)
cat("\nSample-type counts:\n"); print(table(sample_type_label))

sample_df <- data.frame(
  barcode = bc,
  patient = patient,
  sample_type_code = sample_type_code,
  sample_type = sample_type_label,
  stringsAsFactors = FALSE
)
write.csv(sample_df, file.path(OUT_DIR, "sample_manifest.csv"), row.names = FALSE)

paired_patients <- intersect(
  sample_df$patient[sample_df$sample_type == "Primary Tumor"],
  sample_df$patient[sample_df$sample_type == "Solid Tissue Normal"]
)
cat("\nPaired patients (have BOTH 01 and 11 in 450K):", length(paired_patients), "\n")
writeLines(paired_patients, file.path(OUT_DIR, "paired_patient_ids.txt"))

# --- Per-sample missingness ------------------------------------------------
cat("\n=== Per-sample CpG missingness ===\n")
na_per_sample <- colSums(is.na(assay(se))) / nrow(se)
sample_na_df <- data.frame(barcode = colnames(se), frac_na = na_per_sample)
summary(sample_na_df$frac_na) |> print()
write.csv(sample_na_df, file.path(OUT_DIR, "sample_missingness.csv"), row.names = FALSE)

# --- Per-probe missingness ------------------------------------------------
cat("\n=== Per-probe sample missingness ===\n")
na_per_probe <- rowSums(is.na(assay(se))) / ncol(se)
cat("Probes with >0% missing:", sum(na_per_probe > 0), "\n")
cat("Probes with >5% missing:", sum(na_per_probe > 0.05), "\n")
cat("Probes with >20% missing:", sum(na_per_probe > 0.20), "\n")

# --- Clinical metadata completeness ---------------------------------------
cat("\n=== Clinical metadata key fields ===\n")
key_fields <- c(
  "submitter_id", "gender", "race", "ethnicity", "age_at_diagnosis",
  "vital_status", "days_to_death", "days_to_last_follow_up",
  "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
  "primary_diagnosis", "prior_malignancy", "prior_treatment",
  "tumor_grade"
)
present_fields <- intersect(key_fields, colnames(clinical))
missing_stats <- sapply(present_fields, function(f) mean(is.na(clinical[[f]]) | clinical[[f]] == ""))
clin_df <- data.frame(field = present_fields, frac_missing = round(missing_stats, 3))
print(clin_df)
write.csv(clin_df, file.path(OUT_DIR, "clinical_missingness.csv"), row.names = FALSE)

# Overlap: patients with methylation AND clinical
patients_with_meth <- unique(sample_df$patient)
patients_with_clin <- unique(clinical$submitter_id)
cat("\nPatients with methylation:", length(patients_with_meth), "\n")
cat("Patients with clinical:", length(patients_with_clin), "\n")
cat("Intersection:", length(intersect(patients_with_meth, patients_with_clin)), "\n")

# --- Unsupervised look: PCA on top-variance CpGs --------------------------
cat("\n=== PCA on 10k top-variance CpGs (complete rows only) ===\n")
beta <- assay(se)
keep_rows <- rowSums(is.na(beta)) == 0
cat("Probes with no missing:", sum(keep_rows), "\n")
beta_cc <- beta[keep_rows, , drop = FALSE]
vars <- matrixStats::rowVars(beta_cc)
top_idx <- order(vars, decreasing = TRUE)[seq_len(min(10000, length(vars)))]
beta_top <- beta_cc[top_idx, , drop = FALSE]

# PCA on samples (transpose so rows = samples)
pca <- prcomp(t(beta_top), scale. = FALSE, center = TRUE)
pca_df <- data.frame(
  barcode = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3]
)
pca_df$patient <- substr(pca_df$barcode, 1, 12)
pca_df$sample_type <- sample_df$sample_type[match(pca_df$barcode, sample_df$barcode)]
write.csv(pca_df, file.path(OUT_DIR, "pca_top10k.csv"), row.names = FALSE)

var_expl <- pca$sdev^2 / sum(pca$sdev^2)
cat("PC1-5 var explained:", round(var_expl[1:5], 3), "\n")

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = sample_type)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = "TCGA-KIRC 450K methylation — PCA (top 10k CpGs by variance)",
       subtitle = sprintf("PC1 var: %.1f%%  |  PC2 var: %.1f%%",
                          100*var_expl[1], 100*var_expl[2]))
ggsave(file.path(OUT_DIR, "pca_top10k_PC12.png"), p_pca, width = 7, height = 5, dpi = 150)

cat("\n=== EDA done. Outputs in:", OUT_DIR, "===\n")
print(list.files(OUT_DIR))

sink(file.path(OUT_DIR, "sessionInfo_eda.txt"))
print(sessionInfo())
sink()
