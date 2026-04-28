# PC1-PC10 diagnostic for TCGA-KIRC 450K methylation.
#
# For each top-10 PC, test association with:
#   - Sample type (01 vs 11) -- expected PC1 driver
#   - Clinical: gender, race, ethnicity, age, stage, T/N/M, tumor_grade,
#     vital_status, primary_diagnosis
#   - Barcode-derived batch: TSS (tissue source site), plate, center, portion/analyte
#
# For each (PC, covariate):
#   - categorical covariate  -> one-way ANOVA F-test p-value and adj R^2
#   - numeric covariate      -> Pearson correlation (r and p)
#
# Outputs:
#   scree_top20.png  -- scree plot
#   pc_pairs.png     -- PC1-2, 3-4, 5-6 with sample_type coloring
#   pc_covariate_association.csv
#   pc_covariate_heatmap.png

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DATA_DIR <- file.path(WORKDIR, "data_kirc")
OUT_DIR  <- file.path(WORKDIR, "eda_outputs_kirc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Loading cached SE and clinical ...\n")
se <- readRDS(file.path(DATA_DIR, "kirc_450k_SE.rds"))
clinical <- readRDS(file.path(DATA_DIR, "clinical.rds"))

# -- Build sample-level metadata table from barcodes + clinical --
bc <- colnames(se)
meta <- data.frame(
  barcode        = bc,
  patient        = substr(bc, 1, 12),
  tss            = substr(bc, 6, 7),         # tissue source site
  sample_type_cd = substr(bc, 14, 15),
  vial           = substr(bc, 16, 16),
  portion        = substr(bc, 18, 19),
  analyte        = substr(bc, 20, 20),
  plate          = substr(bc, 22, 25),
  center         = substr(bc, 27, 28),
  stringsAsFactors = FALSE
)
meta$sample_type <- dplyr::recode(meta$sample_type_cd,
  "01" = "Primary Tumor", "05" = "Additional New Primary", "11" = "Solid Tissue Normal",
  .default = paste0("Other_", meta$sample_type_cd))

# Merge clinical by patient
clin_keep <- c("submitter_id", "gender", "race", "ethnicity", "age_at_diagnosis",
               "vital_status", "ajcc_pathologic_stage", "ajcc_pathologic_t",
               "ajcc_pathologic_n", "ajcc_pathologic_m", "tumor_grade",
               "primary_diagnosis")
clin_keep <- intersect(clin_keep, colnames(clinical))
clin_small <- clinical[, clin_keep, drop = FALSE]
meta <- dplyr::left_join(meta, clin_small,
                         by = c("patient" = "submitter_id"))

# -- PCA on top-10k variance CpGs, complete rows --
cat("Computing PCA ...\n")
beta <- assay(se)
keep_rows <- rowSums(is.na(beta)) == 0
beta_cc <- beta[keep_rows, , drop = FALSE]
vars <- matrixStats::rowVars(beta_cc)
top_idx <- order(vars, decreasing = TRUE)[seq_len(min(10000, length(vars)))]
beta_top <- beta_cc[top_idx, , drop = FALSE]
pca <- prcomp(t(beta_top), scale. = FALSE, center = TRUE)

var_expl <- pca$sdev^2 / sum(pca$sdev^2)
cat("PC1-PC20 var explained:\n")
print(round(var_expl[1:20] * 100, 2))

# -- Scree plot (top 20) --
scree_df <- data.frame(PC = 1:20, var = 100 * var_expl[1:20])
p_scree <- ggplot(scree_df, aes(PC, var)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", var)), vjust = -0.4, size = 3) +
  scale_x_continuous(breaks = 1:20) +
  labs(y = "Variance explained (%)", title = "Scree plot — KIRC methylation PCA (top 10k variance CpGs)") +
  theme_bw()
ggsave(file.path(OUT_DIR, "scree_top20.png"), p_scree, width = 9, height = 4.5, dpi = 150)

# -- Save PC scores table --
pc_scores <- as.data.frame(pca$x[, 1:10])
colnames(pc_scores) <- paste0("PC", 1:10)
pc_scores$barcode <- rownames(pca$x)
meta <- dplyr::left_join(meta, pc_scores, by = "barcode")
write.csv(meta, file.path(OUT_DIR, "meta_with_pc1_10.csv"), row.names = FALSE)

# -- Pairwise PC plots colored by sample_type --
long_plot_df <- meta
long_plot_df$sample_type <- factor(long_plot_df$sample_type)
g12 <- ggplot(long_plot_df, aes(PC1, PC2, color = sample_type)) + geom_point(alpha = .7) + theme_bw()
g34 <- ggplot(long_plot_df, aes(PC3, PC4, color = sample_type)) + geom_point(alpha = .7) + theme_bw()
g56 <- ggplot(long_plot_df, aes(PC5, PC6, color = sample_type)) + geom_point(alpha = .7) + theme_bw()
# Save combined
library(grid); library(gridExtra)
g_all <- arrangeGrob(g12, g34, g56, ncol = 3)
ggsave(file.path(OUT_DIR, "pc_pairs_sample_type.png"), g_all, width = 15, height = 4.5, dpi = 150)

# -- Association: each PC vs each covariate --
covariates <- c("sample_type", "tss", "plate", "center", "portion", "analyte",
                "gender", "race", "ethnicity", "age_at_diagnosis",
                "vital_status", "ajcc_pathologic_stage", "ajcc_pathologic_t",
                "ajcc_pathologic_n", "ajcc_pathologic_m", "tumor_grade",
                "primary_diagnosis")
covariates <- intersect(covariates, colnames(meta))

pc_names <- paste0("PC", 1:10)

results <- list()
for (pc in pc_names) {
  y <- meta[[pc]]
  for (v in covariates) {
    x <- meta[[v]]
    if (all(is.na(x))) next
    if (is.numeric(x)) {
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) < 5) next
      r <- cor(x[ok], y[ok])
      p <- cor.test(x[ok], y[ok])$p.value
      results[[length(results)+1]] <- data.frame(
        PC = pc, covariate = v, type = "numeric",
        stat = r, p_value = p, r2 = r^2
      )
    } else {
      # categorical: one-way anova of PC ~ x
      df <- data.frame(y = y, x = factor(x))
      df <- df[!is.na(df$y) & !is.na(df$x) & df$x != "", ]
      if (nlevels(droplevels(df$x)) < 2 || nrow(df) < 5) next
      df$x <- droplevels(df$x)
      fit <- lm(y ~ x, data = df)
      s <- summary(fit)
      f <- s$fstatistic
      pval <- if (!is.null(f)) pf(f[1], f[2], f[3], lower.tail = FALSE) else NA
      results[[length(results)+1]] <- data.frame(
        PC = pc, covariate = v, type = "categorical",
        stat = NA, p_value = pval, r2 = s$r.squared
      )
    }
  }
}
res_df <- do.call(rbind, results)
res_df$PC <- factor(res_df$PC, levels = pc_names)
res_df$neglog10p <- -log10(pmax(res_df$p_value, 1e-300))
write.csv(res_df, file.path(OUT_DIR, "pc_covariate_association.csv"), row.names = FALSE)

# -- Heatmap of -log10(p) --
hm_df <- res_df %>% select(PC, covariate, r2) %>% spread(covariate, r2)
mat <- as.matrix(hm_df[, -1])
rownames(mat) <- hm_df$PC

long_hm <- res_df
p_hm <- ggplot(long_hm, aes(x = PC, y = covariate, fill = r2)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r2)), color = "white", size = 2.8) +
  scale_fill_gradient(low = "#eeeeee", high = "#b2182b", limits = c(0,1), name = "R^2") +
  theme_bw() +
  labs(title = "PC ~ covariate R^2 heatmap (top 10 PCs)")
ggsave(file.path(OUT_DIR, "pc_covariate_heatmap.png"), p_hm, width = 9, height = 6.5, dpi = 150)

cat("\nTop associations per PC (by R^2):\n")
res_df %>% group_by(PC) %>% slice_max(r2, n = 3) %>% as.data.frame() %>% print()

cat("\nOutputs in:", OUT_DIR, "\n")
