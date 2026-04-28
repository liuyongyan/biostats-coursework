# 10_vhl_subgroup.R
# Subgroup-aware VHL methylation × expression analysis (Fig 5).
#
# Background: a naive Pearson/Spearman scatter on the full cohort yields
# Pearson r = -0.27 (p = 1e-6) but Spearman ρ = -0.02 (p = 0.69). This
# split is the visual signature of a subgroup-driven pattern — a small
# number of high-methylation patients drive a strong linear signal that
# vanishes under rank correlation because the bulk of patients have
# β ≈ 0 and rank randomly among themselves.
#
# This script formalizes that intuition:
#   1. Define "responder" tumors as β >= 0.2 at cg13672843 (Layer 2 pct_responder threshold)
#   2. Compare VHL expression: responder vs non-responder tumors (Mann-Whitney)
#   3. Paired analysis: for patients with BOTH methylation AND RNA-seq pairs,
#      compute Δβ and Δlog2(VHL TPM); test whether methylation gain → expression loss

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(SummarizedExperiment)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "vhl_expr_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load joined data from previous step
joined <- read.csv(file.path(OUT_DIR, "vhl_meth_expr_joined.csv"), stringsAsFactors = FALSE)
tumor_df  <- joined %>% filter(sample_type == "Primary Tumor")
normal_df <- joined %>% filter(sample_type == "Solid Tissue Normal")

cat("Tumors:", nrow(tumor_df), " Normals:", nrow(normal_df), "\n")

# ---- (1) Subgroup definition: β >= 0.2 ------------------------------------
RESP_THRESH <- 0.2
tumor_df <- tumor_df %>%
  mutate(group = ifelse(beta >= RESP_THRESH, "Responder (β ≥ 0.2)", "Non-responder (β < 0.2)"))
n_resp <- sum(tumor_df$group == "Responder (β ≥ 0.2)")
cat(sprintf("Responder tumors (β >= %.2f): %d / %d = %.1f%%\n",
            RESP_THRESH, n_resp, nrow(tumor_df), 100*n_resp/nrow(tumor_df)))

resp <- tumor_df %>% filter(group == "Responder (β ≥ 0.2)")
nonresp <- tumor_df %>% filter(group == "Non-responder (β < 0.2)")
mw <- wilcox.test(resp$vhl_log2, nonresp$vhl_log2, alternative = "less")
ttest <- t.test(resp$vhl_log2, nonresp$vhl_log2, alternative = "less")
cat(sprintf("\nResponder VHL log2(TPM+1): mean=%.2f median=%.2f n=%d\n",
            mean(resp$vhl_log2), median(resp$vhl_log2), nrow(resp)))
cat(sprintf("Non-resp.  VHL log2(TPM+1): mean=%.2f median=%.2f n=%d\n",
            mean(nonresp$vhl_log2), median(nonresp$vhl_log2), nrow(nonresp)))
cat(sprintf("Wilcoxon (responder < non-resp): p = %.3g\n", mw$p.value))
cat(sprintf("t-test   (responder < non-resp): p = %.3g, mean diff = %.3f\n",
            ttest$p.value, mean(resp$vhl_log2) - mean(nonresp$vhl_log2)))

# Effect size: Cohen's d
sd_pool <- sqrt(((nrow(resp)-1)*var(resp$vhl_log2) + (nrow(nonresp)-1)*var(nonresp$vhl_log2)) /
                (nrow(resp) + nrow(nonresp) - 2))
cohen_d <- (mean(resp$vhl_log2) - mean(nonresp$vhl_log2)) / sd_pool
cat(sprintf("Cohen's d (responder vs non-responder): %.3f\n", cohen_d))

# ---- (2) Paired analysis: Δβ vs Δlog2(VHL TPM) ----------------------------
# A patient must have both tumor + normal aliquots in BOTH methylation AND RNA-seq

paired_pat <- intersect(
  tumor_df$patient,
  normal_df$patient
)
cat(sprintf("\nPatients with paired meth + RNA-seq tumor and normal: %d\n", length(paired_pat)))

paired_df <- data.frame(patient = paired_pat) %>%
  left_join(tumor_df  %>% select(patient, beta_T = beta, vhl_T = vhl_log2), by = "patient") %>%
  left_join(normal_df %>% select(patient, beta_N = beta, vhl_N = vhl_log2), by = "patient") %>%
  mutate(delta_beta = beta_T - beta_N,
         delta_vhl  = vhl_T - vhl_N)

# Pearson + Spearman on paired Δs
cor_pearson  <- cor.test(paired_df$delta_beta, paired_df$delta_vhl, method = "pearson")
cor_spearman <- cor.test(paired_df$delta_beta, paired_df$delta_vhl, method = "spearman", exact = FALSE)
cat(sprintf("\nPaired analysis: n = %d\n", nrow(paired_df)))
cat(sprintf("  Pearson:  r   = %.3f, p = %.3g\n", cor_pearson$estimate, cor_pearson$p.value))
cat(sprintf("  Spearman: rho = %.3f, p = %.3g\n", cor_spearman$estimate, cor_spearman$p.value))

# Δβ > 0.2 paired-responder subgroup
paired_df <- paired_df %>%
  mutate(paired_group = ifelse(delta_beta >= RESP_THRESH, "Δβ ≥ 0.2 (gain)",
                          ifelse(delta_beta <= -RESP_THRESH, "Δβ ≤ -0.2 (loss)", "|Δβ| < 0.2")))
cat("\nPaired group sizes:\n"); print(table(paired_df$paired_group))
gain    <- paired_df %>% filter(paired_group == "Δβ ≥ 0.2 (gain)")
neutral <- paired_df %>% filter(paired_group == "|Δβ| < 0.2")
if (nrow(gain) >= 5) {
  mw_p <- wilcox.test(gain$delta_vhl, neutral$delta_vhl, alternative = "less")
  cat(sprintf("\nPaired test: gain (n=%d, mean Δ-VHL=%.2f) vs neutral (n=%d, mean Δ-VHL=%.2f)\n",
              nrow(gain), mean(gain$delta_vhl), nrow(neutral), mean(neutral$delta_vhl)))
  cat(sprintf("  Wilcoxon (gain Δ-VHL < neutral Δ-VHL): p = %.3g\n", mw_p$p.value))
}

# ---- (3) Save numerical summary -------------------------------------------
summary_table <- data.frame(
  test = c("Tumor Pearson (all 318)",
           "Tumor Spearman (all 318)",
           "Responder vs Non-resp (Wilcoxon, β >= 0.2)",
           "Responder vs Non-resp (t-test)",
           "Paired Pearson (Δβ vs Δlog2-VHL)",
           "Paired Spearman (Δβ vs Δlog2-VHL)",
           "Paired Δβ-gain vs neutral (Wilcoxon)"),
  n = c(nrow(tumor_df), nrow(tumor_df),
        paste0(nrow(resp), "/", nrow(nonresp)),
        paste0(nrow(resp), "/", nrow(nonresp)),
        nrow(paired_df), nrow(paired_df),
        if(nrow(gain) >= 5) paste0(nrow(gain), "/", nrow(neutral)) else "n/a"),
  estimate = c(sprintf("r=%.3f",  cor.test(tumor_df$beta, tumor_df$vhl_log2, method="pearson")$estimate),
               sprintf("ρ=%.3f", cor.test(tumor_df$beta, tumor_df$vhl_log2, method="spearman", exact=FALSE)$estimate),
               sprintf("Δlog2 = %.2f", mean(resp$vhl_log2) - mean(nonresp$vhl_log2)),
               sprintf("Cohen's d = %.2f", cohen_d),
               sprintf("r=%.3f", cor_pearson$estimate),
               sprintf("ρ=%.3f", cor_spearman$estimate),
               if(nrow(gain) >= 5) sprintf("Δlog2 diff = %.2f", mean(gain$delta_vhl) - mean(neutral$delta_vhl)) else "n/a"),
  p_value = c(sprintf("%.3g", cor.test(tumor_df$beta, tumor_df$vhl_log2, method="pearson")$p.value),
              sprintf("%.3g", cor.test(tumor_df$beta, tumor_df$vhl_log2, method="spearman", exact=FALSE)$p.value),
              sprintf("%.3g", mw$p.value),
              sprintf("%.3g", ttest$p.value),
              sprintf("%.3g", cor_pearson$p.value),
              sprintf("%.3g", cor_spearman$p.value),
              if(nrow(gain) >= 5) sprintf("%.3g", mw_p$p.value) else "n/a")
)
print(summary_table)
write.csv(summary_table, file.path(OUT_DIR, "vhl_subgroup_summary.csv"), row.names = FALSE)

# ---- (4) Final figure -----------------------------------------------------
# Panel A: scatter highlighting responder subset
mean_resp    <- mean(resp$vhl_log2)
mean_nonresp <- mean(nonresp$vhl_log2)
x_max        <- max(tumor_df$beta, normal_df$beta, na.rm = TRUE)

p_scatter <- ggplot() +
  # Horizontal mean reference segments — drawn first so points sit on top.
  # The y-axis difference between these two lines is the actual finding.
  geom_segment(aes(x = 0,           xend = RESP_THRESH,
                   y = mean_nonresp, yend = mean_nonresp),
               color = "grey25", linewidth = 0.7) +
  geom_segment(aes(x = RESP_THRESH, xend = x_max + 0.01,
                   y = mean_resp,    yend = mean_resp),
               color = "#d73027", linewidth = 0.9) +
  # Points: normals (triangle, larger) then non-responders then responders
  geom_point(data = normal_df, aes(x = beta, y = vhl_log2),
             color = "#4575b4", size = 2.6, alpha = 0.9, shape = 17) +
  geom_point(data = nonresp, aes(x = beta, y = vhl_log2),
             color = "grey55", size = 1.4, alpha = 0.55) +
  geom_point(data = resp, aes(x = beta, y = vhl_log2),
             color = "#d73027", size = 2.2, alpha = 0.85) +
  geom_vline(xintercept = RESP_THRESH, linetype = "dashed", color = "grey30") +
  annotate("text", x = RESP_THRESH + 0.01, y = 1.85,
           label = "responder threshold (β = 0.2)",
           hjust = 0, size = 2.5, color = "grey30") +
  annotate("text", x = RESP_THRESH - 0.005, y = mean_nonresp,
           label = sprintf("non-responder mean = %.2f", mean_nonresp),
           hjust = 1, vjust = -0.5, size = 2.5, color = "grey25") +
  annotate("text", x = x_max + 0.005, y = mean_resp,
           label = sprintf("responder mean = %.2f", mean_resp),
           hjust = 1, vjust = -0.5, size = 2.5, color = "#d73027") +
  annotate("text", x = x_max, y = 5.7,
           label = sprintf("Tumors (β ≥ 0.2): n = %d (%.0f%%)\nWilcoxon p = %.2g, Cohen's d = %.2f",
                            nrow(resp), 100*nrow(resp)/nrow(tumor_df),
                            mw$p.value,
                            (mean_resp - mean_nonresp) /
                              sqrt(((nrow(resp)-1)*var(resp$vhl_log2) +
                                    (nrow(nonresp)-1)*var(nonresp$vhl_log2)) /
                                   (nrow(resp) + nrow(nonresp) - 2))),
           hjust = 1, vjust = 1, size = 2.7) +
  labs(x = expression("cg13672843 " * beta * " (tumor or normal)"),
       y = expression(log[2]("VHL TPM" + 1))) +
  theme_bw(base_size = 9)

# Panel B: paired Δβ vs Δlog2(VHL)
p_paired <- NULL
if (nrow(paired_df) >= 10) {
  p_paired <- ggplot(paired_df, aes(x = delta_beta, y = delta_vhl)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dotted") +
    geom_vline(xintercept = 0, color = "grey60", linetype = "dotted") +
    geom_vline(xintercept = c(-RESP_THRESH, RESP_THRESH), color = "grey30", linetype = "dashed") +
    geom_point(aes(color = paired_group), size = 1.8, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.4) +
    scale_color_manual(values = c("Δβ ≥ 0.2 (gain)"   = "#d73027",
                                   "Δβ ≤ -0.2 (loss)" = "#1a9850",
                                   "|Δβ| < 0.2"       = "grey60"),
                        name = NULL) +
    labs(x = expression(Delta * beta * " (tumor − normal) at cg13672843"),
         y = expression(Delta * log[2]("VHL TPM" + 1) * " (tumor − normal)"),
         title = "Per-patient methylation gain → expression loss (paired)",
         subtitle = sprintf("n = %d patients with both tumor and normal aliquots (meth + RNA-seq)\nPearson r = %.2f (p = %.2g) | Spearman ρ = %.2f (p = %.2g)",
                            nrow(paired_df), cor_pearson$estimate, cor_pearson$p.value,
                            cor_spearman$estimate, cor_spearman$p.value)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")
}

if (!is.null(p_paired)) {
  combined <- p_scatter / p_paired + plot_layout(heights = c(1, 1))
  ggsave(file.path(OUT_DIR, "vhl_subgroup_figure.png"), combined,
         width = 7, height = 7.5, dpi = 150)
  cat("\nSaved combined figure.\n")
}
# Always also save a single-panel scatter for paper use
ggsave(file.path(OUT_DIR, "vhl_subgroup_scatter_single.png"), p_scatter,
       width = 7.5, height = 3.6, dpi = 150)
cat("Saved single-panel scatter for paper figure use.\n")

cat("\nDone.\n")
