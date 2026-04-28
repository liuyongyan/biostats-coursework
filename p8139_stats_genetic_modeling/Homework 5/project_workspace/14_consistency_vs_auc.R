# Empirical check: are pct_consistent and disc_auc redundant, or complementary?
#
# Tests two claims:
#   (a) If they're perfectly correlated, one filter is redundant and can be dropped
#   (b) If they diverge, they capture different aspects and should both stay

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
HC_DIR <- file.path(WORKDIR, "hitclass_out")
OUT_DIR <- file.path(WORKDIR, "hitclass_out")  # save plot next to existing outputs

tt <- read.csv(file.path(HC_DIR, "all_cpgs_with_summary.csv"),
               stringsAsFactors = FALSE)
cat("Loaded", nrow(tt), "CpGs\n")

# Subset to FDR<0.05 (the only ones that could enter any Layer)
fdr_sig <- tt[tt$adj.P.Val < 0.05, ]
cat("FDR<0.05:", nrow(fdr_sig), "\n")

# --- Correlation ------------------------------------------------------------
cat("\n=== Correlation between pct_consistent and disc_auc (among FDR<0.05) ===\n")
pearson_r  <- cor(fdr_sig$pct_consistent, fdr_sig$disc_auc, method = "pearson")
spearman_r <- cor(fdr_sig$pct_consistent, fdr_sig$disc_auc, method = "spearman")
cat(sprintf("Pearson  r = %.4f\n", pearson_r))
cat(sprintf("Spearman ρ = %.4f\n", spearman_r))

# --- Cross-tabulation of filter agreement -----------------------------------
cat("\n=== Cross-tab of filter pass (among FDR<0.05) ===\n")
fdr_sig$pct_pass  <- fdr_sig$pct_consistent > 0.8
fdr_sig$auc_pass  <- fdr_sig$disc_auc > 0.85
cross <- with(fdr_sig, table(pct_pass, auc_pass))
print(cross)

n_total  <- nrow(fdr_sig)
n_pct    <- sum(fdr_sig$pct_pass)
n_auc    <- sum(fdr_sig$auc_pass)
n_both   <- sum(fdr_sig$pct_pass & fdr_sig$auc_pass)
n_only_p <- sum(fdr_sig$pct_pass & !fdr_sig$auc_pass)
n_only_a <- sum(!fdr_sig$pct_pass & fdr_sig$auc_pass)
cat(sprintf("\nFDR<0.05 total:              %d\n", n_total))
cat(sprintf("  Passes pct_consistent>0.8: %d (%.1f%%)\n", n_pct, 100*n_pct/n_total))
cat(sprintf("  Passes disc_auc>0.85:      %d (%.1f%%)\n", n_auc, 100*n_auc/n_total))
cat(sprintf("  Passes both:               %d (%.1f%%)\n", n_both, 100*n_both/n_total))
cat(sprintf("  Passes pct only (NOT auc): %d (%.1f%%)  ← filtered OUT by AUC\n",
            n_only_p, 100*n_only_p/n_total))
cat(sprintf("  Passes auc only (NOT pct): %d (%.1f%%)  ← filtered OUT by pct_consistent\n",
            n_only_a, 100*n_only_a/n_total))

# --- Jaccard similarity of the two filter sets ------------------------------
jaccard <- n_both / (n_pct + n_auc - n_both)
cat(sprintf("\nJaccard similarity of two filter sets: %.3f\n", jaccard))
cat("  (1.0 = identical filters, 0.5 = half overlap, etc.)\n")

# --- 2D scatter plot --------------------------------------------------------
# Subsample to 50K points for plotting speed
set.seed(8139)
plot_df <- fdr_sig[sample(nrow(fdr_sig), min(50000, nrow(fdr_sig))), ]
plot_df$quadrant <- with(plot_df, ifelse(pct_pass & auc_pass, "both",
                                 ifelse(pct_pass, "pct only",
                                 ifelse(auc_pass, "auc only", "neither"))))

p <- ggplot(plot_df, aes(pct_consistent, disc_auc, color = quadrant)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("both" = "#d73027", "pct only" = "#4575b4",
                                 "auc only" = "#fdae61", "neither" = "grey80")) +
  labs(x = "pct_consistent (within-patient directional consistency)",
       y = "disc_auc (unpaired discrimination)",
       title = "pct_consistent vs disc_auc — are they redundant?",
       subtitle = sprintf("Pearson r = %.3f  Spearman ρ = %.3f  Jaccard = %.3f",
                          pearson_r, spearman_r, jaccard)) +
  theme_bw()
ggsave(file.path(OUT_DIR, "consistency_vs_auc.png"), p,
       width = 8, height = 6, dpi = 150)

# --- Look at specific examples of divergence --------------------------------
cat("\n=== Examples: pct_consistent high but disc_auc low (pct_pass & !auc_pass) ===\n")
only_pct <- fdr_sig %>% filter(pct_pass & !auc_pass) %>%
  arrange(disc_auc) %>%
  select(cpg, chr, gene, pct_consistent, disc_auc, mean_delta_beta, adj.P.Val)
cat("Top 5 'high consistency but low AUC':\n")
print(head(only_pct, 5))

cat("\n=== Examples: disc_auc high but pct_consistent low (auc_pass & !pct_pass) ===\n")
only_auc <- fdr_sig %>% filter(!pct_pass & auc_pass) %>%
  arrange(pct_consistent) %>%
  select(cpg, chr, gene, pct_consistent, disc_auc, mean_delta_beta, adj.P.Val)
cat("Top 5 'high AUC but low consistency':\n")
print(head(only_auc, 5))

cat("\nDone. Plot at", file.path(OUT_DIR, "consistency_vs_auc.png"), "\n")
