# C3 — Threshold sensitivity analysis for the three-layer framework
# Goal: show that V1 literature-biomarker recovery and the VHL cg13672843
#       Layer-2 placement are stable across a reasonable grid of
#       (AUC threshold, pct_consistent threshold) values.
# This directly addresses §4.3 (iv) "principled but not optimized thresholds".

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(reshape2)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
OUT_DIR <- file.path(WORKDIR, "sensitivity_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load: per-CpG metrics + literature biomarker list -----------------------
es <- read.csv(file.path(WORKDIR, "hitclass_out", "all_cpgs_with_summary.csv"),
               stringsAsFactors = FALSE)
cat("Loaded:", nrow(es), "CpGs, columns:", paste(head(colnames(es),16), collapse=", "), "\n")
cat("FDR-sig:", sum(es$adj.P.Val < 0.05), "\n")

bm <- read.csv(file.path(WORKDIR, "validation_out", "literature_biomarker_list.csv"),
               stringsAsFactors = FALSE)
cat("Biomarker list:", nrow(bm), "genes\n")
PROM_FLAGS <- c("TSS1500","TSS200","5'UTR","1stExon")

has_prom <- function(feat_str) {
  if (is.na(feat_str) || feat_str == "") return(FALSE)
  any(unlist(strsplit(feat_str, ";")) %in% PROM_FLAGS)
}
match_gene <- function(gene_str, target) {
  if (is.na(gene_str) || gene_str == "") return(FALSE)
  toks <- toupper(unlist(strsplit(gene_str, ";")))
  # Aliases consistent with the main analysis
  if (any(target == "RASSF1A")) target <- c(target, "RASSF1")
  any(toks %in% toupper(target))
}

# Build per-biomarker summary of all annotated promoter CpGs
es$is_prom <- vapply(es$feat, has_prom, logical(1))
prom_pool  <- es %>% filter(adj.P.Val < 0.05, is_prom)  # FDR-sig + promoter
cat("FDR-sig promoter CpGs:", nrow(prom_pool), "\n")

# --- Sensitivity grid -------------------------------------------------------
auc_grid <- c(0.75, 0.80, 0.85, 0.90, 0.95)
pct_grid <- c(0.70, 0.75, 0.80, 0.85, 0.90)

results <- expand.grid(auc_thr = auc_grid, pct_thr = pct_grid)
results$layer1_size  <- NA_integer_
results$v1_recovery  <- NA_integer_      # genes with >=1 promoter CpG in L1 pool
results$vhl_cg13672843_layer <- NA_character_

for (i in seq_len(nrow(results))) {
  a  <- results$auc_thr[i]
  p_ <- results$pct_thr[i]
  L1_mask <- (es$adj.P.Val < 0.05) & (es$disc_auc > a) & (es$pct_consistent > p_)
  L1 <- es[L1_mask, ]
  results$layer1_size[i] <- nrow(L1)

  # Per-biomarker recovery: gene has >=1 promoter CpG in this Layer-1 pool
  recovered <- sum(vapply(bm$gene, function(g) {
    L1_g <- L1[L1$is_prom & vapply(L1$gene, match_gene, logical(1), target = g), ]
    nrow(L1_g) >= 1
  }, logical(1)))
  results$v1_recovery[i] <- recovered

  # cg13672843 layer placement
  if ("cg13672843" %in% es$cpg) {
    row <- es[es$cpg == "cg13672843", ]
    in_L1 <- (row$adj.P.Val < 0.05) && (row$disc_auc > a) && (row$pct_consistent > p_)
    in_L2 <- (row$adj.P.Val < 0.05) && (!in_L1) && (row$pct_responder >= 0.10)
    if (in_L1) results$vhl_cg13672843_layer[i] <- "L1"
    else if (in_L2) results$vhl_cg13672843_layer[i] <- "L2"
    else if (row$adj.P.Val < 0.05) results$vhl_cg13672843_layer[i] <- "L3"
    else results$vhl_cg13672843_layer[i] <- "n.s."
  }
}
cat("\nGrid results:\n"); print(results)
write.csv(results, file.path(OUT_DIR, "sensitivity_grid.csv"), row.names = FALSE)

# --- Heatmap A: V1 recovery -------------------------------------------------
results$is_default <- (results$auc_thr == 0.85) & (results$pct_thr == 0.80)
p_v1 <- ggplot(results, aes(x = factor(auc_thr), y = factor(pct_thr))) +
  geom_tile(aes(fill = v1_recovery), color = "white", linewidth = 0.4) +
  geom_text(aes(label = v1_recovery), size = 3, color = "black") +
  geom_text(data = subset(results, is_default),
            aes(label = "*"), size = 5, color = "red", vjust = -1.6, hjust = -0.8) +
  scale_fill_gradient2(low = "#4575b4", mid = "#fee090", high = "#d73027",
                       midpoint = 12, limits = c(min(results$v1_recovery)-0.5,
                                                  max(results$v1_recovery)+0.5),
                       name = "V1\nrecovery") +
  labs(x = "AUC threshold", y = "pct_consistent threshold",
       title = "V1 biomarker recovery (out of 20)",
       subtitle = "* = default thresholds (0.85, 0.80)") +
  theme_bw(base_size = 9)

# --- Heatmap B: Layer 1 size ------------------------------------------------
p_size <- ggplot(results, aes(x = factor(auc_thr), y = factor(pct_thr))) +
  geom_tile(aes(fill = log10(layer1_size)), color = "white", linewidth = 0.4) +
  geom_text(aes(label = format(layer1_size, big.mark = ",")), size = 2.4, color = "black") +
  geom_text(data = subset(results, is_default),
            aes(label = "*"), size = 5, color = "red", vjust = -1.6, hjust = -0.8) +
  scale_fill_gradient(low = "#fff7bc", high = "#cc4c02",
                      name = "log10(Layer-1\nCpG count)") +
  labs(x = "AUC threshold", y = "pct_consistent threshold",
       title = "Layer 1 size") +
  theme_bw(base_size = 9)

# --- Combined panel ---------------------------------------------------------
suppressPackageStartupMessages(library(patchwork))
combined <- p_v1 + p_size + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, "threshold_sensitivity.png"),
       combined, width = 8.5, height = 3.4, dpi = 150)

cat("\ncg13672843 layer placement table across grid:\n")
print(table(results$vhl_cg13672843_layer))
cat("\nV1 recovery range across grid:\n")
print(range(results$v1_recovery))
cat("\nDone.\n")
