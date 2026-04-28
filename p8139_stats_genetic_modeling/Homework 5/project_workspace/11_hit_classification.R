# Stage 1-2 of biomarker discovery pipeline (revised 2026-04-23, Option B).
#
# Gate (three filters, all pre-specified):
#   (G1) FDR < 0.05               — statistical significance
#   (G2) disc_auc > 0.85          — per-CpG unpaired AUC, scale-invariant
#                                    effect-size; replaces prior |Δβ|>0.1
#   (G3) pct_consistent > 0.8     — ≥128/160 pairs same direction
#
# Rank (within gate):
#   per-CpG unpaired AUC (tumor vs normal), descending
#
# Select:
#   top-50 with one-CpG-per-gene constraint
#
# Layer taxonomy (revised):
#   Layer 1 = gate-passing
#   Layer 2 = FDR<0.05 AND NOT Layer 1 AND pct_responder (|Δβ|>0.2) >= 0.10
#   Layer 3 = FDR<0.05, NOT Layer 1, NOT Layer 2

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "hitclass_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load -------------------------------------------------------------------
cat("Loading inputs...\n")
M     <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta  <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
anno_ <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
tt_q1 <- read.csv(file.path(WORKDIR, "q1_out", "toptable_q1.csv"),
                  stringsAsFactors = FALSE)

# --- Paired matrices --------------------------------------------------------
tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")

M2beta <- function(m) 2^m / (2^m + 1)
beta_T <- M2beta(M[, bc_T])
beta_N <- M2beta(M[, bc_N])
delta_beta_mat <- beta_T - beta_N      # probes x 160 pairs
cat("Δβ matrix:", nrow(delta_beta_mat), "x", ncol(delta_beta_mat), "\n")

# --- Per-CpG summaries ------------------------------------------------------
cat("\nComputing per-CpG summaries...\n")
n_pairs <- ncol(delta_beta_mat)
mean_db <- rowMeans(delta_beta_mat)
n_pos <- rowSums(delta_beta_mat > 0)
n_neg <- rowSums(delta_beta_mat < 0)
pct_consistent <- pmax(n_pos, n_neg) / n_pairs
n_responder <- rowSums(abs(delta_beta_mat) > 0.2)
pct_responder <- n_responder / n_pairs

# --- Per-CpG unpaired AUC (tumor vs normal on all 484 samples) --------------
cat("Computing per-CpG unpaired AUC...\n")
# Use all primary tumor samples AND all solid normal samples
is_tumor <- meta$sample_type == "Primary Tumor"
is_norm  <- meta$sample_type == "Solid Tissue Normal"
cat("  tumor samples:", sum(is_tumor), " normal samples:", sum(is_norm), "\n")

# AUC via Mann-Whitney U: for each probe,
#   AUC = (sum_rank(tumor) - n_T*(n_T+1)/2) / (n_T * n_N)
# To handle ordering direction (hypermethylated vs hypomethylated), we use
# max(AUC, 1-AUC) so AUC is always >= 0.5 and reflects discrimination magnitude.
fast_auc_per_cpg <- function(M_tumor, M_normal) {
  # M_tumor: probes x n_T, M_normal: probes x n_N
  stopifnot(nrow(M_tumor) == nrow(M_normal))
  n_T <- ncol(M_tumor); n_N <- ncol(M_normal)
  n_all <- n_T + n_N
  # Combine and rank per probe
  M_all <- cbind(M_tumor, M_normal)
  # Rank per row; ties.method "average" (standard)
  ranks <- matrixStats::rowRanks(M_all, ties.method = "average")
  sum_tumor_ranks <- rowSums(ranks[, 1:n_T, drop = FALSE])
  U_tumor <- sum_tumor_ranks - n_T * (n_T + 1) / 2
  auc <- U_tumor / (n_T * n_N)
  # AUC reflects direction; symmetrize to discrimination magnitude
  disc_auc <- pmax(auc, 1 - auc)
  list(auc = auc, disc_auc = disc_auc)
}

auc_res <- fast_auc_per_cpg(M[, is_tumor], M[, is_norm])
tt <- tt_q1 %>%
  mutate(mean_delta_beta = mean_db[cpg],
         pct_consistent   = pct_consistent[cpg],
         pct_responder    = pct_responder[cpg],
         auc              = auc_res$auc[cpg],
         disc_auc         = auc_res$disc_auc[cpg])

# --- Sanity check -----------------------------------------------------------
cat("\nSanity: pct_consistent max =", max(tt$pct_consistent, na.rm=TRUE),
    "; min =", min(tt$pct_consistent, na.rm=TRUE), "\n")
cat("disc_auc range: [", round(min(tt$disc_auc, na.rm=TRUE),3), ",",
    round(max(tt$disc_auc, na.rm=TRUE),3), "]\n")

# --- Apply gate -------------------------------------------------------------
cat("\n=== Gate (FDR<0.05 + disc_auc>0.85 + pct_consistent>0.8) ===\n")
fdr_sig <- tt$adj.P.Val < 0.05
gate_ok <- fdr_sig &
  tt$disc_auc > 0.85 &
  tt$pct_consistent > 0.8

cat("After G1 (FDR<0.05):              ", sum(fdr_sig), "\n")
cat("After G2 (+disc_auc>0.85):        ", sum(fdr_sig & tt$disc_auc > 0.85), "\n")
cat("After G3 (+pct_consistent>0.8):   ", sum(gate_ok), "\n")

candidate_pool <- tt[gate_ok, ] %>%
  arrange(desc(disc_auc))
cat("Candidate pool size:", nrow(candidate_pool), "\n")

# --- Layer taxonomy ---------------------------------------------------------
# Layer 2 lower bound 0.10: ~5x the null expectation of P(|Δβ|>0.2) ≈ 0.01-0.02
tt$layer <- "not_sig"
tt$layer[fdr_sig & gate_ok] <- "L1_biomarker_candidate"
tt$layer[fdr_sig & !gate_ok & tt$pct_responder >= 0.10] <- "L2_subgroup_driven"
tt$layer[fdr_sig & tt$layer == "not_sig"] <- "L3_weak_pervasive"

cat("\nLayer counts:\n"); print(table(tt$layer))

# --- Top-50 panel with one-CpG-per-gene constraint --------------------------
cat("\n=== Top-50 selection with one-per-gene constraint ===\n")
first_gene <- function(s) {
  s[is.na(s) | s == ""] <- NA
  sapply(strsplit(s, ";"), function(x) if (length(x) == 0) NA_character_ else x[1])
}
candidate_pool$gene_primary <- first_gene(candidate_pool$gene)

selected <- list()
seen_genes <- character()
for (i in seq_len(nrow(candidate_pool))) {
  g <- candidate_pool$gene_primary[i]
  if (is.na(g)) {
    # Always allow unassigned probes (may be intergenic CpG islands)
    selected[[length(selected)+1]] <- candidate_pool[i, ]
  } else {
    if (!(g %in% seen_genes)) {
      selected[[length(selected)+1]] <- candidate_pool[i, ]
      seen_genes <- c(seen_genes, g)
    }
  }
  if (length(selected) >= 50) break
}
top50 <- do.call(rbind, selected)
cat("Top-50 panel: n =", nrow(top50), "; unique genes =",
    length(unique(top50$gene_primary[!is.na(top50$gene_primary)])), "\n")

# --- Save outputs -----------------------------------------------------------
write.csv(tt, file.path(OUT_DIR, "all_cpgs_with_summary.csv"), row.names = FALSE)
write.csv(candidate_pool, file.path(OUT_DIR, "candidate_pool.csv"), row.names = FALSE)
write.csv(top50[, c("cpg","gene_primary","chr","pos","feat",
                     "mean_delta_beta","pct_consistent","disc_auc",
                     "logFC","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "top50_panel.csv"), row.names = FALSE)

layer1 <- tt[tt$layer == "L1_biomarker_candidate", ]
layer2 <- tt[tt$layer == "L2_subgroup_driven", ]
layer3 <- tt[tt$layer == "L3_weak_pervasive", ]
write.csv(layer1, file.path(OUT_DIR, "layer1_set.csv"), row.names = FALSE)
write.csv(layer2, file.path(OUT_DIR, "layer2_set.csv"), row.names = FALSE)
write.csv(layer3, file.path(OUT_DIR, "layer3_set.csv"), row.names = FALSE)

# --- Plots ------------------------------------------------------------------
cat("\nBuilding plots...\n")

# 1. Attrition funnel
attrition <- data.frame(
  stage = factor(c("All tested","FDR<0.05","+disc_auc>0.85","+pct_consistent>0.8","Top-50 panel"),
                 levels = c("All tested","FDR<0.05","+disc_auc>0.85","+pct_consistent>0.8","Top-50 panel")),
  n = c(nrow(tt), sum(fdr_sig),
        sum(fdr_sig & tt$disc_auc > 0.85),
        sum(gate_ok), nrow(top50))
)
p_fun <- ggplot(attrition, aes(x = stage, y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = scales::comma(n)), vjust = -0.3, size = 3.5) +
  scale_y_log10(labels = scales::comma) +
  labs(y = "# CpGs (log scale)",
       title = "Biomarker attrition funnel (Q1 → gate → top-50)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT_DIR, "attrition_funnel.png"), p_fun, width = 7, height = 4.5, dpi = 150)

# 2. disc_auc vs pct_consistent, colored by gate pass
tt$gate_pass <- ifelse(gate_ok, "Layer 1 (gate pass)",
                       ifelse(fdr_sig, "FDR only", "n.s."))
p_gate <- ggplot(tt[sample(nrow(tt), min(100000, nrow(tt))), ],
                 aes(x = disc_auc, y = pct_consistent, color = gate_pass)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("n.s." = "grey80", "FDR only" = "dodgerblue",
                                 "Layer 1 (gate pass)" = "firebrick")) +
  labs(x = "disc_auc (per-CpG unpaired discrimination AUC)",
       y = "pct_consistent (majority-direction fraction)",
       title = "Gate visualization — Layer 1 = top-right quadrant",
       subtitle = sprintf("Gate passes: %d CpGs", sum(gate_ok))) +
  theme_bw()
ggsave(file.path(OUT_DIR, "gate_scatter.png"), p_gate, width = 8, height = 5.5, dpi = 150)

# 3. AUC distribution in Layer 1 vs random CpGs
set.seed(8139)
random_sample <- sample(which(!gate_ok), 10000)
auc_df <- rbind(
  data.frame(set = "Layer 1 (gate pass)", disc_auc = tt$disc_auc[gate_ok]),
  data.frame(set = "Random (non-gate)",   disc_auc = tt$disc_auc[random_sample])
)
p_auc <- ggplot(auc_df, aes(x = disc_auc, fill = set)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  labs(x = "Per-CpG discrimination AUC", y = "# CpGs",
       title = "AUC distribution: Layer 1 vs non-Layer-1 CpGs") +
  theme_bw()
ggsave(file.path(OUT_DIR, "auc_distribution.png"), p_auc, width = 8, height = 4.5, dpi = 150)

# --- Print summary ----------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Layer 1 (gate-pass):     ", sum(tt$layer == "L1_biomarker_candidate"), "\n")
cat("Layer 2 (subgroup):      ", sum(tt$layer == "L2_subgroup_driven"), "\n")
cat("Layer 3 (weak pervasive):", sum(tt$layer == "L3_weak_pervasive"), "\n")
cat("Top-50 panel size:       ", nrow(top50), "\n")
cat("Top-50 max disc_auc:     ", round(max(top50$disc_auc), 4), "\n")
cat("Top-50 min disc_auc:     ", round(min(top50$disc_auc), 4), "\n")
cat("Top-50 direction split:  hyper =", sum(top50$mean_delta_beta > 0),
    "; hypo =", sum(top50$mean_delta_beta < 0), "\n")

sink(file.path(OUT_DIR, "sessionInfo_11.txt"))
print(sessionInfo())
sink()
cat("\nDone. Outputs in", OUT_DIR, "\n")
