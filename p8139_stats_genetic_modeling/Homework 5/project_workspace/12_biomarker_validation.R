# Stage 3 — three-way validation of the top-50 biomarker panel.
#
#   V1 Concurrent: compare top-50 panel + Layer-1 pool against 20-biomarker
#                  literature list (pre-specified with expected layers).
#   V2 Predictive: patient-level 5-fold CV with logistic regression;
#                  compare top-50 CV AUC vs L1-baseline (random 50 from
#                  FDR<0.05) and L2-baseline (random 50 from FDR<0.05 ∩
#                  |Δβ|>0.1), 30 replicates each.
#   V3 Biological: GSEA + ORA on top-50 mapped genes.
#
# All thresholds/assays pre-specified in plan.md §5b.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pROC)          # for AUC
  library(fgsea)
  library(msigdbr)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
HC_DIR  <- file.path(WORKDIR, "hitclass_out")
OUT_DIR <- file.path(WORKDIR, "validation_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load -------------------------------------------------------------------
M     <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta  <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
anno_ <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
top50 <- read.csv(file.path(HC_DIR, "top50_panel.csv"), stringsAsFactors = FALSE)
layer1 <- read.csv(file.path(HC_DIR, "layer1_set.csv"), stringsAsFactors = FALSE)
layer2 <- read.csv(file.path(HC_DIR, "layer2_set.csv"), stringsAsFactors = FALSE)
tt_all <- read.csv(file.path(HC_DIR, "all_cpgs_with_summary.csv"),
                   stringsAsFactors = FALSE)

# ============================================================================
# V1: Concurrent validity — literature biomarker list
# ============================================================================

cat("\n===============================================\n")
cat("V1 — Concurrent validity (literature overlap)\n")
cat("===============================================\n")

# Pre-specified biomarker list (see plan.md §5b)
lit <- data.frame(
  gene = c("SFRP1","SFRP2","RASSF1A","GATA5","DKK3",
           "BNC1","COL14A1","PCDH17","SLIT2","SFRP4","SFRP5",
           "VHL","CDKN2A","APAF1","TIMP3","NEFH","UCHL1","PITX2","SPG20","RPRM"),
  tier = c("A","A","A","A","A",
           "B","B","B","B","B","B",
           "A","B","B","B","B","B","B","B","C"),
  expected_layer = c("L1","L1","L1","L1","L1",
                     "L1","L1","L1","L1","L2","L2",
                     "L2","L2","L2","L2","L2","L2","L2","L2","not_detected"),
  stringsAsFactors = FALSE
)
write.csv(lit, file.path(OUT_DIR, "literature_biomarker_list.csv"), row.names = FALSE)
cat("Literature list (n=", nrow(lit), "):\n"); print(lit)

# Note: RASSF1A is typically annotated as just "RASSF1" in 450K annotation
# (the 450K UCSC gene names use RASSF1 which includes RASSF1A isoform).
# CDKN2A also has several aliases (INK4A, p16). We do alias-aware matching.
aliases <- list(
  RASSF1A = c("RASSF1", "RASSF1A")  # RASSF1A is the transcript-specific name; 450K annotation uses "RASSF1"
  # Note: do NOT alias CDKN2A with CDKN2B — they are distinct tumor-suppressor genes at the same locus
  # but biologically and annotation-wise independent.
)

# Helper: does a probe's gene string match any of the target gene names?
match_gene <- function(gene_str, target_genes) {
  if (is.na(gene_str) || gene_str == "") return(FALSE)
  names_in_probe <- unique(unlist(strsplit(gene_str, ";")))
  any(toupper(names_in_probe) %in% toupper(target_genes))
}

# For each literature biomarker, find all CpGs in its promoter region
promoter_flags <- c("TSS1500","TSS200","5'UTR","1stExon")
has_promoter <- function(feat_str) {
  if (is.na(feat_str) || feat_str == "") return(FALSE)
  any(unlist(strsplit(feat_str, ";")) %in% promoter_flags)
}

cat("\nChecking each literature biomarker against Top-50 and Layer 1 pool...\n")
v1 <- data.frame()
for (i in seq_len(nrow(lit))) {
  g <- lit$gene[i]
  target <- if (g %in% names(aliases)) aliases[[g]] else g

  # All probes assigned to this gene
  has_gene <- vapply(tt_all$gene, match_gene, logical(1), target_genes = target)
  has_prom <- vapply(tt_all$feat, has_promoter, logical(1))
  gene_prom_cpgs <- tt_all[has_gene & has_prom, ]

  # Is any of them in top-50?
  in_top50 <- any(gene_prom_cpgs$cpg %in% top50$cpg)
  # Is any in layer 1 (gate-pass pool)?
  in_layer1 <- any(gene_prom_cpgs$layer == "L1_biomarker_candidate")
  # Is any in layer 2?
  in_layer2 <- any(gene_prom_cpgs$layer == "L2_subgroup_driven")
  # Is any reaching FDR<0.05 at all?
  in_fdr <- any(gene_prom_cpgs$adj.P.Val < 0.05, na.rm = TRUE)
  # Direction of most-significant promoter CpG
  best_hit <- gene_prom_cpgs[which.min(gene_prom_cpgs$adj.P.Val), ]
  if (nrow(best_hit) == 0) best_hit <- data.frame(mean_delta_beta = NA,
                                                   adj.P.Val = NA,
                                                   layer = "missing")
  v1 <- rbind(v1, data.frame(
    gene = g,
    tier = lit$tier[i],
    expected_layer = lit$expected_layer[i],
    n_promoter_cpgs = nrow(gene_prom_cpgs),
    best_delta_beta = best_hit$mean_delta_beta[1],
    best_fdr = best_hit$adj.P.Val[1],
    best_cpg_layer = best_hit$layer[1],
    in_top50 = in_top50,
    in_layer1_pool = in_layer1,
    in_layer2_pool = in_layer2,
    any_fdr_sig = in_fdr,
    stringsAsFactors = FALSE
  ))
}
print(v1)

# V1 summary metrics
cat("\n--- V1 summary ---\n")
n_expected_L1 <- sum(v1$expected_layer == "L1")
n_expected_L2 <- sum(v1$expected_layer == "L2")

# Assignment agreement: for each gene, does observed layer match expected?
v1$observed_layer <- ifelse(v1$in_layer1_pool, "L1",
                            ifelse(v1$in_layer2_pool, "L2",
                                   ifelse(v1$any_fdr_sig, "L3", "not_detected")))
v1$prediction_correct <- v1$expected_layer == v1$observed_layer
write.csv(v1, file.path(OUT_DIR, "v1_literature_overlap.csv"), row.names = FALSE)

cat(sprintf("Total literature biomarkers: %d\n", nrow(v1)))
cat(sprintf("  ... any in top-50:              %d / %d\n",
            sum(v1$in_top50), nrow(v1)))
cat(sprintf("  ... any in Layer 1 pool:        %d / %d\n",
            sum(v1$in_layer1_pool), nrow(v1)))
cat(sprintf("  ... any in Layer 2 pool:        %d / %d\n",
            sum(v1$in_layer2_pool), nrow(v1)))
cat(sprintf("  ... any FDR<0.05 (L1+L2+L3):    %d / %d\n",
            sum(v1$any_fdr_sig), nrow(v1)))
cat("\n")
cat(sprintf("By expected layer (predicted vs observed):\n"))
cat(sprintf("  Expected L1 (n=%d), found in Layer 1 pool:  %d\n",
            n_expected_L1, sum(v1$expected_layer=="L1" & v1$in_layer1_pool)))
cat(sprintf("  Expected L2 (n=%d), found in Layer 2 pool:  %d\n",
            n_expected_L2, sum(v1$expected_layer=="L2" & v1$in_layer2_pool)))

cat(sprintf("\nOverall prediction accuracy: %d/%d = %.0f%%\n",
            sum(v1$prediction_correct), nrow(v1),
            100*mean(v1$prediction_correct)))

# Hypergeometric test: is Layer 1 enriched for literature biomarkers?
universe_genes <- unique(unlist(strsplit(paste(tt_all$gene, collapse=";"), ";")))
universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
layer1_genes <- unique(unlist(strsplit(paste(
  tt_all$gene[tt_all$layer == "L1_biomarker_candidate"], collapse=";"), ";")))
layer1_genes <- layer1_genes[!is.na(layer1_genes) & layer1_genes != ""]

# Flatten aliases back to the canonical gene names for the universe test
lit_genes_expanded <- unlist(lapply(seq_len(nrow(lit)), function(i) {
  g <- lit$gene[i]
  if (g %in% names(aliases)) aliases[[g]] else g
}))
lit_in_universe <- intersect(lit_genes_expanded, universe_genes)
lit_in_layer1 <- intersect(lit_genes_expanded, layer1_genes)

N <- length(universe_genes)
K <- length(layer1_genes)
n <- length(lit_in_universe)
k <- length(lit_in_layer1)
hyper_p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

cat(sprintf("\nHypergeometric test:\n"))
cat(sprintf("  Universe genes:        %d\n", N))
cat(sprintf("  Layer 1 genes:         %d\n", K))
cat(sprintf("  Literature in universe: %d\n", n))
cat(sprintf("  Literature in Layer 1:  %d\n", k))
cat(sprintf("  Hypergeometric p:      %.3g\n", hyper_p))

v1_summary <- list(
  # Total counts across all 20 literature biomarkers regardless of expected layer
  n_total_biomarkers = nrow(v1),
  n_total_in_top50 = sum(v1$in_top50),
  n_total_in_layer1_pool = sum(v1$in_layer1_pool),
  n_total_in_layer2_pool = sum(v1$in_layer2_pool),
  n_total_any_fdr = sum(v1$any_fdr_sig),
  # Predicted-layer-matched counts (expected L1 actually in L1, etc.)
  n_expected_L1 = n_expected_L1,
  n_expected_L1_in_top50 = sum(v1$expected_layer=="L1" & v1$in_top50),
  n_expected_L1_in_layer1_pool = sum(v1$expected_layer=="L1" & v1$in_layer1_pool),
  n_expected_L2 = n_expected_L2,
  n_expected_L2_in_layer2_pool = sum(v1$expected_layer=="L2" & v1$in_layer2_pool),
  # Overall layer prediction accuracy
  overall_accuracy = mean(v1$prediction_correct),
  # Hypergeometric enrichment
  hyper_p_layer1 = hyper_p,
  hyper_k_layer1 = k,
  hyper_n_universe = n
)
saveRDS(v1_summary, file.path(OUT_DIR, "v1_summary.rds"))

# ============================================================================
# V2: Predictive validity — 5-fold CV AUC
# ============================================================================

cat("\n===============================================\n")
cat("V2 — Predictive validity (5-fold CV)\n")
cat("===============================================\n")

# Only use paired samples (320 samples, 160 patients) for balanced CV
tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")
paired_bc <- c(bc_T, bc_N)
paired_patient <- c(paired_pat, paired_pat)
paired_label <- c(rep(1, length(bc_T)), rep(0, length(bc_N)))  # 1=tumor, 0=normal

# Also try including unpaired tumors (324-160=164) and any unpaired normals
# For now, stick with paired 320 samples.
cat("CV samples:", length(paired_bc), "(", length(paired_pat), "patients × 2)\n")

# Patient-level 5-fold CV
k_folds <- 5
patient_fold <- sample(rep(1:k_folds, length.out = length(paired_pat)))
# map patient fold → per-sample fold
sample_fold <- patient_fold[match(paired_patient, paired_pat)]

# Run one CV given a CpG panel
run_cv_auc <- function(cpg_panel) {
  M_panel <- M[cpg_panel, paired_bc, drop = FALSE]
  # Each row = sample, columns = CpGs
  X <- t(M_panel)
  preds_all <- numeric(length(paired_bc))
  for (f in 1:k_folds) {
    is_test <- sample_fold == f
    fit <- suppressWarnings(
      glm(paired_label[!is_test] ~ ., data = as.data.frame(X[!is_test, , drop = FALSE]),
          family = binomial())
    )
    preds_all[is_test] <- predict(fit, newdata = as.data.frame(X[is_test, , drop = FALSE]),
                                   type = "response")
  }
  auc_obj <- suppressMessages(roc(paired_label, preds_all, quiet = TRUE))
  as.numeric(pROC::auc(auc_obj))
}

# --- Top-50 panel -----------------------------------------------------------
cat("\nCV on top-50 panel...\n")
top50_cpgs <- top50$cpg
# sanity: all are in M matrix
stopifnot(all(top50_cpgs %in% rownames(M)))
top50_auc <- run_cv_auc(top50_cpgs)
cat(sprintf("Top-50 CV AUC = %.4f\n", top50_auc))

# --- L1 baseline: random 50 from FDR<0.05, 30 reps --------------------------
cat("\nBaseline L1: random 50 from FDR<0.05 (30 reps)...\n")
fdr_pool <- tt_all$cpg[tt_all$adj.P.Val < 0.05]
l1_aucs <- numeric(30)
for (r in 1:30) {
  random50 <- sample(fdr_pool, 50)
  l1_aucs[r] <- run_cv_auc(random50)
  cat(sprintf("  rep %2d: AUC = %.4f\n", r, l1_aucs[r]))
}

# --- L2 baseline: random 50 from FDR<0.05 ∩ |Δβ|>0.1, 30 reps ---------------
cat("\nBaseline L2: random 50 from FDR<0.05 ∩ |Δβ|>0.1 (30 reps)...\n")
strong_pool <- tt_all$cpg[tt_all$adj.P.Val < 0.05 & abs(tt_all$mean_delta_beta) > 0.1]
l2_aucs <- numeric(30)
for (r in 1:30) {
  random50 <- sample(strong_pool, 50)
  l2_aucs[r] <- run_cv_auc(random50)
  cat(sprintf("  rep %2d: AUC = %.4f\n", r, l2_aucs[r]))
}

# --- Summary statistics -----------------------------------------------------
v2_summary <- data.frame(
  set = c("Top-50 (gate+AUC rank)", "L1 baseline (FDR<0.05)", "L2 baseline (+|Δβ|>0.1)"),
  auc_mean = c(top50_auc, mean(l1_aucs), mean(l2_aucs)),
  auc_sd   = c(NA, sd(l1_aucs), sd(l2_aucs)),
  auc_p025 = c(NA, quantile(l1_aucs, 0.025), quantile(l2_aucs, 0.025)),
  auc_p975 = c(NA, quantile(l1_aucs, 0.975), quantile(l2_aucs, 0.975))
)
cat("\n--- V2 summary ---\n")
print(v2_summary)

# Permutation-style p-value: how often did baselines beat top-50?
p_vs_L1 <- mean(l1_aucs >= top50_auc)
p_vs_L2 <- mean(l2_aucs >= top50_auc)
cat(sprintf("\nP(L1 baseline >= top-50): %.3f\n", p_vs_L1))
cat(sprintf("P(L2 baseline >= top-50): %.3f\n", p_vs_L2))

saveRDS(list(top50_auc = top50_auc, l1_aucs = l1_aucs, l2_aucs = l2_aucs,
             p_vs_L1 = p_vs_L1, p_vs_L2 = p_vs_L2, summary = v2_summary),
        file.path(OUT_DIR, "v2_cv_results.rds"))
write.csv(v2_summary, file.path(OUT_DIR, "v2_cv_summary.csv"), row.names = FALSE)

# Plot: AUC distributions
auc_df_long <- rbind(
  data.frame(set = "Top-50 (gate+AUC)", auc = top50_auc),
  data.frame(set = "L1 baseline",       auc = l1_aucs),
  data.frame(set = "L2 baseline",       auc = l2_aucs)
)
p_auc <- ggplot(auc_df_long, aes(x = set, y = auc, color = set)) +
  geom_boxplot(data = auc_df_long[auc_df_long$set != "Top-50 (gate+AUC)", ],
               outlier.size = 0.5) +
  geom_point(size = 3) +
  geom_hline(yintercept = top50_auc, linetype = "dashed", color = "red") +
  labs(y = "5-fold CV AUC",
       title = "V2: Predictive validity of top-50 biomarker panel",
       subtitle = sprintf("Top-50 AUC = %.4f | p vs L1 = %.3f, p vs L2 = %.3f",
                          top50_auc, p_vs_L1, p_vs_L2)) +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "v2_auc_comparison.png"), p_auc, width = 7, height = 5, dpi = 150)

# ============================================================================
# V3: Biological validity — GSEA + ORA on top-50 genes
# ============================================================================

cat("\n===============================================\n")
cat("V3 — Biological validity (GSEA / ORA on top-50)\n")
cat("===============================================\n")

# Build ranked list for GSEA using top-50 genes vs rest of Layer-1
# Actually, for a top-50 panel, ORA (hypergeometric) is more appropriate
# than GSEA. Do both.

first_gene <- function(s) {
  s[is.na(s) | s == ""] <- NA
  sapply(strsplit(s, ";"), function(x) if (length(x) == 0) NA_character_ else x[1])
}

# For top-50: get gene list
top50_genes <- unique(first_gene(top50$gene))
top50_genes <- top50_genes[!is.na(top50_genes)]
cat("Top-50 mapped genes:", length(top50_genes), "\n"); print(top50_genes)

# Universe: all genes in tt_all with gene annotation
tt_all$gene1 <- first_gene(tt_all$gene)
universe_g <- unique(tt_all$gene1[!is.na(tt_all$gene1)])
cat("Universe genes:", length(universe_g), "\n")

# Load Hallmark + KEGG sets
hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
nm_col <- "gs_name"; sym_col <- "gene_symbol"
hallmark_list <- split(hallmark[[sym_col]], hallmark[[nm_col]])

kegg <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                                   subcollection = "CP:KEGG_LEGACY"),
                 error = function(e) NULL)
if (is.null(kegg) || nrow(kegg) == 0) {
  kegg <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                                     subcollection = "CP:KEGG"),
                   error = function(e) NULL)
}
if (!is.null(kegg) && nrow(kegg) > 0) {
  kegg_list <- split(kegg[[sym_col]], kegg[[nm_col]])
  set_list <- c(hallmark_list, kegg_list)
} else {
  set_list <- hallmark_list
}
cat("Gene sets loaded:", length(set_list), "\n")

# ORA (hypergeometric) on top-50 genes
cat("\nRunning ORA on top-50 genes...\n")
ora_res <- do.call(rbind, lapply(names(set_list), function(nm) {
  gs <- intersect(set_list[[nm]], universe_g)
  if (length(gs) < 5) return(NULL)
  k  <- length(intersect(top50_genes, gs))
  K  <- length(gs)
  N  <- length(universe_g)
  n  <- length(top50_genes)
  p  <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(pathway = nm, overlap = k, set_size = K,
             universe = N, panel_size = n, p = p)
}))
ora_res$padj <- p.adjust(ora_res$p, method = "BH")
ora_res <- ora_res[order(ora_res$p), ]
write.csv(ora_res, file.path(OUT_DIR, "v3_ora_top50.csv"), row.names = FALSE)
cat("Top 10 enriched pathways (ORA):\n")
print(head(ora_res[, c("pathway","overlap","set_size","p","padj")], 10))

# Check specifically for KIRC-related pathways
kirc_keywords <- c("RENAL","HYPOXIA","WNT","VHL","HIF","CLEAR_CELL")
kirc_hits <- ora_res[grepl(paste(kirc_keywords, collapse="|"),
                            toupper(ora_res$pathway)), ]
cat("\nKIRC-relevant pathway hits:\n")
print(kirc_hits)

v3_summary <- list(
  n_pathways_fdr05 = sum(ora_res$padj < 0.05, na.rm = TRUE),
  n_pathways_fdr25 = sum(ora_res$padj < 0.25, na.rm = TRUE),
  top_pathway = ora_res$pathway[1],
  top_pathway_padj = ora_res$padj[1],
  kirc_relevant_hits = nrow(kirc_hits),
  kirc_relevant_min_padj = if (nrow(kirc_hits) > 0) min(kirc_hits$padj) else NA
)
cat("\n--- V3 summary ---\n")
print(v3_summary)
saveRDS(v3_summary, file.path(OUT_DIR, "v3_summary.rds"))

# ============================================================================
# Overall decision
# ============================================================================

cat("\n===============================================\n")
cat("Overall decision (per plan.md §5b)\n")
cat("===============================================\n")

V1_pass_strict <- v1_summary$n_total_in_top50 >= 12      # any literature bmk in top-50 >= 12/20
V1_pass_lenient <- v1_summary$n_total_in_layer1_pool >= 12  # any in Layer 1 pool >= 12/20
V2_pass_strict <- top50_auc > quantile(l2_aucs, 0.95)
V3_pass <- v3_summary$kirc_relevant_hits > 0 &&
            !is.na(v3_summary$kirc_relevant_min_padj) &&
            v3_summary$kirc_relevant_min_padj < 0.05

cat(sprintf("V1 pass (strict — top-50 recovers ≥12/20 L1): %s\n", V1_pass_strict))
cat(sprintf("V1 pass (lenient — Layer-1 pool recovers ≥12/20 L1): %s\n", V1_pass_lenient))
cat(sprintf("V2 pass (top-50 AUC > L2-baseline 95%%ile):  %s\n", V2_pass_strict))
cat(sprintf("V3 pass (≥1 KIRC-related pathway padj<0.05): %s\n", V3_pass))

overall <- if (V1_pass_strict && V2_pass_strict && V3_pass) "SUCCESS" else
           if (v1_summary$n_total_in_top50 < 8 && (top50_auc - mean(l2_aucs)) < 0.02) "FAILURE" else
           "PARTIAL"
cat(sprintf("\n>>> OVERALL DECISION: %s <<<\n", overall))

writeLines(c(
  sprintf("Date: %s", Sys.time()),
  sprintf("V1 strict pass:  %s (%d/20 in top-50)",
          V1_pass_strict, v1_summary$n_total_in_top50),
  sprintf("V1 lenient pass: %s (%d/20 in Layer 1 pool; %d/20 in Layer 2 pool)",
          V1_pass_lenient, v1_summary$n_total_in_layer1_pool,
          v1_summary$n_total_in_layer2_pool),
  sprintf("V2 pass:         %s (top-50 AUC=%.4f, L2 mean=%.4f)",
          V2_pass_strict, top50_auc, mean(l2_aucs)),
  sprintf("V3 pass:         %s (KIRC-relevant hits: %d)",
          V3_pass, v3_summary$kirc_relevant_hits),
  sprintf("Overall:         %s", overall)
), file.path(OUT_DIR, "decision.txt"))

sink(file.path(OUT_DIR, "sessionInfo_12.txt"))
print(sessionInfo())
sink()
cat("\nDone. Outputs in", OUT_DIR, "\n")
