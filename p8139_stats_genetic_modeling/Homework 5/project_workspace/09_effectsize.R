# Supplementary: compute beta-unit effect sizes and biologically significant
# thresholds for Q1 and Q3. Statistical significance is not the same as
# biological significance — for methylation, |Delta beta| > 0.1 is a common
# "strong effect" threshold.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "effectsize_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# M -> beta
M2beta <- function(m) 2^m / (2^m + 1)

# ---- Q1 beta-scale effect sizes --------------------------------------------
# Re-load paired matrices to compute mean beta per sample type
M     <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta  <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
tt_q1 <- read.csv(file.path(WORKDIR, "q1_out", "toptable_q1.csv"),
                  stringsAsFactors = FALSE)

tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")
M_T <- M[, bc_T]; M_N <- M[, bc_N]

mean_beta_T <- rowMeans(M2beta(M_T))
mean_beta_N <- rowMeans(M2beta(M_N))
delta_beta  <- mean_beta_T - mean_beta_N

tt_q1$mean_beta_tumor  <- mean_beta_T[tt_q1$cpg]
tt_q1$mean_beta_normal <- mean_beta_N[tt_q1$cpg]
tt_q1$delta_beta       <- delta_beta[tt_q1$cpg]

# Summaries
n <- nrow(tt_q1)
fdr_sig <- tt_q1$adj.P.Val < 0.05
strong  <- fdr_sig & abs(tt_q1$delta_beta) > 0.1
very_strong <- fdr_sig & abs(tt_q1$delta_beta) > 0.2
cat("Q1 summaries:\n")
cat("  total probes:                     ", n, "\n")
cat("  FDR<0.05:                         ", sum(fdr_sig), sprintf("(%.1f%%)\n", 100*mean(fdr_sig)))
cat("  FDR<0.05 AND |delta beta|>0.1:    ", sum(strong),  sprintf("(%.1f%%)\n", 100*mean(strong)))
cat("  FDR<0.05 AND |delta beta|>0.2:    ", sum(very_strong), sprintf("(%.1f%%)\n", 100*mean(very_strong)))
cat("  Bonferroni AND |delta beta|>0.1:  ",
    sum(tt_q1$P.Value < 0.05/n & abs(tt_q1$delta_beta) > 0.1), "\n")
cat("  Among strong hits:\n")
cat("    hyper (delta_beta>0.1):         ",
    sum(fdr_sig & tt_q1$delta_beta > 0.1), "\n")
cat("    hypo  (delta_beta<-0.1):        ",
    sum(fdr_sig & tt_q1$delta_beta < -0.1), "\n")

write.csv(tt_q1, file.path(OUT_DIR, "q1_with_beta.csv"), row.names = FALSE)

# Top 30 by |t| among strong hits, with gene annotation
top_strong <- tt_q1[fdr_sig & abs(tt_q1$delta_beta) > 0.1, ]
top_strong <- top_strong[order(-abs(top_strong$t)), ][1:min(30, nrow(top_strong)), ]
write.csv(top_strong[, c("cpg","chr","pos","gene","feat",
                         "delta_beta","mean_beta_tumor","mean_beta_normal",
                         "logFC","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "q1_top30_strong.csv"), row.names = FALSE)

# Effect-size distribution plot
p_es <- ggplot(tt_q1, aes(delta_beta)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = c(-0.1, 0.1), color = "red", linetype = "dashed") +
  labs(x = "Delta beta (tumor - normal)", y = "# CpGs",
       title = "Q1: effect-size distribution (paired beta difference)") +
  theme_bw()
ggsave(file.path(OUT_DIR, "q1_deltabeta_hist.png"), p_es, width = 7, height = 4, dpi = 150)

# Scatter: Δβ vs -log10(p) — volcano on β scale
p_volc_b <- ggplot(tt_q1, aes(delta_beta, -log10(P.Value),
                              color = ifelse(strong, "strong", ifelse(fdr_sig, "FDR only","n.s.")))) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_vline(xintercept = c(-0.1, 0.1), color = "red", linetype = "dashed") +
  scale_color_manual(values = c("n.s."="grey70","FDR only"="dodgerblue","strong"="firebrick"),
                     name = "") +
  labs(x = "Delta beta", y = "-log10(p)",
       title = "Q1 beta-scale volcano") +
  theme_bw()
ggsave(file.path(OUT_DIR, "q1_betavolcano.png"), p_volc_b, width = 7.5, height = 5.5, dpi = 150)

# ---- Q3 beta-scale effect sizes --------------------------------------------
# For Q3 the coefficient is M per unit of grade. Convert to approximate β
# effect per grade-unit by evaluating the slope at mean β (via derivative).
tt_q3 <- read.csv(file.path(WORKDIR, "q3_out", "toptable_q3.csv"),
                  stringsAsFactors = FALSE)
# AveExpr is mean M across tumor samples in Q3; derive approx mean beta
mean_beta_q3 <- M2beta(tt_q3$AveExpr)
# dβ/dM at a point β is β(1-β)/log2(e) = β(1-β)*ln2/1 ... actually:
# β = 2^M / (2^M+1), dβ/dM = β(1-β)*ln2
# So approximate Δβ per unit grade ≈ logFC * β(1-β) * ln(2)
tt_q3$approx_beta_per_grade <- tt_q3$logFC * mean_beta_q3 * (1 - mean_beta_q3) * log(2)
# Across grades G1 to G4 (range 3) ≈ 3x
tt_q3$approx_total_beta_range <- tt_q3$approx_beta_per_grade * 3

n3 <- nrow(tt_q3)
fdr3 <- tt_q3$adj.P.Val < 0.05
strong3 <- fdr3 & abs(tt_q3$approx_total_beta_range) > 0.1
cat("\nQ3 summaries:\n")
cat("  FDR<0.05:                                 ", sum(fdr3),
    sprintf("(%.1f%%)\n", 100*mean(fdr3)))
cat("  FDR<0.05 AND |Δβ over G1-G4 range|>0.1:   ", sum(strong3), "\n")

write.csv(tt_q3, file.path(OUT_DIR, "q3_with_beta.csv"), row.names = FALSE)
top_strong3 <- tt_q3[strong3, ]
top_strong3 <- top_strong3[order(top_strong3$P.Value), ][1:min(30, nrow(top_strong3)), ]
write.csv(top_strong3[, c("cpg","chr","pos","gene","feat",
                          "approx_total_beta_range","logFC","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "q3_top30_strong.csv"), row.names = FALSE)

cat("\nDone. Outputs in", OUT_DIR, "\n")
