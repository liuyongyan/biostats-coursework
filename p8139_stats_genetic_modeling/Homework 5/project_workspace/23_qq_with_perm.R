#!/usr/bin/env Rscript
# 23_qq_with_perm.R
# Build a Q-Q plot that visually shows observed inflation vs permutation null,
# so the calibration argument doesn't depend on subtitle text alone.
#
# Strategy:
#   - Compute the observed paired-DM moderated p-values on a 30k-probe subset.
#   - Run B = 20 sign-flip permutations on the same subset; collect full p-values.
#   - Plot observed Q-Q points + permutation 95% envelope (per-rank quantiles
#     across permutation reps) + the diagonal.
#
# Output: writes per-rank quantile data to silencing_out/qq_perm_envelope.rds
# and a stand-alone PNG so we can preview it before integrating into Fig 2.

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(ggplot2)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
SIL_DIR  <- file.path(WORKDIR, "silencing_out")
FIG_DIR  <- file.path(WORKDIR, "figures_final")
dir.create(SIL_DIR, showWarnings = FALSE, recursive = TRUE)

M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))

tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")
M_diff_true <- M[, bc_T] - M[, bc_N]
rm(M); gc(verbose = FALSE)

subset_probes <- sample(nrow(M_diff_true), min(30000, nrow(M_diff_true)))
M_diff_sub <- M_diff_true[subset_probes, ]
n <- nrow(M_diff_sub)

pvals_from_diff <- function(Mdiff) {
  fit <- lmFit(Mdiff, design = matrix(1, ncol(Mdiff), 1,
                                      dimnames = list(NULL, "x")))
  fit <- eBayes(fit)
  topTable(fit, coef = 1, number = Inf, sort.by = "none")$P.Value
}

cat("[obs] computing observed p-values on 30k subset...\n")
p_obs <- pvals_from_diff(M_diff_sub)
lambda_obs <- median(qchisq(1 - p_obs, df = 1)) / qchisq(0.5, df = 1)
cat(sprintf("  observed lambda_GC = %.2f\n", lambda_obs))

B <- 20
P_perm <- matrix(NA_real_, nrow = n, ncol = B)
for (b in seq_len(B)) {
  flip <- sample(c(-1, 1), ncol(M_diff_sub), replace = TRUE)
  M_perm <- sweep(M_diff_sub, 2, flip, `*`)
  P_perm[, b] <- pvals_from_diff(M_perm)
  lam_b <- median(qchisq(1 - P_perm[, b], df = 1)) / qchisq(0.5, df = 1)
  cat(sprintf("  perm %2d: lambda = %.3f\n", b, lam_b))
}

# Build per-rank quantiles across permutations.
# For each permutation, sort -log10(p) descending; that gives the observed
# values vs the same expected-rank axis.
sorted_obs  <- sort(-log10(p_obs), decreasing = TRUE)
sorted_perm <- apply(-log10(P_perm), 2, sort, decreasing = TRUE)
expected    <- -log10(ppoints(n))

# Trim to a manageable grid for plotting (every ~30th rank for plotting,
# all kept in rds for analysis).
qq_dt <- data.table(
  expected     = expected,
  observed     = sorted_obs,
  perm_med     = apply(sorted_perm, 1, median),
  perm_lo      = apply(sorted_perm, 1, quantile, probs = 0.025),
  perm_hi      = apply(sorted_perm, 1, quantile, probs = 0.975)
)
saveRDS(qq_dt,  file.path(SIL_DIR, "qq_perm_envelope.rds"))
saveRDS(list(p_obs = p_obs, P_perm = P_perm, lambda_obs = lambda_obs),
        file.path(SIL_DIR, "perm_pvalues.rds"))

# Standalone preview plot
theme_paper <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_text(size = base - 1, color = "grey25")
    )
}
qq_thin <- qq_dt[seq(1, .N, by = max(1, floor(.N / 4000)))]
lambda_perm_mean <- mean(apply(P_perm, 2, function(p)
  median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)))

p_qq <- ggplot(qq_thin, aes(x = expected)) +
  geom_ribbon(aes(ymin = perm_lo, ymax = perm_hi),
              fill = "#a6cee3", alpha = 0.55) +
  geom_line(aes(y = perm_med), color = "#1f78b4", linewidth = 0.4) +
  geom_point(aes(y = observed), size = 0.5, alpha = 0.8, color = "grey15") +
  geom_abline(slope = 1, intercept = 0,
              color = "grey50", linewidth = 0.3, linetype = "dashed") +
  labs(x = expression(expected ~ -log[10](italic(p))),
       y = expression(observed ~ -log[10](italic(p))),
       title = "Q-Q plot: observed vs permutation null",
       subtitle = sprintf(
         "Observed λ_GC = %.1f (real labels); permutation-null λ_GC = %.2f (mean of %d sign-flip reps); 30k-probe subset.",
         lambda_obs, lambda_perm_mean, B)) +
  theme_paper(10)
ggsave(file.path(FIG_DIR, "fig02b_qq_preview.png"), p_qq,
       width = 6.5, height = 4.5, dpi = 220)

cat(sprintf("\nWritten:\n  %s\n  %s\n",
            file.path(SIL_DIR, "qq_perm_envelope.rds"),
            file.path(FIG_DIR, "fig02b_qq_preview.png")))
