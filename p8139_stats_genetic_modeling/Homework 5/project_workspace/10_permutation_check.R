# Permutation-based null calibration check for Q1.
# Randomly flip the sign of each patient's within-pair difference (tumor-normal
# or normal-tumor) and re-run the one-sample t-test per CpG. Under the null of
# "no systematic tumor-vs-normal effect", this should yield a uniform p-value
# distribution (lambda approximately 1). Any deviation would suggest confounding
# or methodology issues.
#
# We use 100 permutations for each of 20 permutation rounds -> compute mean
# lambda.

suppressPackageStartupMessages({
  library(limma)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "audit_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))

tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")
M_diff_true <- M[, bc_T] - M[, bc_N]

# Subset to 30k probes for speed (lambda is stable on a random subset)
subset_probes <- sample(nrow(M_diff_true), min(30000, nrow(M_diff_true)))
M_diff_sub <- M_diff_true[subset_probes, ]

lambda_from_diff <- function(Mdiff) {
  fit <- lmFit(Mdiff, design = matrix(1, ncol(Mdiff), 1,
                                      dimnames = list(NULL, "x")))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  median(qchisq(1 - tt$P.Value, df = 1)) / qchisq(0.5, df = 1)
}

# True lambda on the 30k subset
lambda_true <- lambda_from_diff(M_diff_sub)
cat(sprintf("Observed lambda on 30k subset (true labels): %.2f\n", lambda_true))

# Permutation: flip sign of each pair randomly, recompute lambda
B <- 20
lambdas_perm <- numeric(B)
for (b in seq_len(B)) {
  flip <- sample(c(-1, 1), ncol(M_diff_sub), replace = TRUE)
  M_perm <- sweep(M_diff_sub, 2, flip, `*`)
  lambdas_perm[b] <- lambda_from_diff(M_perm)
  cat(sprintf("  perm %2d: lambda = %.3f\n", b, lambdas_perm[b]))
}

cat("\nSummary:\n")
cat(sprintf("  True-label lambda (30k subset): %.2f\n", lambda_true))
cat(sprintf("  Permutation lambdas: mean=%.3f, sd=%.3f, range=[%.2f, %.2f]\n",
            mean(lambdas_perm), sd(lambdas_perm),
            min(lambdas_perm), max(lambdas_perm)))
cat("\nInterpretation:\n")
cat("  If the permutation lambdas cluster near 1.0, the methodology is\n")
cat("  well-calibrated and the observed high lambda under true labels\n")
cat("  reflects real biological signal, not confounding.\n")

writeLines(c(
  sprintf("true_label_lambda_subset = %.2f", lambda_true),
  sprintf("perm_lambda_mean = %.3f", mean(lambdas_perm)),
  sprintf("perm_lambda_sd = %.3f", sd(lambdas_perm)),
  sprintf("perm_lambda_range = [%.2f, %.2f]",
          min(lambdas_perm), max(lambdas_perm)),
  sprintf("n_permutations = %d", B),
  sprintf("n_probes_subset = %d", length(subset_probes))
), file.path(OUT_DIR, "permutation_lambda.txt"))

saveRDS(list(lambda_true = lambda_true, lambdas_perm = lambdas_perm),
        file.path(OUT_DIR, "permutation_lambda.rds"))
cat("\nDone.\n")
