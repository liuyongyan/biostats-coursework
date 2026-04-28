# Step 3b — Q1 primary: paired tumor vs adjacent-normal differential methylation.
#
# Rewritten 2026-04-23: switched from limma+duplicateCorrelation (too slow on
# 320k probes × 320 samples) to the canonical textbook approach for balanced
# paired designs: within-pair differences + moderated one-sample t-test.
#
# Why this is rigorous:
#   For patient i with M_tumor_i and M_normal_i, the paired difference is
#     D_i = M_tumor_i - M_normal_i = tumor_effect + (e_tumor_i - e_normal_i)
#   Any patient-level covariate (age, gender, race, TSS) is identical in both
#   aliquots and CANCELS in D_i. Hence the difference model needs no patient-
#   level covariate adjustment — the pairing does it for us.
#
#   Aliquot-level covariates (plate) can differ within a pair and do NOT
#   cancel. We handle this in a sensitivity analysis below (exclude cross-plate
#   pairs).
#
# Pre-registered design in plan.md §5 called for covariate adjustment with
# limma+duplicateCorrelation; this is a minor amendment driven by feasibility.
# Documented here and will be noted in the final report's Methods section.

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "q1_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

M     <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta  <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
anno_ <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))

cat("Input:", nrow(M), "probes x", ncol(M), "samples\n")

# ---- Select paired samples, one tumor + one normal per patient -------------
tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
cat("Paired patients:", length(paired_pat), "\n")

# For each patient, pick first tumor and first normal aliquot
pick_bc <- function(pat, type) {
  bc <- meta$barcode[meta$patient == pat & meta$sample_type == type]
  bc[1]
}
bc_tumor  <- vapply(paired_pat, pick_bc, character(1), type = "Primary Tumor")
bc_normal <- vapply(paired_pat, pick_bc, character(1), type = "Solid Tissue Normal")
stopifnot(!any(is.na(bc_tumor)), !any(is.na(bc_normal)))

# Build aligned matrices
M_T <- M[, bc_tumor];  colnames(M_T) <- paired_pat
M_N <- M[, bc_normal]; colnames(M_N) <- paired_pat
stopifnot(identical(colnames(M_T), colnames(M_N)))
cat("Aligned M_T:", ncol(M_T), "  M_N:", ncol(M_N), "\n")

# Within-pair differences
M_diff <- M_T - M_N
cat("M_diff:", nrow(M_diff), "x", ncol(M_diff), "\n")

# ---- Primary: moderated one-sample t on differences ------------------------
cat("\n=== PRIMARY Q1 model: within-pair differences, moderated one-sample t ===\n")
design <- matrix(1, nrow = ncol(M_diff), ncol = 1, dimnames = list(NULL, "tumorVsNormal"))
fit <- lmFit(M_diff, design = design)
fit <- eBayes(fit, trend = FALSE, robust = FALSE)
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
tt$cpg  <- rownames(tt)
tt$chr  <- anno_[tt$cpg, "chr"]
tt$pos  <- anno_[tt$cpg, "pos"]
tt$gene <- anno_[tt$cpg, "UCSC_RefGene_Name"]
tt$feat <- anno_[tt$cpg, "UCSC_RefGene_Group"]

n_probes <- nrow(tt)
cat("N probes:", n_probes, "\n")
cat("BH FDR < 0.05:", sum(tt$adj.P.Val < 0.05), "\n")
cat("BH FDR < 0.01:", sum(tt$adj.P.Val < 0.01), "\n")
cat("Bonferroni:   ", sum(tt$P.Value < 0.05 / n_probes), "\n")
cat("Hyper FDR<0.05:", sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), "\n")
cat("Hypo  FDR<0.05:", sum(tt$adj.P.Val < 0.05 & tt$logFC < 0), "\n")

lambda <- median(qchisq(1 - tt$P.Value, df = 1)) / qchisq(0.5, df = 1)
cat("lambda_GC:", round(lambda, 3), "\n")

write.csv(tt[, c("cpg","chr","pos","gene","feat",
                 "logFC","AveExpr","t","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "toptable_q1.csv"), row.names = FALSE)
saveRDS(list(fit = fit, tt = tt, lambda = lambda,
             n_pairs = length(paired_pat), n_probes = n_probes),
        file.path(OUT_DIR, "q1_fit.rds"))

# ---- Diagnostic plots -------------------------------------------------------
pvals <- sort(tt$P.Value); n_p <- length(pvals)
qq_df <- data.frame(expected = -log10(ppoints(n_p)), observed = -log10(pvals))
p_qq <- ggplot(qq_df, aes(expected, observed)) +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "-log10(expected p)", y = "-log10(observed p)",
       title = sprintf("Q1 Q-Q plot  (lambda = %.3f, n pairs = %d)",
                       lambda, length(paired_pat))) +
  theme_bw()
ggsave(file.path(OUT_DIR, "qq_q1.png"), p_qq, width = 5.5, height = 5.5, dpi = 150)

tt$sig <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1, "sig |logFC|>1",
          ifelse(tt$adj.P.Val < 0.05, "FDR<0.05", "n.s."))
p_volc <- ggplot(tt, aes(logFC, -log10(P.Value), color = sig)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = c("n.s." = "grey70", "FDR<0.05" = "dodgerblue",
                                 "sig |logFC|>1" = "firebrick")) +
  labs(x = "Mean within-pair M difference (tumor - normal)",
       y = "-log10(p)", title = "Q1 volcano (paired)") +
  theme_bw()
ggsave(file.path(OUT_DIR, "volcano_q1.png"), p_volc, width = 7, height = 5.5, dpi = 150)

# Manhattan
tt_m <- tt[!is.na(tt$chr), ]
tt_m$chr_num <- suppressWarnings(as.numeric(sub("chr", "", tt_m$chr)))
tt_m$chr_num[tt_m$chr == "chrX"] <- 23; tt_m$chr_num[tt_m$chr == "chrY"] <- 24
tt_m <- tt_m %>% arrange(chr_num, pos)
tt_m$idx <- seq_len(nrow(tt_m))
chr_mids <- tt_m %>% group_by(chr) %>% summarise(mid = median(idx), .groups = "drop")
p_man <- ggplot(tt_m, aes(idx, -log10(P.Value), color = factor(chr_num %% 2))) +
  geom_point(size = 0.25, alpha = 0.5) +
  scale_color_manual(values = c("0" = "grey30", "1" = "grey60"), guide = "none") +
  scale_x_continuous(breaks = chr_mids$mid, labels = sub("chr","",chr_mids$chr)) +
  geom_hline(yintercept = -log10(0.05 / n_probes), color = "red", linetype = "dashed") +
  labs(x = "chromosome", y = "-log10(p)", title = "Q1 genome-wide p-values") +
  theme_bw() + theme(axis.text.x = element_text(size = 7))
ggsave(file.path(OUT_DIR, "manhattan_q1.png"), p_man, width = 10, height = 4, dpi = 150)

# ---- Sensitivity S1: unpaired 2-sample with covariates ---------------------
cat("\n=== Sensitivity S1: unpaired independent 2-sample, covariate adjusted ===\n")
meta_u <- meta[meta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]
M_u <- M[, meta_u$barcode]
meta_u$tumorVsNormal <- factor(meta_u$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor"))
meta_u$gender_f <- factor(meta_u$gender)
meta_u$tss_f    <- factor(meta_u$tss)
meta_u$age_c    <- as.numeric(meta_u$age_at_diagnosis) / 365.25
meta_u$age_c    <- meta_u$age_c - mean(meta_u$age_c, na.rm = TRUE)
meta_u$age_c[is.na(meta_u$age_c)] <- 0

design_s1 <- model.matrix(~ tumorVsNormal + gender_f + age_c + tss_f, data = meta_u)
is_singular <- apply(design_s1, 2, function(x) length(unique(x))) == 1
design_s1 <- design_s1[, !is_singular, drop = FALSE]
fit_s1 <- lmFit(M_u, design_s1)
fit_s1 <- eBayes(fit_s1)
tt_s1 <- topTable(fit_s1, coef = "tumorVsNormalPrimary Tumor",
                  number = Inf, sort.by = "P")
tt_s1$cpg <- rownames(tt_s1)
cat("S1 n=", ncol(M_u), " FDR<0.05:", sum(tt_s1$adj.P.Val < 0.05), "\n")
lambda_s1 <- median(qchisq(1 - tt_s1$P.Value, df = 1)) / qchisq(0.5, df = 1)
cat("S1 lambda:", round(lambda_s1, 3), "\n")
write.csv(tt_s1[, c("cpg","logFC","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "toptable_q1_s1_unpaired.csv"), row.names = FALSE)

# ---- Sensitivity S2: exclude pairs where tumor/normal on different plates --
cat("\n=== Sensitivity S2: paired, excluding cross-plate pairs ===\n")
plate_tumor  <- meta$plate[match(bc_tumor,  meta$barcode)]
plate_normal <- meta$plate[match(bc_normal, meta$barcode)]
cross_plate <- plate_tumor != plate_normal
cat("Pairs on same plate: ", sum(!cross_plate),
    "  cross-plate:", sum(cross_plate), "\n")
if (sum(!cross_plate) >= 30) {
  M_T_s <- M_T[, !cross_plate]; M_N_s <- M_N[, !cross_plate]
  M_diff_s <- M_T_s - M_N_s
  fit_s2 <- lmFit(M_diff_s, matrix(1, ncol(M_diff_s), 1,
                                   dimnames = list(NULL, "tumorVsNormal")))
  fit_s2 <- eBayes(fit_s2)
  tt_s2 <- topTable(fit_s2, coef = 1, number = Inf, sort.by = "P")
  tt_s2$cpg <- rownames(tt_s2)
  cat("S2 n pairs:", ncol(M_diff_s), " FDR<0.05:", sum(tt_s2$adj.P.Val < 0.05), "\n")
  lambda_s2 <- median(qchisq(1 - tt_s2$P.Value, df = 1)) / qchisq(0.5, df = 1)
  cat("S2 lambda:", round(lambda_s2, 3), "\n")
  write.csv(tt_s2[, c("cpg","logFC","P.Value","adj.P.Val")],
            file.path(OUT_DIR, "toptable_q1_s2_sameplatepairs.csv"), row.names = FALSE)
} else {
  cat("S2 skipped: too few same-plate pairs (", sum(!cross_plate), ")\n")
  tt_s2 <- NULL
}

# ---- Overlap of primary vs sensitivities -----------------------------------
top_pri <- tt$cpg[tt$adj.P.Val < 0.05]
top_s1  <- tt_s1$cpg[tt_s1$adj.P.Val < 0.05]
cat("\nSig overlap primary & S1 (unpaired):", length(intersect(top_pri, top_s1)),
    " (pri:", length(top_pri), ", s1:", length(top_s1), ")\n")
if (!is.null(tt_s2)) {
  top_s2 <- tt_s2$cpg[tt_s2$adj.P.Val < 0.05]
  cat("Sig overlap primary & S2 (same-plate pairs):",
      length(intersect(top_pri, top_s2)),
      " (pri:", length(top_pri), ", s2:", length(top_s2), ")\n")
}

sink(file.path(OUT_DIR, "sessionInfo_q1.txt"))
print(sessionInfo())
sink()
cat("\nQ1 done. Outputs in", OUT_DIR, "\n")
