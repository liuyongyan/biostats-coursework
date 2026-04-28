# Step 3c — Q3 secondary: within-tumor CpG associations with tumor grade.
#
# Pre-registered (plan.md §5):
#   Input:  primary tumor samples (-01), n~324.
#   Test:   limma linear model M ~ ordinal(grade) + gender + age + TSS.
#           Coefficient on grade (treated as numeric ordinal, G1=1..G4=4).
#   MCC:    BH FDR q < 0.05.
#   Sensitivity: Kendall's tau per CpG on a subsample (too slow genome-wide),
#                or repeat with grade as factor (categorical) for comparison.

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "q3_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

M      <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta   <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
anno_  <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))

# Subset to primary tumors
keep <- meta$sample_type == "Primary Tumor"
meta_t <- meta[keep, ]
M_t <- M[, keep]
cat("Primary tumor samples:", ncol(M_t), "\n")

# ---- Grade encoding ---------------------------------------------------------
cat("\nRaw tumor_grade distribution:\n"); print(table(meta_t$tumor_grade, useNA = "ifany"))

# Recode to ordinal numeric: G1=1, G2=2, G3=3, G4=4
grade_map <- c("G1" = 1, "G2" = 2, "G3" = 3, "G4" = 4)
meta_t$grade_num <- grade_map[meta_t$tumor_grade]
drop_missing <- is.na(meta_t$grade_num)
cat("Dropping", sum(drop_missing), "samples with unknown grade (GX / NA)\n")
meta_t <- meta_t[!drop_missing, ]
M_t    <- M_t[, !drop_missing]
cat("After grade filter:", ncol(M_t), "samples\n")
print(table(meta_t$grade_num))

# ---- Covariate prep ---------------------------------------------------------
meta_t$gender_f <- factor(meta_t$gender)
meta_t$tss_f    <- factor(meta_t$tss)
meta_t$age_c    <- as.numeric(meta_t$age_at_diagnosis) / 365.25
meta_t$age_c    <- meta_t$age_c - mean(meta_t$age_c, na.rm = TRUE)
meta_t$age_c[is.na(meta_t$age_c)] <- 0

# ---- Primary: limma linear -------------------------------------------------
cat("\n=== PRIMARY Q3: limma M ~ ordinal(grade) + gender + age + TSS ===\n")
design <- model.matrix(~ grade_num + gender_f + age_c + tss_f, data = meta_t)
is_singular <- apply(design, 2, function(x) length(unique(x))) == 1
design <- design[, !is_singular, drop = FALSE]
cat("Design rank:", qr(design)$rank, "  ncol:", ncol(design), "\n")

fit <- lmFit(M_t, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "grade_num", number = Inf, sort.by = "P")
tt$cpg  <- rownames(tt)
tt$chr  <- anno_[tt$cpg, "chr"]
tt$pos  <- anno_[tt$cpg, "pos"]
tt$gene <- anno_[tt$cpg, "UCSC_RefGene_Name"]
tt$feat <- anno_[tt$cpg, "UCSC_RefGene_Group"]

n_probes <- nrow(tt)
cat("N probes:", n_probes, "\n")
cat("BH FDR<0.05:", sum(tt$adj.P.Val < 0.05), "\n")
cat("BH FDR<0.01:", sum(tt$adj.P.Val < 0.01), "\n")
cat("Bonferroni:", sum(tt$P.Value < 0.05/n_probes), "\n")
cat("Positive (methylation increases with grade):",
    sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), "\n")
cat("Negative:", sum(tt$adj.P.Val < 0.05 & tt$logFC < 0), "\n")

lambda <- median(qchisq(1 - tt$P.Value, df = 1)) / qchisq(0.5, df = 1)
cat("lambda_GC:", round(lambda, 3), "\n")

write.csv(tt[, c("cpg","chr","pos","gene","feat",
                 "logFC","AveExpr","t","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "toptable_q3.csv"), row.names = FALSE)
saveRDS(list(fit = fit, tt = tt, lambda = lambda),
        file.path(OUT_DIR, "q3_fit.rds"))

# ---- Diagnostic plots -------------------------------------------------------
pvals <- sort(tt$P.Value); n_p <- length(pvals)
qq_df <- data.frame(expected = -log10(ppoints(n_p)),
                    observed = -log10(pvals))
p_qq <- ggplot(qq_df, aes(expected, observed)) +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "-log10 expected p", y = "-log10 observed p",
       title = sprintf("Q3 Q-Q plot   lambda=%.3f", lambda)) +
  theme_bw()
ggsave(file.path(OUT_DIR, "qq_q3.png"), p_qq, width = 5.5, height = 5.5, dpi = 150)

tt$sig <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.2, "FDR+|effect|",
          ifelse(tt$adj.P.Val < 0.05, "FDR<0.05", "n.s."))
p_volc <- ggplot(tt, aes(logFC, -log10(P.Value), color = sig)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = c("n.s." = "grey70", "FDR<0.05" = "dodgerblue",
                                 "FDR+|effect|" = "firebrick")) +
  labs(x = "logFC (M per unit of grade)", y = "-log10(p)",
       title = "Q3 volcano (tumor-only, grade trend)") +
  theme_bw()
ggsave(file.path(OUT_DIR, "volcano_q3.png"), p_volc, width = 7, height = 5.5, dpi = 150)

# ---- Sensitivity: grade as categorical factor (ANOVA F-test) ---------------
cat("\n=== Sensitivity: grade as factor (F-test across all 4 levels) ===\n")
meta_t$grade_f <- factor(meta_t$grade_num, levels = 1:4,
                        labels = c("G1","G2","G3","G4"))
design_f <- model.matrix(~ grade_f + gender_f + age_c + tss_f, data = meta_t)
is_singular <- apply(design_f, 2, function(x) length(unique(x))) == 1
design_f <- design_f[, !is_singular, drop = FALSE]
fit_f <- lmFit(M_t, design_f); fit_f <- eBayes(fit_f)
# Test any of the three grade contrasts (G2, G3, G4 each vs G1)
grade_coefs <- grep("^grade_f", colnames(design_f), value = TRUE)
cat("Grade factor coefs:", paste(grade_coefs, collapse=", "), "\n")
tt_f <- topTable(fit_f, coef = grade_coefs, number = Inf, sort.by = "F")
tt_f$cpg <- rownames(tt_f)
cat("Factor-model FDR<0.05:", sum(tt_f$adj.P.Val < 0.05), "\n")
# Overlap with ordinal top hits
top_ord <- rownames(tt)[tt$adj.P.Val < 0.05]
top_fac <- rownames(tt_f)[tt_f$adj.P.Val < 0.05]
cat("Overlap ordinal & factor FDR<0.05:", length(intersect(top_ord, top_fac)),
    " (ord:", length(top_ord), ", fac:", length(top_fac), ")\n")
write.csv(tt_f[, c("cpg","F","P.Value","adj.P.Val")],
          file.path(OUT_DIR, "toptable_q3_factorGrade.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "sessionInfo_q3.txt"))
print(sessionInfo())
sink()
cat("\nQ3 done. Outputs in", OUT_DIR, "\n")
