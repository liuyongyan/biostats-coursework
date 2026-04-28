# Step 3a — Preprocessing for Q1 / Q3.
# Pre-registered pipeline from plan.md §5:
#   1. Drop chr X / chr Y probes (sex dominates PC2 per EDA)
#   2. Drop cross-reactive probes (Chen 2013 list via DMRcate::rmSNPandCH)
#   3. Drop SNP-overlap probes (same call)
#   4. Drop probes with >20% sample missingness
#   5. Keep complete-case probes
#   6. beta -> M = log2(beta/(1-beta)) with clipping [0.001, 0.999]
# Also build sample metadata table and drop the 1 Additional New Primary sample.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(dplyr)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DATA_DIR <- file.path(WORKDIR, "data_kirc")
OUT_DIR  <- file.path(WORKDIR, "preprocess_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Load cached SE + clinical ===\n")
se       <- readRDS(file.path(DATA_DIR, "kirc_450k_SE.rds"))
clinical <- readRDS(file.path(DATA_DIR, "clinical.rds"))
cat("SE: ", nrow(se), " probes x ", ncol(se), " samples\n")

# ---- Build sample metadata ---------------------------------------------------
bc <- colnames(se)
meta <- data.frame(
  barcode        = bc,
  patient        = substr(bc, 1, 12),
  tss            = substr(bc, 6, 7),
  sample_type_cd = substr(bc, 14, 15),
  plate          = substr(bc, 22, 25),
  center         = substr(bc, 27, 28),
  stringsAsFactors = FALSE
)
meta$sample_type <- dplyr::recode(meta$sample_type_cd,
  "01" = "Primary Tumor", "05" = "Additional New Primary", "11" = "Solid Tissue Normal",
  .default = paste0("Other_", meta$sample_type_cd))

clin_keep <- c("submitter_id", "gender", "race", "age_at_diagnosis",
               "ajcc_pathologic_stage", "ajcc_pathologic_t",
               "tumor_grade", "vital_status", "primary_diagnosis")
clin_keep <- intersect(clin_keep, colnames(clinical))
meta <- dplyr::left_join(meta, clinical[, clin_keep], by = c("patient" = "submitter_id"))

cat("\nSample type breakdown:\n"); print(table(meta$sample_type))

# Drop the 1 Additional New Primary (pre-registered decision)
keep_sample <- meta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
cat("Dropping", sum(!keep_sample), "samples (non primary/normal):\n")
print(meta$barcode[!keep_sample])
meta <- meta[keep_sample, ]
se   <- se[, keep_sample]
stopifnot(all(colnames(se) == meta$barcode))
cat("After sample filter:", ncol(se), "samples\n")

# ---- Probe filter: chr X / Y --------------------------------------------------
cat("\n=== Probe filter: chr X/Y ===\n")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Intersect probe IDs (SE set may differ from annotation set by a few probes)
common_probes <- intersect(rownames(se), rownames(anno))
cat("SE probes:", nrow(se), "  Annotation probes:", nrow(anno),
    "  Intersection:", length(common_probes), "\n")
if (length(common_probes) < 0.95 * nrow(se)) {
  stop("Too few common probes between SE and annotation — check platform mismatch")
}
se <- se[common_probes, ]
anno <- anno[common_probes, ]
stopifnot(identical(rownames(anno), rownames(se)))

is_sex <- anno$chr %in% c("chrX", "chrY")
cat("probes on chrX:", sum(anno$chr == "chrX"), "\n")
cat("probes on chrY:", sum(anno$chr == "chrY"), "\n")
se_ <- se[!is_sex, ]
anno_ <- anno[!is_sex, ]
cat("after chrX/Y drop:", nrow(se_), "\n")

# ---- Probe filter: cross-reactive + SNP (DMRcate) ----------------------------
cat("\n=== Probe filter: cross-reactive + SNP overlap (DMRcate) ===\n")
if (requireNamespace("DMRcate", quietly = TRUE)) {
  beta_mat <- assay(se_)
  # rmSNPandCH expects a beta or M matrix with probe rownames
  # mafcut = 0.05 and dist = 2 are defaults; rmcrosshyb = TRUE drops cross-reactive
  kept <- DMRcate::rmSNPandCH(beta_mat, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE)
  cat("DMRcate probes kept:", nrow(kept), " (dropped:", nrow(beta_mat) - nrow(kept), ")\n")
  keep_probes <- rownames(kept)
  se_ <- se_[keep_probes, ]
  anno_ <- anno_[match(keep_probes, rownames(anno_)), ]
} else {
  cat("!!! DMRcate unavailable — falling back to SNP-only filter via annotation Probe_maf > 0.05\n")
  has_snp <- !is.na(anno_$Probe_rs) & !is.na(anno_$Probe_maf) & anno_$Probe_maf > 0.05
  cat("probes with MAF>0.05 SNP overlap:", sum(has_snp), "\n")
  se_ <- se_[!has_snp, ]
  anno_ <- anno_[!has_snp, ]
  cat("(cross-reactive filter not applied — documented as a limitation)\n")
}

# ---- Probe filter: missingness ----------------------------------------------
cat("\n=== Probe filter: >20% missing -> drop; then require complete case ===\n")
beta <- assay(se_)
na_per_probe <- rowSums(is.na(beta)) / ncol(beta)
cat("probes with >20% missing:", sum(na_per_probe > 0.20), "\n")
se_ <- se_[na_per_probe <= 0.20, ]
anno_ <- anno_[na_per_probe <= 0.20, ]
cat("after >20% missing drop:", nrow(se_), "\n")

beta <- assay(se_)
complete_rows <- rowSums(is.na(beta)) == 0
cat("probes with zero missing:", sum(complete_rows), "\n")
se_ <- se_[complete_rows, ]
anno_ <- anno_[complete_rows, ]
cat("final probe count:", nrow(se_), "\n")

# ---- beta -> M --------------------------------------------------------------
cat("\n=== beta -> M transform with clipping ===\n")
beta <- assay(se_)
beta_clip <- pmin(pmax(beta, 0.001), 0.999)
M <- log2(beta_clip / (1 - beta_clip))
cat("M matrix:", nrow(M), "x", ncol(M), "  any NA:", any(is.na(M)), "\n")

# ---- Save -------------------------------------------------------------------
stopifnot(identical(colnames(M), meta$barcode))
saveRDS(M,     file.path(OUT_DIR, "M_matrix.rds"))
saveRDS(anno_, file.path(OUT_DIR, "probe_annotation.rds"))
saveRDS(meta,  file.path(OUT_DIR, "sample_metadata.rds"))

cat("\nWritten to", OUT_DIR, ":\n")
print(list.files(OUT_DIR))

sink(file.path(OUT_DIR, "sessionInfo_preprocess.txt"))
print(sessionInfo())
sink()
cat("\nDone.\n")
