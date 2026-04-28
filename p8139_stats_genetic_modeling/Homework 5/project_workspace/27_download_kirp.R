#!/usr/bin/env Rscript
# 27_download_kirp.R
# Download TCGA-KIRP 450K methylation + STAR-Counts RNA-seq + clinical.
# Mirrors 01_download_kirc.R but for project = TCGA-KIRP.
# Idempotent: skips already-downloaded data.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

set.seed(8139)
WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DATA_DIR <- file.path(WORKDIR, "data_kirp")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()
cat(sprintf("=== Start: %s ===\n", t_start))

# ---- Methylation 450K -------------------------------------------------------
cat("\n[1/3] Querying KIRP 450K methylation...\n")
query_meth <- GDCquery(
  project = "TCGA-KIRP",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
meth_files <- getResults(query_meth)
cat(sprintf("  450K files: %d  |  unique cases: %d\n",
            nrow(meth_files), length(unique(meth_files$cases.submitter_id))))
cat("  Sample types:\n")
print(table(meth_files$sample_type))
saveRDS(query_meth, file.path(DATA_DIR, "query_meth_450k.rds"))
write.csv(meth_files, file.path(DATA_DIR, "query_meth_450k_manifest.csv"),
          row.names = FALSE)

cat("  Downloading methylation files...\n")
GDCdownload(query = query_meth, method = "api", files.per.chunk = 50,
            directory = file.path(DATA_DIR, "GDCdata"))

cat("  Preparing methylation SE...\n")
meth_se <- GDCprepare(query_meth, directory = file.path(DATA_DIR, "GDCdata"))
saveRDS(meth_se, file.path(DATA_DIR, "kirp_450k_SE.rds"))
cat(sprintf("  Methylation SE: %d CpGs × %d samples\n",
            nrow(meth_se), ncol(meth_se)))

# ---- RNA-seq STAR-Counts -----------------------------------------------------
cat("\n[2/3] Querying KIRP RNA-seq STAR-Counts...\n")
query_rna <- GDCquery(
  project = "TCGA-KIRP",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
rna_files <- getResults(query_rna)
cat(sprintf("  RNA-seq files: %d  |  unique cases: %d\n",
            nrow(rna_files), length(unique(rna_files$cases.submitter_id))))
saveRDS(query_rna, file.path(DATA_DIR, "query_rna.rds"))
cat("  Downloading RNA-seq files...\n")
GDCdownload(query = query_rna, method = "api", files.per.chunk = 50,
            directory = file.path(DATA_DIR, "GDCdata"))
cat("  Preparing RNA-seq SE...\n")
rna_se <- GDCprepare(query_rna, directory = file.path(DATA_DIR, "GDCdata"))
saveRDS(rna_se, file.path(DATA_DIR, "kirp_rnaseq_SE.rds"))
cat(sprintf("  RNA-seq SE: %d genes × %d samples\n",
            nrow(rna_se), ncol(rna_se)))

# ---- Clinical ----------------------------------------------------------------
cat("\n[3/3] Querying KIRP clinical metadata...\n")
clinical <- GDCquery_clinic(project = "TCGA-KIRP", type = "clinical")
saveRDS(clinical, file.path(DATA_DIR, "clinical.rds"))
write.csv(clinical, file.path(DATA_DIR, "clinical.csv"), row.names = FALSE)
cat(sprintf("  Clinical: %d patients × %d cols\n", nrow(clinical), ncol(clinical)))

t_end <- Sys.time()
cat(sprintf("\n=== Done: %s (elapsed %.1f min) ===\n",
            t_end, as.numeric(difftime(t_end, t_start, units = "mins"))))
print(list.files(DATA_DIR))

sink(file.path(DATA_DIR, "sessionInfo_download.txt"))
print(sessionInfo())
sink()
