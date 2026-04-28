# Download TCGA-KIRC 450K methylation + clinical data.
#
# Strategy:
#   1. Query all open-access 450K methylation beta-value files for TCGA-KIRC.
#   2. Query clinical + biospecimen metadata.
#   3. Download beta files via TCGAbiolinks::GDCdownload (method = "api").
#   4. Save the GDCquery objects as RDS so later steps can resume without
#      re-hitting the API.
#
# This script is idempotent: re-running skips already-downloaded files.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DATA_DIR <- file.path(WORKDIR, "data_kirc")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Step 1: query methylation files ===\n")
query_meth <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
meth_files <- getResults(query_meth)
cat("450K methylation files:", nrow(meth_files), "\n")
cat("Unique cases:", length(unique(meth_files$cases.submitter_id)), "\n")
cat("Sample types:\n")
print(table(meth_files$sample_type))

saveRDS(query_meth, file.path(DATA_DIR, "query_meth_450k.rds"))
write.csv(meth_files, file.path(DATA_DIR, "query_meth_450k_manifest.csv"), row.names = FALSE)

cat("\n=== Step 2: query clinical ===\n")
# GDCquery_clinic pulls the harmonized clinical table.
clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical")
cat("Clinical rows:", nrow(clinical), "  cols:", ncol(clinical), "\n")
saveRDS(clinical, file.path(DATA_DIR, "clinical.rds"))
write.csv(clinical, file.path(DATA_DIR, "clinical.csv"), row.names = FALSE)

cat("\n=== Step 3: download methylation files ===\n")
GDCdownload(
  query = query_meth,
  method = "api",
  files.per.chunk = 50,
  directory = file.path(DATA_DIR, "GDCdata")
)

cat("\n=== Done ===\n")
cat("DATA_DIR:", DATA_DIR, "\n")
print(list.files(DATA_DIR))

sink(file.path(DATA_DIR, "sessionInfo_download.txt"))
print(sessionInfo())
sink()
