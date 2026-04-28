# Install the R environment for the TCGA methylation project.
# Run once. Writes sessionInfo() to sessionInfo_after_install.txt on success.

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "minfi",
  "sesame",
  "sesameData",
  "limma",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "missMethyl"
)

cran_pkgs <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "data.table", "pheatmap", "qqman"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

cat("\n--- sessionInfo ---\n")
sink("sessionInfo_after_install.txt")
print(sessionInfo())
sink()
cat("Done. Wrote sessionInfo_after_install.txt\n")
