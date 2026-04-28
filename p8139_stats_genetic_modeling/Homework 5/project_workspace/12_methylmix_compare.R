#!/usr/bin/env Rscript
# 12_methylmix_compare.R
# Continue from where 11_methylmix_benchmark.R crashed: load saved MethylMix
# result and our silencing list, compute overlap statistics.

suppressPackageStartupMessages({
  library(data.table)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
SIL_DIR <- file.path(WORKDIR, "silencing_out")
OUT_DIR <- file.path(WORKDIR, "methylmix_out")

# Load MethylMix result
mm <- readRDS(file.path(OUT_DIR, "methylmix_result.rds"))
mm_inputs <- readRDS(file.path(OUT_DIR, "methylmix_inputs.rds"))

mm_genes <- as.character(mm$MethylationDrivers)
cat(sprintf("MethylationDrivers: %d genes\n", length(mm_genes)))
cat(sprintf("First 20: %s\n", paste(head(mm_genes, 20), collapse = ", ")))
writeLines(mm_genes, file.path(OUT_DIR, "methylmix_genes.txt"))

# Load our silencing list
sil_tab <- fread(file.path(SIL_DIR, "silencing_table.csv"))
our_silencing <- unique(sil_tab$gene[sil_tab$silencing == TRUE])
N_back <- nrow(mm_inputs$METcancer)  # background = MethylMix's tested set
N_ours <- length(our_silencing)
N_mm   <- length(mm_genes)

# Restrict our silencing to the genes tested by MethylMix (fair comparison)
mm_tested <- rownames(mm_inputs$METcancer)
our_tested <- intersect(our_silencing, mm_tested)
N_ours_in_mmbg <- length(our_tested)

N_both <- length(intersect(our_tested, mm_genes))
N_only_ours <- length(setdiff(our_tested, mm_genes))
N_only_mm   <- length(setdiff(mm_genes, our_tested))

cat(sprintf("\n=== Background (MethylMix-tested genes): %d ===\n", N_back))
cat(sprintf("Our silencing (intersect MethylMix bg): %d\n", N_ours_in_mmbg))
cat(sprintf("MethylMix MethylationDrivers: %d\n", N_mm))
cat(sprintf("Overlap: %d  | Only ours: %d  | Only MethylMix: %d\n",
            N_both, N_only_ours, N_only_mm))

jaccard <- N_both / length(union(our_tested, mm_genes))
cat(sprintf("Jaccard: %.3f\n", jaccard))

# Hypergeometric on overlap given ours as foreground
hyp_p <- phyper(N_both - 1, N_mm, N_back - N_mm, N_ours_in_mmbg, lower.tail = FALSE)
cat(sprintf("Hypergeometric p (overlap > random): %.2e\n", hyp_p))

# Recovery of pre-registered TSGs
TSG18 <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
           "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
           "TIMP3","NEFH","UCHL1","PITX2")
mm_tsg   <- intersect(TSG18, mm_genes)
ours_tsg <- intersect(TSG18, our_silencing)
cat(sprintf("\n=== 18-TSG recovery ===\n"))
cat(sprintf("ours    : %2d/18  %s\n", length(ours_tsg),
            paste(ours_tsg, collapse = ", ")))
cat(sprintf("MethylMix: %2d/18  %s\n", length(mm_tsg),
            paste(mm_tsg, collapse = ", ")))
cat(sprintf("intersection: %2d  %s\n",
            length(intersect(ours_tsg, mm_tsg)),
            paste(intersect(ours_tsg, mm_tsg), collapse = ", ")))
cat(sprintf("only-ours : %2d  %s\n",
            length(setdiff(ours_tsg, mm_tsg)),
            paste(setdiff(ours_tsg, mm_tsg), collapse = ", ")))
cat(sprintf("only-mm   : %2d  %s\n",
            length(setdiff(mm_tsg, ours_tsg)),
            paste(setdiff(mm_tsg, ours_tsg), collapse = ", ")))

# Write a summary CSV
summary_dt <- data.frame(
  metric = c("background_genes_tested_by_methylmix",
             "ours_silencing_in_mm_bg", "mm_drivers", "overlap",
             "ours_only", "mm_only",
             "jaccard", "hypergeometric_p",
             "ours_TSG_recovered", "mm_TSG_recovered", "both_TSG_recovered"),
  value  = c(N_back,
             N_ours_in_mmbg, N_mm, N_both,
             N_only_ours, N_only_mm,
             round(jaccard, 4), formatC(hyp_p, format = "e", digits = 2),
             length(ours_tsg), length(mm_tsg),
             length(intersect(mm_tsg, ours_tsg))))
fwrite(summary_dt, file.path(OUT_DIR, "benchmark_summary.csv"))
cat(sprintf("\nWritten: %s\n", file.path(OUT_DIR, "benchmark_summary.csv")))
