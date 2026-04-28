#!/usr/bin/env Rscript
# 08_silencing_fdr_grid.R
# Sensitivity analysis for the silencing screen across the (q1_FDR, assoc_FDR)
# threshold grid. Produces silencing_out/sensitivity_fdr_grid.csv, which the
# final report reads in §3.7.
#
# Method
# ------
#   For each (q1_FDR_t, assoc_FDR_t) in {0.01, 0.05, 0.10}^2:
#     1. Re-threshold silencing_table.csv (output of 07_methylation_silencing.R)
#        using the strict-silencing definition:
#            q1_FDR < q1_FDR_t  AND  assoc_FDR < assoc_FDR_t
#            AND  q1_logFC > 0  AND  beta_coef < 0
#     2. Count silencing (CpG, gene) pairs and unique silencing genes.
#     3. Run hypergeometric ORA against Hallmark + KEGG_LEGACY (same
#        pathway sources and background as 09_silencing_pathway.R), record
#        the top-5 pathways by raw p.
#
# Output columns: q1_FDR, assoc_FDR, n_pairs, n_genes, top5_ORA.
# (The default cell q1_FDR=0.05, assoc_FDR=0.05 reproduces the headline
# 12,922 pairs / 4,483 genes reported in §3.3.)

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
SIL_DIR <- file.path(WORKDIR, "silencing_out")

cat("[load] silencing_table.csv\n")
out <- fread(file.path(SIL_DIR, "silencing_table.csv"))

cat("[load] Hallmark + KEGG_LEGACY gene sets via msigdbr\n")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg_leg <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:KEGG_LEGACY")
pathway_dt <- rbind(
  data.table(name = hallmark$gs_name, gene = hallmark$gene_symbol),
  data.table(name = kegg_leg$gs_name, gene = kegg_leg$gene_symbol)
)

# Hypergeometric ORA: foreground = silencing genes for this grid cell,
# background = all tested genes (same as 09_silencing_pathway.R).
ora_top5 <- function(fg_genes, bg_genes) {
  N_bg <- length(bg_genes); N_fg <- length(fg_genes)
  if (N_fg == 0) return("")
  res <- pathway_dt[, {
    set_in_bg <- intersect(unique(gene), bg_genes)
    K <- length(set_in_bg)
    if (K < 5) .(p = NA_real_) else {
      k <- length(intersect(set_in_bg, fg_genes))
      .(p = phyper(k - 1, K, N_bg - K, N_fg, lower.tail = FALSE))
    }
  }, by = name]
  setorder(res, p)
  paste(res[!is.na(p), name][1:5], collapse = " | ")
}

bg_all <- unique(out$gene)
grid <- CJ(q1_FDR = c(0.01, 0.05, 0.10),
           assoc_FDR = c(0.01, 0.05, 0.10))

cat("[grid] running 9 cells...\n")
results <- grid[, {
  sil_idx <- !is.na(q1_FDR.x <- out$q1_FDR) & !is.na(out$assoc_FDR) &
              out$q1_FDR    < q1_FDR &
              out$assoc_FDR < assoc_FDR &
              out$q1_logFC > 0 & out$beta_coef < 0
  fg_genes <- unique(out$gene[sil_idx])
  cat(sprintf("  q1_FDR=%.2f assoc_FDR=%.2f -> %d pairs, %d genes\n",
              q1_FDR, assoc_FDR, sum(sil_idx), length(fg_genes)))
  .(n_pairs = sum(sil_idx),
    n_genes = length(fg_genes),
    top5_ORA = ora_top5(fg_genes, bg_all))
}, by = .(q1_FDR, assoc_FDR)]

fwrite(results, file.path(SIL_DIR, "sensitivity_fdr_grid.csv"))
cat(sprintf("\n[done] wrote %s\n", file.path(SIL_DIR, "sensitivity_fdr_grid.csv")))
print(results)
