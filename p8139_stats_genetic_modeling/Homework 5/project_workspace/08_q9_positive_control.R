# Step 3e — Q9 positive control: VHL/HIF pipeline sanity check.
#
# Pre-registered assertions (plan.md §5):
#   A1. At least one VHL promoter CpG (probes mapping to VHL TSS +/- 1500 bp)
#       appears hypermethylated in tumor vs normal with FDR < 0.05.
#   A2. HALLMARK_HYPOXIA and/or KEGG_RENAL_CELL_CARCINOMA enriches in Q5 GSEA
#       (FDR < 0.25).
#   A3. VHL promoter methylation direction is positive (tumor > normal).
#
# If any assertion fails, PAUSE — audit the pipeline.

suppressPackageStartupMessages({
  library(dplyr)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
OUT_DIR <- file.path(WORKDIR, "q9_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

tt_q1 <- read.csv(file.path(WORKDIR, "q1_out", "toptable_q1.csv"),
                  stringsAsFactors = FALSE)

# ---- A1 / A3: VHL promoter check --------------------------------------------
# UCSC_RefGene_Group flags the feature region; promoter-relevant:
#   TSS1500, TSS200, 5'UTR, 1stExon
promoter_flags <- c("TSS1500", "TSS200", "5'UTR", "1stExon")

has_vhl <- function(gene_str) {
  if (is.na(gene_str)) return(FALSE)
  any(toupper(unlist(strsplit(gene_str, ";"))) == "VHL")
}
has_promoter <- function(feat_str) {
  if (is.na(feat_str) || feat_str == "") return(FALSE)
  any(unlist(strsplit(feat_str, ";")) %in% promoter_flags)
}

tt_q1$is_vhl <- vapply(tt_q1$gene, has_vhl, logical(1))
tt_q1$is_vhl_promoter <- tt_q1$is_vhl & vapply(tt_q1$feat, has_promoter, logical(1))

vhl_rows <- tt_q1[tt_q1$is_vhl, ]
vhl_prom <- tt_q1[tt_q1$is_vhl_promoter, ]
cat("VHL-associated probes total:", nrow(vhl_rows), "\n")
cat("VHL promoter probes:", nrow(vhl_prom), "\n")
if (nrow(vhl_prom) > 0) {
  print(vhl_prom[, c("cpg","chr","pos","feat","logFC","P.Value","adj.P.Val")])
}
write.csv(vhl_prom, file.path(OUT_DIR, "vhl_promoter_probes.csv"), row.names = FALSE)
write.csv(vhl_rows, file.path(OUT_DIR, "vhl_all_probes.csv"), row.names = FALSE)

# Assertion A1: any VHL promoter CpG with FDR<0.05 and logFC>0
A1_pass <- any(vhl_prom$adj.P.Val < 0.05 & vhl_prom$logFC > 0)
# Assertion A3: among VHL promoter CpGs with FDR<0.05, median logFC > 0
sig_vhl <- vhl_prom[vhl_prom$adj.P.Val < 0.05, ]
A3_pass <- if (nrow(sig_vhl) > 0) median(sig_vhl$logFC) > 0 else FALSE

cat("\nA1 (VHL promoter CpG FDR<0.05 & logFC>0):", A1_pass, "\n")
cat("A3 (VHL promoter median logFC positive):", A3_pass, "\n")

# ---- A2: Pathway enrichment check -------------------------------------------
gsea_q1_path <- file.path(WORKDIR, "q5_out", "gsea_q1.csv")
if (file.exists(gsea_q1_path)) {
  gsea_q1 <- read.csv(gsea_q1_path, stringsAsFactors = FALSE)
  target_sets <- c("HALLMARK_HYPOXIA", "KEGG_RENAL_CELL_CARCINOMA")
  hits <- gsea_q1[gsea_q1$pathway %in% target_sets, ]
  cat("\nTarget GSEA hits (Q1):\n"); print(hits)
  A2_pass <- any(hits$padj < 0.25, na.rm = TRUE)
} else {
  cat("GSEA Q1 results not found — A2 cannot be evaluated\n")
  A2_pass <- NA
}
cat("A2 (HYPOXIA or RENAL_CELL_CARCINOMA padj<0.25):", A2_pass, "\n")

# ---- Summary report ---------------------------------------------------------
result <- data.frame(
  assertion = c("A1_vhl_promoter_sig", "A2_pathway_enrichment", "A3_vhl_direction"),
  pass = c(A1_pass, A2_pass, A3_pass),
  note = c(
    sprintf("%d/%d VHL promoter CpGs FDR<0.05 & logFC>0",
            sum(vhl_prom$adj.P.Val < 0.05 & vhl_prom$logFC > 0), nrow(vhl_prom)),
    "HYPOXIA or RCC in Q1 GSEA padj<0.25",
    if (nrow(sig_vhl) > 0)
      sprintf("median logFC across FDR<0.05 VHL promoter = %.3f", median(sig_vhl$logFC))
    else "no FDR-sig VHL promoter CpG"
  )
)
print(result)
write.csv(result, file.path(OUT_DIR, "q9_assertions.csv"), row.names = FALSE)

overall <- all(result$pass, na.rm = TRUE) && !any(is.na(result$pass))
cat(sprintf("\n=== Overall positive-control verdict: %s ===\n",
            if (isTRUE(overall)) "PASS" else "FAIL (investigate)"))
cat("\nQ9 done. Outputs in", OUT_DIR, "\n")
