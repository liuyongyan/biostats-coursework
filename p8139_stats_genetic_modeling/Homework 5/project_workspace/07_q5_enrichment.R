# Step 3d — Q5 interpretation: GSEA + ORA on Q1 and Q3 CpG lists.
#
# Plan:
#   1. Map each CpG to a single gene (take first name in UCSC_RefGene_Name,
#      skipping CpGs with no gene).
#   2. Per gene, keep the CpG with the largest |t| (within each of Q1 and Q3).
#   3. Rank genes by signed -log10(p) * sign(logFC) for GSEA.
#   4. Run fgsea against MSigDB Hallmark (50 sets) + KEGG renal carcinoma set.
#   5. Run ORA (hypergeometric) on FDR<0.05 gene list vs all tested genes.

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

set.seed(8139)

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
OUT_DIR <- file.path(WORKDIR, "q5_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

tt_q1 <- read.csv(file.path(WORKDIR, "q1_out", "toptable_q1.csv"), stringsAsFactors = FALSE)
tt_q3 <- read.csv(file.path(WORKDIR, "q3_out", "toptable_q3.csv"), stringsAsFactors = FALSE)

first_gene <- function(s) {
  s[is.na(s) | s == ""] <- NA
  sapply(strsplit(s, ";"), function(x) if (length(x) == 0) NA_character_ else x[1])
}
tt_q1$gene1 <- first_gene(tt_q1$gene)
tt_q3$gene1 <- first_gene(tt_q3$gene)

collapse_by_gene <- function(tt) {
  tt2 <- tt[!is.na(tt$gene1), ]
  tt2$rank_metric <- -log10(pmax(tt2$P.Value, 1e-300)) * sign(tt2$logFC)
  tt2 %>%
    group_by(gene1) %>%
    arrange(desc(abs(rank_metric))) %>%
    slice(1) %>% ungroup() %>%
    as.data.frame()
}
g_q1 <- collapse_by_gene(tt_q1)
g_q3 <- collapse_by_gene(tt_q3)
cat("Unique genes Q1:", nrow(g_q1), "  Q3:", nrow(g_q3), "\n")

# ---- Gene sets --------------------------------------------------------------
cat("\n=== Build gene set collections (MSigDB Hallmark + KEGG renal) ===\n")
# msigdbr API can be either msigdbr() or msigdbr_df()
if ("msigdbr" %in% getNamespaceExports("msigdbr")) {
  hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
} else {
  hallmark <- msigdbr::msigdbr_df(species = "Homo sapiens", collection = "H")
}
if ("gs_name" %in% colnames(hallmark)) {
  nm_col <- "gs_name"; sym_col <- "gene_symbol"
} else if ("gs_name" %in% colnames(hallmark)) {
  nm_col <- "gs_name"; sym_col <- "gene_symbol"
} else {
  # newer API uses "gs_name" and "gene_symbol" still
  nm_col <- grep("^gs_name$|name$", colnames(hallmark), value = TRUE)[1]
  sym_col <- grep("symbol", colnames(hallmark), value = TRUE)[1]
}
hallmark_list <- split(hallmark[[sym_col]], hallmark[[nm_col]])
cat("Hallmark sets:", length(hallmark_list), "\n")

# KEGG renal cell carcinoma (C2:CP:KEGG)
kegg <- tryCatch({
  if ("msigdbr" %in% getNamespaceExports("msigdbr")) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  } else {
    msigdbr::msigdbr_df(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  }
}, error = function(e) NULL)
if (is.null(kegg) || nrow(kegg) == 0) {
  kegg <- tryCatch({
    if ("msigdbr" %in% getNamespaceExports("msigdbr")) {
      msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG")
    } else {
      msigdbr::msigdbr_df(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG")
    }
  }, error = function(e) NULL)
}
if (!is.null(kegg) && nrow(kegg) > 0) {
  renal_sets <- unique(kegg[[nm_col]][grepl("RENAL", kegg[[nm_col]], ignore.case = TRUE)])
  cat("KEGG renal sets:", paste(renal_sets, collapse = ", "), "\n")
  kegg_list <- split(kegg[[sym_col]], kegg[[nm_col]])
  # Merge
  set_list <- c(hallmark_list, kegg_list)
} else {
  cat("KEGG not available — using Hallmark only\n")
  set_list <- hallmark_list
}
cat("Total gene sets:", length(set_list), "\n")

# ---- GSEA -------------------------------------------------------------------
run_gsea <- function(df, label) {
  r <- df$rank_metric; names(r) <- df$gene1
  r <- r[!duplicated(names(r)) & !is.na(r)]
  r <- sort(r, decreasing = TRUE)
  res <- fgsea::fgsea(pathways = set_list, stats = r, minSize = 10, maxSize = 500)
  res <- res[order(pval), ]
  write.csv(as.data.frame(res[, c("pathway","pval","padj","NES","size","ES")]),
            file.path(OUT_DIR, paste0("gsea_", label, ".csv")), row.names = FALSE)
  cat("GSEA", label, ": significant padj<0.25:", sum(res$padj < 0.25, na.rm = TRUE), "\n")
  cat("Top 10 by padj:\n")
  print(head(res[, c("pathway","padj","NES","size")], 10))
  res
}

cat("\n=== GSEA on Q1 ===\n")
gsea_q1 <- run_gsea(g_q1, "q1")
cat("\n=== GSEA on Q3 ===\n")
gsea_q3 <- run_gsea(g_q3, "q3")

# ---- ORA (hypergeometric) on FDR<0.05 gene list ----------------------------
run_ora <- function(df, label) {
  sig_genes <- unique(df$gene1[df$adj.P.Val < 0.05])
  universe <- unique(df$gene1)
  cat("ORA", label, ": foreground", length(sig_genes), "/", length(universe), "\n")
  if (length(sig_genes) == 0) {
    cat("  no sig genes — skipping ORA\n")
    return(NULL)
  }
  res <- do.call(rbind, lapply(names(set_list), function(nm) {
    gs <- intersect(set_list[[nm]], universe)
    k <- length(intersect(sig_genes, gs))
    K <- length(gs); N <- length(universe); n <- length(sig_genes)
    if (K < 5) return(NULL)
    p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    data.frame(pathway = nm, k = k, K = K, n = n, N = N, p = p)
  }))
  res$padj <- p.adjust(res$p, method = "BH")
  res <- res[order(res$p), ]
  write.csv(res, file.path(OUT_DIR, paste0("ora_", label, ".csv")), row.names = FALSE)
  cat("  ORA top hits:\n"); print(head(res, 10))
  res
}

cat("\n=== ORA on Q1 ===\n")
ora_q1 <- run_ora(g_q1, "q1")
cat("\n=== ORA on Q3 ===\n")
ora_q3 <- run_ora(g_q3, "q3")

sink(file.path(OUT_DIR, "sessionInfo_q5.txt"))
print(sessionInfo())
sink()
cat("\nQ5 done. Outputs in", OUT_DIR, "\n")
