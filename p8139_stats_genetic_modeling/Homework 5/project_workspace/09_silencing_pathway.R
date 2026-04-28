#!/usr/bin/env Rscript
# 09_silencing_pathway.R
# Downstream of 07_methylation_silencing.R:
#   - Pattern partition figure (pct_responder histogram with 0.3 cut)
#   - Silencing volcano (beta_coef vs -log10 assoc_p, colored by silencing/activation)
#   - Multi-gene methylation-expression scatter panel (VHL, SFRP1, RASSF1, GATA5, DKK3, PCDH17)
#   - GSEA on silencing-strict gene list (gene-level: best CpG per gene)
#   - Comparison table: silencing-restricted GSEA vs raw Q1 GSEA

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR  <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
SIL_DIR  <- file.path(WORKDIR, "silencing_out")
Q5_DIR   <- file.path(WORKDIR, "q5_out")
OUT_DIR  <- SIL_DIR
set.seed(8139)

theme_paper <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = base + 1, face = "bold"),
          plot.subtitle = element_text(size = base - 1, color = "grey30"))
}

# === Load silencing table =====================================================
cat("Loading silencing_table.csv...\n")
out <- fread(file.path(SIL_DIR, "silencing_table.csv"))
cat(sprintf("  %d (CpG, gene) pairs\n", nrow(out)))
cat(sprintf("  Silencing: %d (in %d unique genes)\n",
            sum(out$silencing), uniqueN(out$gene[out$silencing])))

# === Load β / log2tpm for scatter panel =======================================
cat("Loading β / log2tpm matrices for scatter panel...\n")
M <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
M2beta <- function(m) 2^m / (2^m + 1)

# Aggregate methylation to (patient, sample_type)
meth_keys_all <- paste(meta$patient, meta$sample_type, sep = "|")
beta_full <- M2beta(M); rm(M); gc(verbose = FALSE)
unique_meth_keys <- unique(meth_keys_all)
beta_meth <- if (length(unique_meth_keys) == ncol(beta_full)) {
  colnames(beta_full) <- meth_keys_all; beta_full
} else {
  agg <- vapply(unique_meth_keys, function(k) {
    ix <- which(meth_keys_all == k)
    if (length(ix) == 1) beta_full[, ix] else rowMeans(beta_full[, ix, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(beta_full)))
  colnames(agg) <- unique_meth_keys; agg
}
rm(beta_full); gc(verbose = FALSE)
meta_unique <- meta[!duplicated(meth_keys_all), ]
rownames(meta_unique) <- paste(meta_unique$patient, meta_unique$sample_type, sep = "|")

# RNA-seq
rna_se <- readRDS(file.path(DATA_DIR, "kirc_rnaseq_SE.rds"))
rna_cd <- as.data.frame(colData(rna_se))[, c("barcode","patient","sample_type")]
keep <- rna_cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
rna_se <- rna_se[, keep]; rna_cd <- rna_cd[keep, , drop = FALSE]
gene_names <- as.character(rowData(rna_se)$gene_name)
tpm <- assay(rna_se, "tpm_unstrand")
rna_keys <- paste(rna_cd$patient, rna_cd$sample_type, sep = "|")
unique_keys <- unique(rna_keys)
tpm_avg <- vapply(unique_keys, function(k) {
  ix <- which(rna_keys == k)
  if (length(ix) == 1) tpm[, ix] else rowMeans(tpm[, ix, drop = FALSE])
}, numeric(nrow(tpm)))
colnames(tpm_avg) <- unique_keys
keep_g <- !is.na(gene_names) & gene_names != ""
tpm_avg <- tpm_avg[keep_g, , drop = FALSE]
gene_lab <- gene_names[keep_g]
log2tpm <- log2(rowsum(tpm_avg, group = gene_lab, reorder = FALSE) + 1)

matched_keys <- intersect(unique_meth_keys, colnames(log2tpm))
match_meta <- meta_unique[matched_keys, , drop = FALSE]
beta_mat <- beta_meth[, matched_keys]
log2tpm  <- log2tpm[, matched_keys]
cat(sprintf("  Matched: %d samples (%d tumor, %d normal)\n",
            length(matched_keys),
            sum(match_meta$sample_type == "Primary Tumor"),
            sum(match_meta$sample_type == "Solid Tissue Normal")))

# === Figure A: silencing volcano ==============================================
cat("[Figure A] silencing volcano...\n")
plot_dt <- copy(out)
plot_dt[, color_class := fcase(
  silencing == TRUE,  "silencing",
  activation == TRUE, "activation",
  default = "neither"
)]
# Subsample 'neither' for plot file size
set.seed(1)
neither_ix <- which(plot_dt$color_class == "neither")
if (length(neither_ix) > 30000) {
  drop_ix <- sample(neither_ix, length(neither_ix) - 30000)
  plot_dt <- plot_dt[-drop_ix]
}

p_volcano <- ggplot(plot_dt, aes(x = beta_coef, y = -log10(assoc_p),
                                  color = color_class, alpha = color_class,
                                  size = color_class)) +
  geom_point() +
  scale_color_manual(values = c(silencing = "#d73027",
                                 activation = "#4575b4",
                                 neither = "grey70"),
                     labels = c(silencing = "Silencing (hyper -> mRNA down)",
                                activation = "Activation (hypo -> mRNA up)",
                                neither = "Not coupled"),
                     name = NULL) +
  scale_alpha_manual(values = c(silencing = 0.7, activation = 0.7, neither = 0.25),
                     guide = "none") +
  scale_size_manual(values = c(silencing = 0.5, activation = 0.5, neither = 0.25),
                    guide = "none") +
  labs(x = expression("β-mRNA association coefficient"),
       y = expression(-log[10](italic(p))),
       title = "Methylation–expression coupling at promoter CpGs",
       subtitle = sprintf("%d (CpG, gene) pairs; %d silencing (red), %d activation (blue)",
                          nrow(out), sum(out$silencing), sum(out$activation))) +
  theme_paper(10) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "silencing_volcano.png"), p_volcano,
       width = 7, height = 5, dpi = 200)

# === Figure B: pct_responder histogram (silencing CpGs only) ==================
cat("[Figure B] pct_responder histogram for silencing CpGs...\n")
sil_dt <- out[silencing == TRUE]
p_hist <- ggplot(sil_dt, aes(x = pct_responder, fill = pattern)) +
  geom_histogram(bins = 40, alpha = 0.85) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(subgroup = "#fc8d62", cohort_wide = "#66c2a5"),
                     name = NULL) +
  annotate("text", x = 0.31, y = Inf, hjust = 0, vjust = 1.5,
           label = "threshold = 0.30", size = 3, color = "black") +
  labs(x = expression("pct_responder = fraction of tumors with " * "|" * beta[tumor] - bar(beta)[normal] * "|" > 0.2),
       y = "# silencing (CpG, gene) pairs",
       title = "Distribution of responder fraction among silencing CpGs",
       subtitle = sprintf("Subgroup: %d pairs (pct_resp < 0.3) | Cohort-wide: %d pairs (pct_resp >= 0.3)",
                          sum(sil_dt$pattern == "subgroup"),
                          sum(sil_dt$pattern == "cohort_wide"))) +
  theme_paper(10) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "silencing_pct_responder_hist.png"), p_hist,
       width = 7, height = 4, dpi = 200)

# === Figure C: multi-gene methylation-expression scatter panel ================
cat("[Figure C] multi-gene scatter panel...\n")
focal_genes <- c("VHL", "SFRP1", "DKK3", "GATA5", "PCDH17", "RASSF1")
make_scatter <- function(g) {
  rows <- out[gene == g & cpg %in% rownames(beta_mat)]
  if (nrow(rows) == 0) return(NULL)
  # pick best-p silencing CpG; if none silencing, best-p assoc CpG
  sil_rows <- rows[silencing == TRUE]
  best <- if (nrow(sil_rows) > 0) sil_rows[which.min(assoc_p)]
          else rows[which.min(assoc_p)]
  cpg <- best$cpg
  beta_x <- beta_mat[cpg, ]
  log2tpm_y <- log2tpm[g, ]
  df <- data.frame(beta = as.numeric(beta_x),
                   y    = as.numeric(log2tpm_y),
                   sample_type = match_meta$sample_type)
  rho <- cor(df$beta, df$y, method = "spearman")
  status <- if (best$silencing) "silencing"
            else if (best$activation) "activation" else "not coupled"
  ggplot(df, aes(x = beta, y = y, color = sample_type)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_manual(values = c("Primary Tumor" = "#d73027",
                                   "Solid Tissue Normal" = "#4575b4"),
                       guide = "none") +
    labs(x = sprintf("%s β", cpg),
         y = bquote(log[2]("TPM"+1)),
         title = bquote(italic(.(g)) ~ " (" * .(best$promoter_group) * ")"),
         subtitle = sprintf("ρ = %+.2f; assoc p = %.1e; pct_resp = %.2f; %s",
                            rho, best$assoc_p, best$pct_responder, status)) +
    theme_paper(8.5)
}
plots <- lapply(focal_genes, make_scatter)
p_panel <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = "Promoter methylation × gene expression at representative loci",
    subtitle = paste0("Red = tumor; blue = adjacent normal. ",
                      "VHL/SFRP1/DKK3/GATA5/PCDH17 = silencing; RASSF1 = ccRCC-specific non-silencing"),
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )
ggsave(file.path(OUT_DIR, "silencing_scatter_panel.png"), p_panel,
       width = 9, height = 6, dpi = 200)

rm(beta_mat, log2tpm, beta_meth); gc(verbose = FALSE)

# === GSEA on silencing-restricted gene list ===================================
cat("Building gene-level rank statistic for silencing GSEA...\n")
# Continuous rank over all tested genes:
#   for each gene, take the best CpG by assoc_p (most-significant β-mRNA coupling),
#   rank_stat = sign(beta_coef) * -log10(assoc_p)
# Convention: negative rank = silencing direction (high β -> low mRNA);
#             positive rank = activation direction (high β -> high mRNA).
# Ties are broken naturally by magnitude, no zero-pile.
gene_best <- out[, .SD[which.min(assoc_p)], by = gene,
                 .SDcols = c("assoc_p", "beta_coef", "q1_logFC")]
gene_best[, rank_stat := sign(beta_coef) * -log10(assoc_p)]
ranks_silencing <- setNames(gene_best$rank_stat, gene_best$gene)
ranks_silencing <- ranks_silencing[is.finite(ranks_silencing)]
cat(sprintf("Ranked %d genes (continuous; sign(beta_coef) * -log10 best assoc p per gene)\n",
            length(ranks_silencing)))
cat(sprintf("  rank range: [%.1f, %.1f], median %.2f\n",
            min(ranks_silencing), max(ranks_silencing), median(ranks_silencing)))

# Build pathway lists (Hallmark + KEGG_LEGACY)
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg_leg <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:KEGG_LEGACY")
pathways <- split(c(hallmark$gene_symbol, kegg_leg$gene_symbol),
                   c(hallmark$gs_name,    kegg_leg$gs_name))

# fgsea
fg <- fgsea(pathways = pathways, stats = ranks_silencing,
            minSize = 10, maxSize = 500, nproc = 1)
fg <- fg[order(padj)]
fg_dt <- as.data.table(fg)
fg_dt[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"),
                               character(1))]
fwrite(fg_dt, file.path(OUT_DIR, "gsea_silencing.csv"))
cat("Top 15 silencing-restricted GSEA hits:\n")
print(head(fg_dt[, .(pathway, NES, pval, padj, size)], 15))

# Also load q5_out/gsea_q1.csv for comparison table
q1_gsea <- fread(file.path(Q5_DIR, "gsea_q1.csv"))
setnames(q1_gsea, c("pathway","pval","padj","NES","size","ES"))
q1_gsea <- q1_gsea[order(padj)]
top_q1 <- head(q1_gsea[, .(pathway, NES = round(NES, 2),
                            padj = signif(padj, 2),
                            size)], 15)
top_sil <- head(fg_dt[, .(pathway, NES = round(NES, 2),
                          padj = signif(padj, 2),
                          size)], 15)
fwrite(top_q1, file.path(OUT_DIR, "gsea_compare_top15_q1.csv"))
fwrite(top_sil, file.path(OUT_DIR, "gsea_compare_top15_silencing.csv"))

# Comparison: which top-15 Q1 hits drop out of the silencing-restricted ranking?
cat("\n=== GSEA comparison: Q1 top 15 vs silencing top 15 ===\n")
cat("\n[Q1 raw DM] top 15:\n"); print(top_q1)
cat("\n[Silencing-restricted] top 15:\n"); print(top_sil)

# Specifically inspect the immune pathways
immune_paths <- grep("INFLAMMATORY|COMPLEMENT|ALLOGRAFT|INTERFERON|IL[0-9]|IMMUNE",
                     fg_dt$pathway, value = TRUE)
cat("\nImmune-related pathways in silencing GSEA (padj):\n")
print(fg_dt[pathway %in% immune_paths,
            .(pathway, NES = round(NES, 2), padj = signif(padj, 2))][order(padj)])

# Specifically check KEGG_RENAL_CELL_CARCINOMA
cat("\nKEGG_RENAL_CELL_CARCINOMA in silencing GSEA:\n")
print(fg_dt[pathway == "KEGG_RENAL_CELL_CARCINOMA",
            .(pathway, NES = round(NES, 2), pval = signif(pval, 2),
              padj = signif(padj, 2))])
cat("KEGG_RENAL_CELL_CARCINOMA in Q1 GSEA:\n")
print(q1_gsea[pathway == "KEGG_RENAL_CELL_CARCINOMA",
              .(pathway, NES = round(NES, 2), pval = signif(pval, 2),
                padj = signif(padj, 2))])

# === ORA on silencing-strict gene list ========================================
# Hypergeometric test: foreground = silencing-strict genes; background = all tested genes
cat("\n[ORA] silencing-strict foreground vs all-tested background...\n")
fg_genes <- unique(out$gene[out$silencing])
bg_genes <- unique(out$gene)
N_bg <- length(bg_genes)
N_fg <- length(fg_genes)
cat(sprintf("  Foreground silencing genes: %d\n  Background tested genes: %d\n", N_fg, N_bg))

ora_test <- function(set_genes, set_name) {
  set_in_bg <- intersect(set_genes, bg_genes)
  if (length(set_in_bg) < 5) return(NULL)
  k <- length(intersect(set_in_bg, fg_genes))
  K <- length(set_in_bg)
  data.table(pathway = set_name,
             k = k, K = K,
             foreground_size = N_fg,
             background_size = N_bg,
             p = phyper(k - 1, K, N_bg - K, N_fg, lower.tail = FALSE))
}

pathway_dt <- rbind(
  data.table(name = hallmark$gs_name, gene = hallmark$gene_symbol),
  data.table(name = kegg_leg$gs_name, gene = kegg_leg$gene_symbol)
)
ora_results <- pathway_dt[, ora_test(unique(gene), name[1]), by = name]
ora_results[, padj := p.adjust(p, method = "BH")]
setorder(ora_results, p)
fwrite(ora_results, file.path(OUT_DIR, "ora_silencing.csv"))

cat("\n[ORA] Top 20 silencing pathways:\n")
print(head(ora_results[, .(pathway = name, k, K, p = signif(p, 2),
                            padj = signif(padj, 2))], 20))

cat("\n[ORA] KEGG_RENAL_CELL_CARCINOMA:\n")
print(ora_results[name == "KEGG_RENAL_CELL_CARCINOMA",
                  .(pathway = name, k, K, p = signif(p, 2),
                    padj = signif(padj, 2))])

cat("\n[ORA] Immune-related sets (sanity: should NOT be enriched in silencing FG):\n")
immune_sets <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_COMPLEMENT",
                 "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_IL2_STAT5_SIGNALING",
                 "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
print(ora_results[name %in% immune_sets,
                  .(pathway = name, k, K, p = signif(p, 2),
                    padj = signif(padj, 2))])

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_20.txt"))
cat("\nDone. Outputs in", OUT_DIR, "\n")
