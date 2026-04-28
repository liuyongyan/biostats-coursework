#!/usr/bin/env Rscript
# 20_figure_main.R
# Re-render journal-quality figures for the silencing-focused final report.
# Outputs to project_workspace/figures_final/

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
Q1_DIR  <- file.path(WORKDIR, "q1_out")
Q5_DIR  <- file.path(WORKDIR, "q5_out")
SIL_DIR <- file.path(WORKDIR, "silencing_out")
EDA_DIR <- file.path(WORKDIR, "eda_outputs_kirc")
OUT_DIR <- file.path(WORKDIR, "figures_final")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(8139)

theme_paper <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_text(size = base - 1, color = "grey25"),
      axis.title = element_text(size = base),
      legend.text = element_text(size = base - 1),
      legend.title = element_text(size = base - 1, face = "bold")
    )
}

DPI <- 220

# === Figure 1: PCA ============================================================
cat("[Fig 1] PCA two-panel (sample type / sex)...\n")
meta_pc <- fread(file.path(EDA_DIR, "meta_with_pc1_10.csv"))
# Exclude the single -05A "Additional New Primary" aliquot to match analysis cohort
meta_pc <- meta_pc[sample_type_cd != "05"]
# PC1 / PC2 variance is fixed by the upstream PCA — read from existing PNG title or recompute
# Variance explained values from EDA log (02_eda_kirc.log / 03_pc_diagnostic.log)
ve1 <- 32.5; ve2 <- 19.6

p_pca_st <- ggplot(meta_pc, aes(PC1, PC2, color = sample_type)) +
  geom_point(alpha = 0.75, size = 1.4) +
  scale_color_manual(values = c("Primary Tumor" = "#1b9e77",
                                 "Solid Tissue Normal" = "#7570b3"),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", ve1),
       y = sprintf("PC2 (%.1f%%)", ve2),
       subtitle = expression("PC1" %~~% "tumor vs normal")) +
  theme_paper(9) +
  theme(legend.position = "bottom")

# Sex panel — drop samples without recorded gender for a clean comparison
meta_sex <- meta_pc[!is.na(gender) & gender %in% c("male","female")]
p_pca_sex <- ggplot(meta_sex, aes(PC1, PC2, color = gender)) +
  geom_point(alpha = 0.75, size = 1.4) +
  scale_color_manual(values = c("female" = "#d95f02", "male" = "#1f78b4"),
                     name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", ve1),
       y = sprintf("PC2 (%.1f%%)", ve2),
       subtitle = expression("PC2" %~~% "sex")) +
  theme_paper(9) +
  theme(legend.position = "bottom")

p_pca <- p_pca_st + p_pca_sex +
  plot_annotation(title = "TCGA-KIRC 450K methylation — PCA on top 10k variance CpGs",
                  subtitle = sprintf("n = %d aliquots after exclusion of the -05A 'Additional New Primary'", nrow(meta_pc)),
                  theme = theme(plot.title = element_text(size = 10, face = "bold"),
                                plot.subtitle = element_text(size = 9, color = "grey25")))

ggsave(file.path(OUT_DIR, "fig01_pca.png"), p_pca,
       width = 9, height = 4.2, dpi = DPI)

# === Figure 2: superseded by 21_figure_dm_diag.R (p-hist + Δβ-hist + volcano) =
# The Manhattan/QQ/M-scale-volcano version below is retained for reference but
# disabled — running it overwrites fig02_q1_diag.png with the inferior layout.
if (FALSE) {
cat("[Fig 2] Q1 paired DM diagnostic...\n")
q1 <- fread(file.path(Q1_DIR, "toptable_q1.csv"))
q1[, chr_num := as.numeric(sub("chr", "", chr))]
# Drop chr X/Y/M just for Manhattan (already excluded but be safe)
q1m <- q1[!is.na(chr_num) & chr_num <= 22]
q1m[, neglog10p := -log10(P.Value)]
# Subsample non-significant for plot file size (Manhattan)
q1m_sig <- q1m[adj.P.Val < 0.05]
q1m_ns  <- q1m[adj.P.Val >= 0.05][sample(.N, min(.N, 50000))]
q1m_plot <- rbind(q1m_sig, q1m_ns)
# Compute cumulative position for Manhattan-style x-axis
chr_lens <- q1m_plot[, .(maxpos = max(pos)), by = chr_num][order(chr_num)]
chr_lens[, offset := c(0, cumsum(as.numeric(maxpos))[-.N])]
q1m_plot <- merge(q1m_plot, chr_lens[, .(chr_num, offset)], by = "chr_num")
q1m_plot[, cumpos := pos + offset]
chr_centers <- q1m_plot[, .(center = mean(cumpos)), by = chr_num][order(chr_num)]

# Manhattan
p_man <- ggplot(q1m_plot, aes(x = cumpos, y = neglog10p,
                               color = factor(chr_num %% 2))) +
  geom_point(size = 0.18, alpha = 0.55) +
  geom_hline(yintercept = -log10(0.05 / nrow(q1)), linetype = "dashed",
             color = "red", linewidth = 0.3) +
  scale_color_manual(values = c("0" = "grey25", "1" = "grey55"), guide = "none") +
  scale_x_continuous(breaks = chr_centers$center, labels = chr_centers$chr_num,
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "chromosome", y = expression(-log[10](italic(p))),
       title = "Q1 paired tumor–normal DM: genome-wide p-values",
       subtitle = sprintf(
         "%s tested CpGs; dashed line = Bonferroni (p < 0.05/%s); 222{,}386 FDR<0.05",
         format(nrow(q1), big.mark = ","), format(nrow(q1), big.mark = ","))) +
  theme_paper(9) +
  theme(plot.subtitle = element_text(size = 8))

# QQ — recompute observed/expected on the full table
qq <- q1[order(P.Value)]
qq[, expected := -log10(ppoints(.N))]
qq[, observed := -log10(P.Value)]
qq_sub <- rbind(qq[1:50000], qq[sample(.N, 30000)])
lambda_q1 <- median(qchisq(1 - q1$P.Value, df = 1)) / qchisq(0.5, df = 1)
p_qq <- ggplot(qq_sub, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.3) +
  geom_point(size = 0.4, alpha = 0.6, color = "grey25") +
  labs(x = expression(expected ~ -log[10](italic(p))),
       y = expression(observed ~ -log[10](italic(p))),
       title = "Q-Q plot",
       subtitle = sprintf("λ_GC = %.1f (perm null λ ≈ 1.0)", lambda_q1)) +
  theme_paper(9)

# Volcano (β-scale)
vol <- q1[!is.na(logFC)]
vol[, dbeta := (2^logFC - 1) / (2^logFC + 1) * sign(logFC) * 0]  # placeholder
# Actually compute from M-value mean differences -> Δβ approx via mean β. Use |logFC| > 0.5849 = 2-fold M
# But we want delta-beta on β scale. logFC here is on M-scale for Q1 — convert
# M-difference to approximate Δβ assuming small effects: Δβ ≈ Δm * 0.25 (rough)
# Better: use our own β-scale Δβ if available; otherwise display logFC directly.
# Display M-scale logFC; relabel axis accordingly.
vol[, color_class := fcase(
  adj.P.Val < 0.05 & abs(logFC) > 1, "strong",
  adj.P.Val < 0.05, "FDR only",
  default = "n.s."
)]
# Subsample
ix_keep <- c(which(vol$color_class != "n.s."),
             sample(which(vol$color_class == "n.s."),
                    min(20000, sum(vol$color_class == "n.s."))))
vol_plot <- vol[ix_keep]
p_vol <- ggplot(vol_plot, aes(x = logFC, y = -log10(P.Value), color = color_class)) +
  geom_point(size = 0.25, alpha = 0.55) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(strong = "#d73027", `FDR only` = "#4575b4",
                                n.s. = "grey70"),
                     name = NULL,
                     breaks = c("strong", "FDR only", "n.s.")) +
  labs(x = expression("M-value logFC (tumor − normal)"),
       y = expression(-log[10](italic(p))),
       title = "Volcano plot",
       subtitle = sprintf(
         "red: |logFC|>1 & FDR<0.05 (n=%s); blue: FDR-only",
         format(sum(vol$color_class == "strong"), big.mark = ","))) +
  theme_paper(9) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0))

p_diag <- (p_man) / (p_qq | p_vol) + plot_layout(heights = c(1, 1.4))
ggsave(file.path(OUT_DIR, "fig02_q1_diag.png"), p_diag,
       width = 9, height = 6.5, dpi = DPI)
}  # end of disabled Fig 2 block — see 21_figure_dm_diag.R for the live version

# === Figure 3: methylation-expression coupling overview =======================
cat("[Fig 3] silencing volcano + pct_responder histogram...\n")
out <- fread(file.path(SIL_DIR, "silencing_table.csv"))
plot_dt <- copy(out)
plot_dt[, color_class := fcase(
  silencing == TRUE,  "Silencing (hyper -> mRNA down)",
  activation == TRUE, "Activation (hypo -> mRNA up)",
  default = "Not coupled"
)]
neither_ix <- which(plot_dt$color_class == "Not coupled")
if (length(neither_ix) > 30000) {
  drop_ix <- sample(neither_ix, length(neither_ix) - 30000)
  plot_dt <- plot_dt[-drop_ix]
}
# Layered approach: hexbin background for "Not coupled", scatter on top for the
# two coupled classes. This handles the dense overplotting in the central region
# while keeping the silencing/activation tails visually crisp.
plot_nc  <- plot_dt[color_class == "Not coupled"]
plot_act <- plot_dt[color_class == "Activation (hypo -> mRNA up)"]
plot_sil <- plot_dt[color_class == "Silencing (hyper -> mRNA down)"]

p_volc <- ggplot() +
  geom_hex(data = plot_nc,
           aes(x = beta_coef, y = -log10(assoc_p)),
           bins = 90, fill = "grey60", color = NA, alpha = 0.55) +
  geom_point(data = plot_act,
             aes(x = beta_coef, y = -log10(assoc_p),
                 color = color_class),
             size = 0.5, alpha = 0.75) +
  geom_point(data = plot_sil,
             aes(x = beta_coef, y = -log10(assoc_p),
                 color = color_class),
             size = 0.5, alpha = 0.75) +
  scale_color_manual(
    values = c("Silencing (hyper -> mRNA down)" = "#d73027",
               "Activation (hypo -> mRNA up)" = "#4575b4"),
    name = NULL,
    breaks = c("Silencing (hyper -> mRNA down)",
               "Activation (hypo -> mRNA up)")
  ) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(x = expression("β–mRNA coefficient (log"[2]*"(TPM+1) per unit β)"),
       y = expression(-log[10](italic(p))),
       title = "Methylation–expression coupling",
       subtitle = sprintf(
         "%s pairs; silencing %s, activation %s; grey hexbin = not coupled",
         format(nrow(out), big.mark = ","),
         format(sum(out$silencing), big.mark = ","),
         format(sum(out$activation), big.mark = ","))) +
  theme_paper(10) +
  theme(legend.position = "bottom")

sil_dt <- out[silencing == TRUE]
p_hist <- ggplot(sil_dt, aes(x = pct_responder, fill = pattern)) +
  geom_histogram(bins = 50, alpha = 0.85, boundary = 0) +
  geom_vline(xintercept = 0.3, linetype = "dashed",
             color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(subgroup = "#fc8d62", cohort_wide = "#66c2a5"),
                    name = NULL,
                    breaks = c("subgroup", "cohort_wide"),
                    labels = c("Subgroup (pct < 0.3)", "Cohort-wide (pct ≥ 0.3)")) +
  annotate("text", x = 0.31, y = Inf, hjust = 0, vjust = 1.3,
           label = "0.30", size = 3, color = "black") +
  labs(x = expression("pct_responder = #" * "{|"*beta[t] - bar(beta)[n]*"|>0.2}/" * n[tumor]),
       y = "# silencing pairs",
       title = "Responder fraction among silencing CpGs",
       subtitle = sprintf("Subgroup %s | Cohort-wide %s",
                          format(sum(sil_dt$pattern == "subgroup"), big.mark = ","),
                          format(sum(sil_dt$pattern == "cohort_wide"), big.mark = ","))) +
  theme_paper(10) + theme(legend.position = "bottom")

p_fig3 <- p_volc + p_hist + plot_layout(widths = c(1.2, 1))
ggsave(file.path(OUT_DIR, "fig03_silencing_overview.png"), p_fig3,
       width = 11, height = 4.8, dpi = DPI)

# === Figure 4: multi-gene methylation-expression scatter panel ================
cat("[Fig 4] multi-gene scatter panel...\n")
# Need β and log2tpm matrices loaded
M <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
M2beta <- function(m) 2^m / (2^m + 1)
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
rm(beta_meth); gc(verbose = FALSE)

# 6 representative genes:
#   VHL  (subgroup silencing, the case study target)
#   SFRP1 (cohort-wide silencing, Wnt antagonist)
#   DKK3 (cohort-wide silencing, Wnt antagonist)
#   GATA5 (subgroup silencing, ccRCC TSG)
#   PCDH17 (cohort-wide silencing, neural cadherin)
#   RASSF1 (negative case: hypomethylation in tumor; classic literature TSG that
#           does NOT show silencing direction in TCGA-KIRC)
focal_genes <- c("VHL","SFRP1","DKK3","GATA5","PCDH17","RASSF1")
make_scatter <- function(g) {
  rows <- out[gene == g & cpg %in% rownames(beta_mat)]
  if (nrow(rows) == 0) return(NULL)
  sil_rows <- rows[silencing == TRUE]
  best <- if (nrow(sil_rows) > 0) sil_rows[which.min(assoc_p)] else rows[which.min(assoc_p)]
  cpg <- best$cpg
  df <- data.frame(beta = as.numeric(beta_mat[cpg, ]),
                   y    = as.numeric(log2tpm[g, ]),
                   sample_type = match_meta$sample_type)
  df_t <- df[df$sample_type == "Primary Tumor", ]
  df_n <- df[df$sample_type == "Solid Tissue Normal", ]
  rho <- cor(df$beta, df$y, method = "spearman")
  status <- fcase(best$silencing == TRUE, "silencing",
                  best$activation == TRUE, "activation",
                  default = "not coupled")
  ann_color <- c("silencing" = "#d73027", "activation" = "#4575b4",
                 "not coupled" = "grey40")[status]
  # Place stat label below the title, outside the data area, so it never
  # overlaps points. Tumour points first, normals on top with a larger
  # triangle so they're visible despite small n.
  ggplot() +
    geom_point(data = df_t, aes(x = beta, y = y),
               color = "#d73027", size = 1.2, alpha = 0.55) +
    geom_point(data = df_n, aes(x = beta, y = y),
               color = "#1f78b4", size = 2.6, alpha = 0.95, shape = 17) +
    labs(x = bquote(.(cpg) ~ beta ~ "(" * .(best$promoter_group) * ")"),
         y = bquote(log[2]("TPM" + 1)),
         title = bquote(italic(.(g))),
         subtitle = sprintf("ρ=%+.2f, p=%.1g, π_resp=%.2f — %s",
                            rho, best$assoc_p, best$pct_responder, status)) +
    theme_paper(9) +
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          plot.subtitle = element_text(color = ann_color, size = 7.5))
}
plots <- lapply(focal_genes, make_scatter)
p_panel <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    subtitle = "Red circles = primary tumor (n=318); blue triangles = adjacent normal (n=24).",
    theme = theme(plot.subtitle = element_text(size = 9, color = "grey25"))
  )
ggsave(file.path(OUT_DIR, "fig04_scatter_panel.png"), p_panel,
       width = 10, height = 6.5, dpi = DPI)

# === Figure 5: VHL/cg13672843 mechanistic case ================================
cat("[Fig 5] VHL mechanistic case (copy existing)...\n")
file.copy(file.path(WORKDIR, "vhl_expr_out", "vhl_subgroup_scatter_single.png"),
          file.path(OUT_DIR, "fig05_vhl_case.png"),
          overwrite = TRUE)

# === Figure 6: pathway comparison Q1 vs silencing-restricted ==================
cat("[Fig 6] pathway Q1 vs silencing comparison...\n")
q1_gsea <- fread(file.path(Q5_DIR, "gsea_q1.csv"))
setnames(q1_gsea, c("pathway","pval","padj","NES","size","ES"))
ora_sil <- fread(file.path(SIL_DIR, "ora_silencing.csv"))
setnames(ora_sil, "name", "pathway")

# Top 10 from raw-DM GSEA
top_q1 <- head(q1_gsea[order(padj)], 10)
top_q1[, source := "Raw DM (GSEA)"]
# Top 10 from silencing ORA
top_sil <- head(ora_sil[order(padj)], 10)
top_sil[, source := "Silencing-restricted (ORA)"]
top_sil[, NES := NA_real_]  # ORA has no NES
# Build a unified plot
plot_q1  <- top_q1[, .(pathway, padj, source, signed_log = -log10(padj) * sign(NES))]
plot_sil <- top_sil[, .(pathway, padj, source, signed_log = -log10(padj))]
plot_combo <- rbind(plot_q1, plot_sil)
plot_combo[, pathway_clean := sub("^HALLMARK_", "", pathway)]
plot_combo[, pathway_clean := sub("^KEGG_", "KEGG ", pathway_clean)]
plot_combo[, pathway_clean := gsub("_", " ", pathway_clean)]
# Mark the inflammatory-response pathway in the raw-DM panel — used in §3.5
# as the canonical example of an immune signal that disappears under silencing
plot_combo[, highlight_immune := pathway %in%
              c("HALLMARK_INFLAMMATORY_RESPONSE",
                "HALLMARK_COMPLEMENT",
                "HALLMARK_ALLOGRAFT_REJECTION",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "KEGG_HEMATOPOIETIC_CELL_LINEAGE",
                "HALLMARK_IL6_JAK_STAT3_SIGNALING")]

p_pathway <- ggplot(plot_combo,
                     aes(x = -log10(padj),
                         y = reorder(pathway_clean, -log10(padj)),
                         fill = source,
                         alpha = highlight_immune)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  annotate("text", x = -log10(0.05) + 0.05, y = 0.6,
           label = "FDR = 0.05", hjust = 0, vjust = 0,
           size = 2.5, color = "grey30") +
  facet_wrap(~ source, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Raw DM (GSEA)" = "#7570b3",
                                "Silencing-restricted (ORA)" = "#1b9e77"),
                    guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 0.95),
                     guide = "none") +
  labs(x = expression(-log[10]("FDR-adjusted p")),
       y = NULL,
       title = "Top pathway enrichment: raw DM (GSEA) vs silencing-restricted (ORA)",
       subtitle = "Saturated bars in the left panel = immune/composition pathways that drop to non-significance on the right.") +
  theme_paper(9) +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8, color = "grey25"))
ggsave(file.path(OUT_DIR, "fig06_pathway_compare.png"), p_pathway,
       width = 11, height = 5, dpi = DPI)

# === Supplementary fig: immune dropout direct comparison ======================
# For paper, useful as inset / supplementary table. Build a lookup.
immune_paths <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_COMPLEMENT",
                   "HALLMARK_ALLOGRAFT_REJECTION",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                   "KEGG_HEMATOPOIETIC_CELL_LINEAGE")
imm <- merge(
  q1_gsea[pathway %in% immune_paths,
          .(pathway, q1_padj = padj, q1_NES = NES)],
  ora_sil[pathway %in% immune_paths,
          .(pathway, sil_padj = padj, sil_p = p, k, K)],
  by = "pathway")
imm <- imm[order(q1_padj)]
fwrite(imm, file.path(OUT_DIR, "supp_immune_dropout.csv"))
cat("Immune dropout table:\n")
print(imm)

cat("\nDone. Figures in", OUT_DIR, "\n")
list.files(OUT_DIR)
