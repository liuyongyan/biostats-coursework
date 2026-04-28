#!/usr/bin/env Rscript
# 22_refigure_3_4.R
# Redo Figs 3 and 4 with cleaner aesthetics.
# - Fig 3: winsorize volcano x-axis, lighter "not coupled" points.
# - Fig 4: stats as in-panel annotation (no subtitle clipping); larger facets;
#          β = 0.2 reference line; readable point sizes.

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
SIL_DIR <- file.path(WORKDIR, "silencing_out")
OUT_DIR <- file.path(WORKDIR, "figures_final")
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

out <- fread(file.path(SIL_DIR, "silencing_table.csv"))

# === Figure 3 ================================================================
cat("[Fig 3] silencing volcano + pct_responder histogram...\n")
plot_dt <- copy(out)
plot_dt[, color_class := fcase(
  silencing == TRUE,  "Silencing (hyper -> mRNA down)",
  activation == TRUE, "Activation (hypo -> mRNA up)",
  default = "Not coupled"
)]
# winsorize coefficients for plot only
XLIM <- c(-12, 12)
plot_dt[, beta_coef_w := pmin(pmax(beta_coef, XLIM[1]), XLIM[2])]
neither_ix <- which(plot_dt$color_class == "Not coupled")
if (length(neither_ix) > 30000) {
  drop_ix <- sample(neither_ix, length(neither_ix) - 30000)
  plot_dt <- plot_dt[-drop_ix]
}
# Plot order: not-coupled first (background), then activation, then silencing on top
plot_dt[, color_class := factor(color_class,
  levels = c("Not coupled", "Activation (hypo -> mRNA up)",
             "Silencing (hyper -> mRNA down)"))]
plot_dt <- plot_dt[order(color_class)]

p_volc <- ggplot(plot_dt, aes(x = beta_coef_w, y = -log10(assoc_p),
                              color = color_class, alpha = color_class,
                              size = color_class)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50",
             linewidth = 0.3) +
  scale_color_manual(
    values = c("Silencing (hyper -> mRNA down)" = "#d73027",
               "Activation (hypo -> mRNA up)" = "#4575b4",
               "Not coupled" = "grey80"),
    name = NULL,
    breaks = c("Silencing (hyper -> mRNA down)",
               "Activation (hypo -> mRNA up)", "Not coupled")
  ) +
  scale_alpha_manual(values = c("Not coupled" = 0.20,
                                "Activation (hypo -> mRNA up)" = 0.75,
                                "Silencing (hyper -> mRNA down)" = 0.75),
                     guide = "none") +
  scale_size_manual(values = c("Not coupled" = 0.16,
                               "Activation (hypo -> mRNA up)" = 0.55,
                               "Silencing (hyper -> mRNA down)" = 0.55),
                    guide = "none") +
  scale_x_continuous(limits = XLIM, expand = expansion(mult = 0.02)) +
  labs(x = expression("β–mRNA coefficient (log"[2]*"(TPM+1) per unit β; winsorized to ±12)"),
       y = expression(-log[10](italic(p))),
       title = "Methylation–expression coupling",
       subtitle = sprintf(
         "%s pairs; silencing %s, activation %s",
         format(nrow(out), big.mark = ","),
         format(sum(out$silencing), big.mark = ","),
         format(sum(out$activation), big.mark = ","))) +
  theme_paper(10) + theme(legend.position = "bottom")

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
cat("  written:", file.path(OUT_DIR, "fig03_silencing_overview.png"), "\n")

# === Figure 4 ================================================================
cat("[Fig 4] multi-gene scatter panel...\n")
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
  rho <- cor(df$beta, df$y, method = "spearman", use = "complete.obs")
  status <- fcase(best$silencing == TRUE, "silencing",
                  best$activation == TRUE, "activation",
                  default = "not coupled")
  ann_color <- c("silencing" = "#b30000", "activation" = "#08519c",
                 "not coupled" = "grey25")[status]
  status_label <- c("silencing" = "silencing",
                    "activation" = "activation",
                    "not coupled" = "not coupled")[status]
  ann_text <- sprintf("ρ = %+.2f\np_assoc = %.1g\nπ_resp = %.2f\n%s",
                      rho, best$assoc_p, best$pct_responder, status_label)
  # Place annotation at top-right of plot area; expand y a bit
  ymax <- max(df$y, na.rm = TRUE)
  ymin <- min(df$y, na.rm = TRUE)
  ypad <- 0.18 * (ymax - ymin)
  xrng <- range(df$beta, na.rm = TRUE)
  ggplot(df, aes(x = beta, y = y, color = sample_type)) +
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "grey50",
               linewidth = 0.3) +
    geom_point(size = 1.0, alpha = 0.55) +
    scale_color_manual(values = c("Primary Tumor" = "#d73027",
                                   "Solid Tissue Normal" = "#1f78b4"),
                       guide = "none") +
    annotate("label", x = xrng[2], y = ymax + ypad,
             label = ann_text,
             hjust = 1, vjust = 1, size = 2.7, label.size = 0,
             fill = scales::alpha("white", 0.85),
             color = ann_color, lineheight = 0.95) +
    coord_cartesian(ylim = c(ymin - 0.05 * (ymax - ymin),
                             ymax + ypad), clip = "off") +
    labs(x = bquote(.(cpg) ~ beta ~ "(" * .(best$promoter_group) * ")"),
         y = bquote(log[2]("TPM" + 1)),
         title = bquote(italic(.(g)))) +
    theme_paper(9) +
    theme(plot.title = element_text(size = 11, face = "bold.italic", hjust = 0),
          plot.margin = margin(6, 8, 4, 6))
}
plots <- lapply(focal_genes, make_scatter)
p_panel <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = "Promoter methylation × gene expression at representative loci",
    subtitle = "Red = primary tumor (n = 318); blue = adjacent normal (n = 24). Dotted line at β = 0.2 (responder threshold).",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9, color = "grey25"))
  )
ggsave(file.path(OUT_DIR, "fig04_scatter_panel.png"), p_panel,
       width = 10, height = 6.8, dpi = DPI)
cat("  written:", file.path(OUT_DIR, "fig04_scatter_panel.png"), "\n")

cat("\nDone.\n")
