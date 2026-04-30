#!/usr/bin/env Rscript
# 18_dv_dm_figure.R
#
# Two-panel figure for the DV section of the report.
# Panel A: DM logFC (M-scale) vs DV log2(var_tumor / var_normal), colored by
#   silencing/activation/uncoupled. Visualizes DV and DM as orthogonal but
#   correlated axes; cohort-wide silencing concentrates in upper-right
#   (DM up + DV up); subgroup-driven silencing tilts toward upper-middle
#   (small DM, but tumor variance > normal); RASSF1 anomalies in upper-left
#   (DM down + DV up).
# Panel B: pct_responder vs |t_DV|, demonstrating the 0.65 Spearman corr
#   underlying the "responder fraction is a thresholded non-parametric DV
#   instantiation" claim.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
DV_DIR  <- file.path(WORKDIR, "dv_out")
FIGDIR  <- file.path(WORKDIR, "figures_final")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

dt <- fread(file.path(DV_DIR, "silencing_table_with_dv.csv"))
cat("Pairs:", nrow(dt), "\n")

# Class label
dt[, class := "uncoupled"]
dt[silencing == TRUE & pattern == "subgroup",    class := "silencing (subgroup)"]
dt[silencing == TRUE & pattern == "cohort_wide", class := "silencing (cohort-wide)"]
dt[activation == TRUE,                            class := "activation"]
dt[, class := factor(class,
   levels = c("uncoupled","activation","silencing (subgroup)","silencing (cohort-wide)"))]

# Anchor genes for panel A — chosen for diagnostic value across the four
# DM x DV quadrants. Best CpG per gene by DV |t| (so the labelled point is
# the gene's strongest variability evidence, which is the figure's subject).
ANCHOR <- c("VHL","SFRP1","DKK3","GATA5",
            "RASSF1","NEFH","PITX2","CDKN2A","APAF1")
labels_dt <- dt[gene %in% ANCHOR][order(-abs(dv_t)), .SD[1], by = gene]

cls_pal <- c("uncoupled" = "grey80",
             "activation" = "#3B82F6",
             "silencing (subgroup)" = "#F59E0B",
             "silencing (cohort-wide)" = "#DC2626")

# ---- Panel A: DM vs DV scatter ---------------------------------------------
# Subsample uncoupled to avoid choking ggplot
set.seed(8139)
dt_unc <- dt[class == "uncoupled"]
n_keep <- min(50000, nrow(dt_unc))
dt_unc_sub <- dt_unc[sample.int(nrow(dt_unc), n_keep)]
dt_plot <- rbind(dt_unc_sub, dt[class != "uncoupled"])
dt_plot[, class := factor(class,
   levels = c("uncoupled","activation","silencing (subgroup)","silencing (cohort-wide)"))]

# Cap extreme values for plotting
dt_plot[, q1_logFC_cap   := pmin(pmax(q1_logFC, -3), 3)]
dt_plot[, dv_logRatio_cap := pmin(pmax(dv_logRatio, -4), 8)]

pA <- ggplot(dt_plot, aes(x = q1_logFC_cap, y = dv_logRatio_cap, color = class)) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.3) +
  geom_point(data = dt_plot[class == "uncoupled"], size = 0.18, alpha = 0.20) +
  geom_point(data = dt_plot[class == "activation"], size = 0.35, alpha = 0.55) +
  geom_point(data = dt_plot[class == "silencing (subgroup)"], size = 0.35, alpha = 0.55) +
  geom_point(data = dt_plot[class == "silencing (cohort-wide)"], size = 0.35, alpha = 0.65) +
  geom_point(data = labels_dt,
             aes(x = pmin(pmax(q1_logFC, -3), 3),
                 y = pmin(pmax(dv_logRatio, -4), 8)),
             color = "black", size = 1.5, shape = 21, fill = "white", stroke = 0.7,
             inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = labels_dt,
             aes(x = pmin(pmax(q1_logFC, -3), 3),
                 y = pmin(pmax(dv_logRatio, -4), 8),
                 label = gene),
             color = "black", size = 2.6, fontface = "italic",
             min.segment.length = 0, segment.size = 0.2, box.padding = 0.25,
             max.overlaps = 30, inherit.aes = FALSE) +
  scale_color_manual(values = cls_pal, name = NULL) +
  scale_x_continuous(breaks = c(-3, -1.5, 0, 1.5, 3),
                     labels = c("≤-3","-1.5","0","1.5","≥3")) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8),
                     labels = c("≤-4","-2","0","2","4","6","≥8")) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "DM log fold-change (M-scale, tumor − normal)",
       y = "DV log2(var tumor / var normal)",
       title = "A. DM vs DV at promoter CpGs",
       subtitle = "labelled CpGs = best per-gene CpG of 18-biomarker benchmark") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7.5),
        legend.margin = margin(t = -2),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5, color = "grey30"))

# ---- Panel B: pct_responder vs |t_DV| --------------------------------------
ok <- !is.na(dt$dv_t) & !is.na(dt$pct_responder)
rho_all <- cor(dt$pct_responder[ok], abs(dt$dv_t[ok]), method = "spearman")

dt_plot_B <- dt[ok]
dt_plot_B[, class := factor(class,
   levels = c("uncoupled","activation","silencing (subgroup)","silencing (cohort-wide)"))]
dt_plot_B[, dv_t_abs := pmin(abs(dv_t), 30)]   # cap for plotting
set.seed(8139)
dt_plot_B_unc <- dt_plot_B[class == "uncoupled"]
keepB <- min(50000, nrow(dt_plot_B_unc))
dt_plot_B_unc <- dt_plot_B_unc[sample.int(nrow(dt_plot_B_unc), keepB)]
dt_plot_B_show <- rbind(dt_plot_B_unc, dt_plot_B[class != "uncoupled"])

pB <- ggplot(dt_plot_B_show, aes(pct_responder, dv_t_abs, color = class)) +
  geom_point(data = dt_plot_B_show[class == "uncoupled"], size = 0.18, alpha = 0.18) +
  geom_point(data = dt_plot_B_show[class == "activation"], size = 0.30, alpha = 0.5) +
  geom_point(data = dt_plot_B_show[class == "silencing (subgroup)"], size = 0.30, alpha = 0.5) +
  geom_point(data = dt_plot_B_show[class == "silencing (cohort-wide)"], size = 0.30, alpha = 0.6) +
  scale_color_manual(values = cls_pal, name = NULL, guide = "none") +
  annotate("text", x = 0.02, y = 28, hjust = 0, size = 2.8,
           label = sprintf("Spearman rho = %.2f (n = %s)",
                           rho_all, format(sum(ok), big.mark = ","))) +
  scale_y_continuous(breaks = c(0, 10, 20, 30),
                     labels = c("0","10","20","≥30")) +
  labs(x = expression(paste("Responder fraction ", hat(pi)[resp])),
       y = expression(paste("|", italic(t)[DV], "|  (varFit moderated DV statistic)")),
       title = expression(bold("B. ") * bold("Responder fraction vs DV |t|"))) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 7.5, color = "grey30"))

# ---- Combine ----------------------------------------------------------------
out <- pA + pB + plot_layout(widths = c(1.15, 1))
ggsave(file.path(FIGDIR, "fig07_dv_dm.png"), out,
       width = 10.5, height = 4.6, dpi = 200, bg = "white")
cat("Wrote", file.path(FIGDIR, "fig07_dv_dm.png"), "\n")
