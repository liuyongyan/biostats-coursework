#!/usr/bin/env Rscript
# 24_refig2.R (v2)
# Replace the GWAS-style Manhattan + Q-Q + volcano with three diagnostics
# that actually carry information for paired tumor-vs-normal cancer methylation:
#   A. p-value histogram with permutation-null overlay
#       (visualizes calibration AND signal pervasiveness)
#   B. paired Δβ histogram with biological-significance threshold lines
#       (visualizes effect-size distribution)
#   C. volcano with paired Δβ axis and 5 anchor TSGs labelled
# Layout: (A | B) on top, C on the bottom.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

set.seed(8139)
WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
Q1_DIR  <- file.path(WORKDIR, "q1_out")
SIL_DIR <- file.path(WORKDIR, "silencing_out")
OUT_DIR <- file.path(WORKDIR, "figures_final")

theme_paper <- function(base = 9) {
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

# === Load data ===============================================================
cat("[load] Q1 toptable + permutation p-values...\n")
q1 <- fread(file.path(Q1_DIR, "toptable_q1.csv"))
perm <- readRDS(file.path(SIL_DIR, "perm_pvalues.rds"))

# === Panel A: p-value histogram + permutation overlay ========================
cat("[Panel A] p-value histogram\n")
n_obs <- nrow(q1)
lambda_obs <- median(qchisq(1 - q1$P.Value, df = 1)) / qchisq(0.5, df = 1)
lambda_perm <- median(qchisq(1 - as.vector(perm$P_perm), df = 1)) / qchisq(0.5, df = 1)

# Bin uniformly on (0,1)
breaks <- seq(0, 1, by = 0.025)
real_hist <- hist(q1$P.Value, breaks = breaks, plot = FALSE)
perm_hist <- hist(as.vector(perm$P_perm), breaks = breaks, plot = FALSE)
hist_dt <- data.table(
  mid       = real_hist$mids,
  real_dens = real_hist$density,
  perm_dens = perm_hist$density
)

# Cap real density at 12 for plotting; flag the spike
Y_CAP_HIST <- 12
hist_dt[, real_capped := real_dens > Y_CAP_HIST]
hist_dt[, real_plot   := pmin(real_dens, Y_CAP_HIST)]

p_hist <- ggplot(hist_dt) +
  geom_col(aes(x = mid, y = real_plot, fill = "Real labels"),
           width = 0.024, alpha = 0.85) +
  geom_step(aes(x = mid - 0.0125, y = perm_dens,
                color = "Permutation null (sign-flip)"),
            linewidth = 0.6, direction = "hv") +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_text(data = hist_dt[real_capped == TRUE],
            aes(x = mid, y = Y_CAP_HIST + 0.6,
                label = sprintf("%.0f", real_dens)),
            size = 2.4, color = "#b30000", vjust = 0) +
  scale_fill_manual(values = c("Real labels" = "#fcae91"), name = NULL) +
  scale_color_manual(values = c("Permutation null (sign-flip)" = "#1f78b4"),
                     name = NULL) +
  scale_y_continuous(limits = c(0, Y_CAP_HIST + 1.5), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "p-value (paired moderated t)", y = "density",
       title = "p-value distribution: real vs permutation null",
       subtitle = sprintf(
         "Real density capped at %d (text = uncapped); λ_GC: real %.1f vs perm %.2f",
         Y_CAP_HIST, lambda_obs, lambda_perm)) +
  theme_paper(9) +
  theme(legend.position = c(0.62, 0.95),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                          color = NA),
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(8, "pt"),
        plot.subtitle = element_text(size = 7.5))

# === Panel B: Δβ histogram with thresholds ===================================
cat("[Panel B] Δβ histogram (computing paired mean Δβ)\n")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
M2beta <- function(m) 2^m / (2^m + 1)
beta_full <- M2beta(M); rm(M); gc(verbose = FALSE)

tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")

delta_beta <- rowMeans(beta_full[, bc_T] - beta_full[, bc_N], na.rm = TRUE)
rm(beta_full); gc(verbose = FALSE)
delta_beta <- delta_beta[!is.na(delta_beta)]

# Counts beyond |Δβ| > 0.2 thresholds
n_total  <- length(delta_beta)
n_hyper  <- sum(delta_beta >  0.2)
n_hypo   <- sum(delta_beta < -0.2)
n_strong <- n_hyper + n_hypo

p_dbeta <- ggplot(data.frame(dbeta = delta_beta), aes(x = dbeta)) +
  geom_histogram(binwidth = 0.01, fill = "#74add1", color = NA, alpha = 0.85) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6),
                     breaks = c(-0.5, -0.2, 0, 0.2, 0.5)) +
  scale_y_sqrt(breaks = c(100, 1000, 5000, 20000, 50000, 100000),
               labels = scales::comma, expand = c(0, 0)) +
  annotate("text", x = -0.55, y = sqrt(40000), hjust = 0, vjust = 0,
           label = sprintf("%s with Δβ < −0.2", format(n_hypo, big.mark = ",")),
           size = 2.7, color = "#1f78b4") +
  annotate("text", x = 0.55, y = sqrt(40000), hjust = 1, vjust = 0,
           label = sprintf("%s with Δβ > 0.2", format(n_hyper, big.mark = ",")),
           size = 2.7, color = "#b30000") +
  labs(x = expression("paired mean " * Delta * beta * " (tumor − normal)"),
       y = "# CpGs (sqrt scale)",
       title = "Effect-size distribution",
       subtitle = sprintf(
         "%s tested CpGs; %s (%.1f%%) clear |Δβ|>0.2 biological cutoff",
         format(n_total, big.mark = ","),
         format(n_strong, big.mark = ","),
         100 * n_strong / n_total)) +
  theme_paper(9) +
  theme(plot.subtitle = element_text(size = 7.5))

# === Panel C: Volcano with anchors ===========================================
cat("[Panel C] volcano with anchors\n")
q1[, dbeta := delta_beta[match(cpg, names(delta_beta))]]
vol <- q1[!is.na(dbeta) & !is.na(P.Value)]

vol[, color_class := fcase(
  adj.P.Val < 0.05 & abs(dbeta) > 0.1, "strong",
  adj.P.Val < 0.05, "FDR only",
  default = "n.s."
)]
ix_keep <- c(which(vol$color_class != "n.s."),
             sample(which(vol$color_class == "n.s."),
                    min(20000, sum(vol$color_class == "n.s."))))
vol_plot <- vol[ix_keep]

anchors <- data.table(
  cpg  = c("cg13672843", "cg17816908", "cg21620540", "cg12664464", "cg06174454"),
  gene = c("VHL", "SFRP1", "DKK3", "GATA5", "RASSF1")
)
anch <- merge(anchors, vol[, .(cpg, dbeta, P.Value)],
              by = "cpg", all.x = TRUE)

Y_CAP_V <- 50
vol_plot[, neglogp := pmin(-log10(P.Value), Y_CAP_V)]
vol_plot[, capped  := -log10(P.Value) > Y_CAP_V]
anch[, neglogp := pmin(-log10(P.Value), Y_CAP_V)]

p_vol <- ggplot(vol_plot, aes(x = dbeta, y = neglogp, color = color_class)) +
  geom_point(size = 0.25, alpha = 0.55) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_point(data = vol_plot[capped == TRUE],
             aes(x = dbeta, y = Y_CAP_V), color = "black",
             shape = 17, size = 0.6) +
  geom_point(data = anch, aes(x = dbeta, y = neglogp),
             color = "black", size = 1.5, inherit.aes = FALSE) +
  geom_text_repel(data = anch, aes(x = dbeta, y = neglogp, label = gene),
                  inherit.aes = FALSE, size = 3, fontface = "italic",
                  min.segment.length = 0, segment.size = 0.3,
                  box.padding = 0.5, max.overlaps = Inf,
                  seed = 1) +
  scale_color_manual(values = c(strong = "#d73027", `FDR only` = "#4575b4",
                                n.s. = "grey70"),
                     name = NULL,
                     breaks = c("strong", "FDR only", "n.s."),
                     labels = c(expression("|"*Delta*beta*"| > 0.1 & FDR < 0.05"),
                                expression("FDR < 0.05 only"),
                                "n.s.")) +
  scale_x_continuous(limits = c(-0.6, 0.6),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
  scale_y_continuous(limits = c(0, Y_CAP_V * 1.02), expand = c(0, 0)) +
  labs(x = expression("paired mean " * Delta * beta * " (tumor − normal)"),
       y = expression(-log[10](italic(p))),
       title = "Joint view: significance vs effect size",
       subtitle = "y capped at 50 (▲ = beyond cap); labelled ● = anchor TSGs (5)") +
  theme_paper(9) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0),
        plot.subtitle = element_text(size = 7.5))

# === Stitch ==================================================================
p_diag <- (p_hist | p_dbeta) / p_vol +
  plot_layout(heights = c(1, 1.3))

ggsave(file.path(OUT_DIR, "fig02_q1_diag.png"), p_diag,
       width = 11, height = 6.8, dpi = DPI)

cat(sprintf("\nWritten %s\n", file.path(OUT_DIR, "fig02_q1_diag.png")))
