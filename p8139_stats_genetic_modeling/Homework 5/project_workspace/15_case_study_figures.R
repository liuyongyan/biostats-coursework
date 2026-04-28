# Generate case-study figures for the final report:
#   - VHL 7 promoter CpGs: per-patient Δβ strip + paired tumor/normal boxplots
#     → illustrates the "subgroup-driven with mixed direction" pattern
#   - SFRP1 top promoter CpG: per-patient Δβ strip + paired boxes
#     → illustrates the "tissue-wide hypermethylation" pattern

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(SummarizedExperiment)
  library(patchwork)
})

WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
OUT_DIR <- file.path(WORKDIR, "casestudy_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load -------------------------------------------------------------------
M      <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta   <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
anno_  <- readRDS(file.path(PRE_DIR, "probe_annotation.rds"))
tt_q1  <- read.csv(file.path(WORKDIR, "q1_out", "toptable_q1.csv"),
                   stringsAsFactors = FALSE)

# Build paired matrices
M2beta <- function(m) 2^m / (2^m + 1)
tumor_pat  <- unique(meta$patient[meta$sample_type == "Primary Tumor"])
normal_pat <- unique(meta$patient[meta$sample_type == "Solid Tissue Normal"])
paired_pat <- intersect(tumor_pat, normal_pat)
pick_bc <- function(p, type) meta$barcode[meta$patient == p & meta$sample_type == type][1]
bc_T <- vapply(paired_pat, pick_bc, character(1), "Primary Tumor")
bc_N <- vapply(paired_pat, pick_bc, character(1), "Solid Tissue Normal")
beta_T <- M2beta(M[, bc_T]); colnames(beta_T) <- paired_pat
beta_N <- M2beta(M[, bc_N]); colnames(beta_N) <- paired_pat

# ---- Gene-specific CpG finder ---------------------------------------------
has_gene <- function(gene_str, target) {
  if (is.na(gene_str) || gene_str == "") return(FALSE)
  any(toupper(unlist(strsplit(gene_str, ";"))) %in% toupper(target))
}
promoter_flags <- c("TSS1500","TSS200","5'UTR","1stExon")
has_prom <- function(feat_str) {
  if (is.na(feat_str) || feat_str == "") return(FALSE)
  any(unlist(strsplit(feat_str, ";")) %in% promoter_flags)
}

# ---- Build per-gene dataframe for plotting --------------------------------
build_gene_df <- function(gene_name, target_names = NULL) {
  if (is.null(target_names)) target_names <- gene_name
  tt_q1$has_gene <- vapply(tt_q1$gene, has_gene, logical(1), target = target_names)
  tt_q1$has_prom <- vapply(tt_q1$feat, has_prom, logical(1))
  gene_cpgs <- tt_q1[tt_q1$has_gene & tt_q1$has_prom, ]
  if (nrow(gene_cpgs) == 0) return(NULL)
  # For each CpG, each patient's β_T and β_N
  rows <- list()
  for (i in seq_len(nrow(gene_cpgs))) {
    cpg <- gene_cpgs$cpg[i]
    if (!cpg %in% rownames(beta_T)) next
    bT <- beta_T[cpg, ]; bN <- beta_N[cpg, ]
    rows[[length(rows)+1]] <- data.frame(
      cpg = cpg,
      feat = gene_cpgs$feat[i],
      patient = paired_pat,
      beta_T = bT, beta_N = bN,
      delta_beta = bT - bN,
      fdr = gene_cpgs$adj.P.Val[i],
      logFC = gene_cpgs$logFC[i],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

cat("Building VHL dataframe...\n")
vhl_df <- build_gene_df("VHL")
cat("VHL: ", length(unique(vhl_df$cpg)), "CpG rows, ", nrow(vhl_df), "patient-CpG pairs\n")

cat("Building SFRP1 dataframe...\n")
sfrp1_df <- build_gene_df("SFRP1")
cat("SFRP1:", length(unique(sfrp1_df$cpg)), "CpG rows, ", nrow(sfrp1_df), "patient-CpG pairs\n")

# ---- Order CpGs left-to-right by median Δβ for readability ----------------
order_cpgs <- function(df) {
  med_order <- df %>% group_by(cpg) %>% summarise(med = median(delta_beta)) %>%
    arrange(desc(med))
  df$cpg <- factor(df$cpg, levels = med_order$cpg)
  df
}
vhl_df   <- order_cpgs(vhl_df)
sfrp1_df <- order_cpgs(sfrp1_df)

# ---- Panel per gene: per-CpG strip of Δβ ---------------------------------
make_panel <- function(df, title, subtitle) {
  # Annotate summary stats per CpG
  ann <- df %>% group_by(cpg) %>%
    summarise(mean_db = mean(delta_beta),
              pct_pos = mean(delta_beta > 0),
              fdr = min(fdr), .groups = "drop")
  ann$label <- sprintf("mean Δβ = %.2f\n%.0f%% pos", ann$mean_db, 100*ann$pct_pos)

  ggplot(df, aes(x = cpg, y = delta_beta)) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_hline(yintercept = c(-0.2, 0.2), color = "grey80", linetype = "dashed") +
    geom_jitter(aes(color = delta_beta > 0), width = 0.15, size = 0.7, alpha = 0.55) +
    stat_summary(fun = median, geom = "crossbar", width = 0.6,
                 color = "black", size = 0.3) +
    geom_text(data = ann, aes(x = cpg, y = 1.25, label = label),
              size = 2.7, inherit.aes = FALSE, vjust = 1) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                       labels = c("FALSE" = "tumor lower", "TRUE" = "tumor higher"),
                       name = "Direction") +
    scale_y_continuous(limits = c(-1.0, 1.35), breaks = seq(-1, 1, 0.5)) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = expression(Delta * beta ~ "(tumor − normal)")) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

p_vhl <- make_panel(vhl_df,
  "VHL promoter CpGs (n = 7) — subgroup-driven, mixed direction",
  "Each dot = one paired patient's Δβ. Crossbar = median. Dashed lines at |Δβ| = 0.2.")

p_sfrp1 <- make_panel(sfrp1_df,
  "SFRP1 promoter CpGs (top 6 by effect size) — tissue-wide hypermethylation",
  "Each dot = one paired patient's Δβ. Crossbar = median.")

# Limit SFRP1 panel to top 6 by effect to keep figure readable
sfrp1_top <- sfrp1_df %>% group_by(cpg) %>%
  summarise(abs_median = abs(median(delta_beta)), .groups = "drop") %>%
  arrange(desc(abs_median)) %>% head(6) %>% pull(cpg)
sfrp1_df_top <- sfrp1_df[sfrp1_df$cpg %in% sfrp1_top, ]
sfrp1_df_top$cpg <- factor(sfrp1_df_top$cpg, levels = sfrp1_top)
p_sfrp1 <- make_panel(sfrp1_df_top,
  "SFRP1 promoter CpGs (top 6 by |median Δβ|) — tissue-wide hypermethylation",
  "Each dot = one paired patient's Δβ. Crossbar = median.")

# Combine
combined <- p_vhl / p_sfrp1 + plot_layout(heights = c(1, 1))
ggsave(file.path(OUT_DIR, "case_study_vhl_sfrp1.png"),
       combined, width = 8.5, height = 7.5, dpi = 150)

# Also save individual for flexibility
ggsave(file.path(OUT_DIR, "case_study_vhl.png"), p_vhl,
       width = 8, height = 4, dpi = 150)
ggsave(file.path(OUT_DIR, "case_study_sfrp1.png"), p_sfrp1,
       width = 8, height = 4, dpi = 150)

# Also save the underlying data for the final report
write.csv(vhl_df,   file.path(OUT_DIR, "vhl_perpatient_deltabeta.csv"),   row.names = FALSE)
write.csv(sfrp1_df, file.path(OUT_DIR, "sfrp1_perpatient_deltabeta.csv"), row.names = FALSE)

cat("\nDone. Figures in", OUT_DIR, "\n")
