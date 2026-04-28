#!/usr/bin/env Rscript
# 25_methylmix_benchmark.R
# Run MethylMix on the same TCGA-KIRC data and compare its
# MethylationDrivenGenes list to our 4,483-gene silencing list.
#
# Strategy:
#   - Aggregate β to gene-level by averaging over promoter CpGs (TSS1500/TSS200/
#     5'UTR/1stExon), matching the (CpG, gene) annotation our pipeline uses.
#   - Build METcancer (gene × tumor), METnormal (gene × normal), GEcancer
#     (gene × tumor), with tumor samples in matching column order between
#     METcancer and GEcancer.
#   - MethylMix() uses BiocParallel; we use MulticoreParam.
#   - Save the gene list, runtime, and overlap stats.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(MethylMix)
  library(BiocParallel)
  library(data.table)
})

set.seed(8139)
WORKDIR <- "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
setwd(WORKDIR)
PRE_DIR <- file.path(WORKDIR, "preprocess_out")
DATA_DIR <- file.path(WORKDIR, "data_kirc")
SIL_DIR <- file.path(WORKDIR, "silencing_out")
OUT_DIR <- file.path(WORKDIR, "methylmix_out")
dir.create(OUT_DIR, showWarnings = FALSE)

t_start <- Sys.time()

# ---- 1. load M-matrix and convert to β -------------------------------------
cat("[1/7] loading M matrix...\n")
M    <- readRDS(file.path(PRE_DIR, "M_matrix.rds"))
meta <- readRDS(file.path(PRE_DIR, "sample_metadata.rds"))
M2beta <- function(m) 2^m / (2^m + 1)
beta_full <- M2beta(M); rm(M); gc(verbose = FALSE)
cat(sprintf("  beta matrix: %d CpGs × %d aliquots\n",
            nrow(beta_full), ncol(beta_full)))

# Aggregate replicate aliquots within (patient, sample_type)
keys <- paste(meta$patient, meta$sample_type, sep = "|")
unique_keys <- unique(keys)
beta_agg <- vapply(unique_keys, function(k) {
  ix <- which(keys == k)
  if (length(ix) == 1) beta_full[, ix] else
    rowMeans(beta_full[, ix, drop = FALSE], na.rm = TRUE)
}, numeric(nrow(beta_full)))
colnames(beta_agg) <- unique_keys
rm(beta_full); gc(verbose = FALSE)

# Build sample-type lookup
sample_type <- vapply(unique_keys, function(k) strsplit(k, "\\|")[[1]][2],
                     character(1))
tumor_keys  <- unique_keys[sample_type == "Primary Tumor"]
normal_keys <- unique_keys[sample_type == "Solid Tissue Normal"]
cat(sprintf("  unique keys: %d tumor + %d normal\n",
            length(tumor_keys), length(normal_keys)))

# ---- 2. aggregate CpG → gene (promoter mean β) -----------------------------
cat("[2/7] aggregating CpGs to gene level (promoter mean β)...\n")
# Use the (CpG, gene) mapping already built in silencing_out/silencing_table.csv
sil_tab <- fread(file.path(SIL_DIR, "silencing_table.csv"))
cpg2gene <- unique(sil_tab[, .(cpg, gene)])
cat(sprintf("  promoter (CpG, gene) mappings: %d unique pairs, %d unique genes\n",
            nrow(cpg2gene), uniqueN(cpg2gene$gene)))

# Restrict to CpGs present in beta matrix
cpg2gene <- cpg2gene[cpg %in% rownames(beta_agg)]
genes <- unique(cpg2gene$gene)
cat(sprintf("  after intersection with beta matrix: %d genes\n", length(genes)))

# Build gene × sample matrix by averaging β across each gene's promoter CpGs.
# Use rowsum() for vectorized aggregation — orders of magnitude faster than
# per-gene subsetting.
keep_cpgs <- cpg2gene$cpg
gene_factor <- cpg2gene$gene
beta_sub <- beta_agg[keep_cpgs, , drop = FALSE]   # CpGs (with potential repeats if a CpG maps to multiple genes)
gene_sum   <- rowsum(beta_sub, group = gene_factor, reorder = FALSE, na.rm = TRUE)
gene_count <- rowsum(matrix(as.numeric(!is.na(beta_sub)),
                             nrow = nrow(beta_sub), ncol = ncol(beta_sub),
                             dimnames = list(rownames(beta_sub), colnames(beta_sub))),
                     group = gene_factor, reorder = FALSE)
gene_beta <- gene_sum / pmax(gene_count, 1)
gene_beta[gene_count == 0] <- NA_real_
rm(beta_agg, beta_sub, gene_sum, gene_count); gc(verbose = FALSE)
cat(sprintf("  gene-level β matrix: %d genes × %d samples\n",
            nrow(gene_beta), ncol(gene_beta)))

# ---- 3. load gene expression and align -------------------------------------
cat("[3/7] loading RNA-seq, aligning samples...\n")
rna_se <- readRDS(file.path(DATA_DIR, "kirc_rnaseq_SE.rds"))
rna_cd <- as.data.frame(colData(rna_se))[, c("barcode","patient","sample_type")]
keep <- rna_cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
rna_se <- rna_se[, keep]; rna_cd <- rna_cd[keep, , drop = FALSE]
gene_names <- as.character(rowData(rna_se)$gene_name)
tpm <- assay(rna_se, "tpm_unstrand")
rna_keys <- paste(rna_cd$patient, rna_cd$sample_type, sep = "|")
unique_rna_keys <- unique(rna_keys)
tpm_avg <- vapply(unique_rna_keys, function(k) {
  ix <- which(rna_keys == k)
  if (length(ix) == 1) tpm[, ix] else rowMeans(tpm[, ix, drop = FALSE])
}, numeric(nrow(tpm)))
colnames(tpm_avg) <- unique_rna_keys
keep_g <- !is.na(gene_names) & gene_names != ""
tpm_avg <- tpm_avg[keep_g, , drop = FALSE]
gene_lab <- gene_names[keep_g]
log2tpm <- log2(rowsum(tpm_avg, group = gene_lab, reorder = FALSE) + 1)

# Subset both matrices to common genes
common_genes <- intersect(rownames(gene_beta), rownames(log2tpm))
cat(sprintf("  common genes between methylation and RNA-seq: %d\n",
            length(common_genes)))
gene_beta <- gene_beta[common_genes, , drop = FALSE]
log2tpm   <- log2tpm[common_genes, , drop = FALSE]

# Build the three MethylMix matrices.
# Tumor samples: keys present in BOTH gene_beta tumor cols AND log2tpm.
tumor_keys_both <- intersect(tumor_keys, colnames(log2tpm))
normal_keys_meth <- normal_keys  # all methylation normals
cat(sprintf("  paired tumor (methylation × RNA-seq): %d, methylation normals: %d\n",
            length(tumor_keys_both), length(normal_keys_meth)))

METcancer <- gene_beta[, tumor_keys_both, drop = FALSE]
METnormal <- gene_beta[, normal_keys_meth, drop = FALSE]
GEcancer  <- log2tpm[, tumor_keys_both, drop = FALSE]

# Sanity
stopifnot(identical(colnames(METcancer), colnames(GEcancer)))
stopifnot(identical(rownames(METcancer), rownames(GEcancer)))
stopifnot(identical(rownames(METcancer), rownames(METnormal)))
cat(sprintf("  METcancer: %d × %d, METnormal: %d × %d, GEcancer: %d × %d\n",
            nrow(METcancer), ncol(METcancer),
            nrow(METnormal), ncol(METnormal),
            nrow(GEcancer),  ncol(GEcancer)))

# Drop rows with any NA in METcancer/METnormal (MethylMix requires complete rows)
ok <- complete.cases(METcancer) & complete.cases(METnormal) & complete.cases(GEcancer)
cat(sprintf("  dropping %d genes with NA; %d retained\n",
            sum(!ok), sum(ok)))
METcancer <- METcancer[ok, , drop = FALSE]
METnormal <- METnormal[ok, , drop = FALSE]
GEcancer  <- GEcancer[ok, , drop = FALSE]

# Save the input matrices for reference / re-runs
saveRDS(list(METcancer = METcancer, METnormal = METnormal, GEcancer = GEcancer),
        file.path(OUT_DIR, "methylmix_inputs.rds"))

# ---- 4. run MethylMix ------------------------------------------------------
cat("[4/7] running MethylMix...\n")
n_cores <- min(8L, parallel::detectCores())
cat(sprintf("  using %d cores\n", n_cores))
register(MulticoreParam(workers = n_cores))
t0 <- Sys.time()
mm <- MethylMix(METcancer = METcancer,
                GEcancer  = GEcancer,
                METnormal = METnormal,
                listOfGenes = NULL,
                filter = TRUE,
                NoNormalMode = FALSE,
                OutputRoot = "")
t1 <- Sys.time()
runtime <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("  MethylMix runtime: %.1f sec (%.2f min)\n",
            runtime, runtime / 60))

# ---- 5. save MethylMix output ----------------------------------------------
cat("[5/7] saving MethylMix results...\n")
saveRDS(mm, file.path(OUT_DIR, "methylmix_result.rds"))
mm_genes <- mm$MethylationDrivenGenes
cat(sprintf("  MethylationDrivenGenes: %d\n", length(mm_genes)))
writeLines(mm_genes, file.path(OUT_DIR, "methylmix_genes.txt"))

# ---- 6. overlap with our silencing list ------------------------------------
cat("[6/7] computing overlap with our silencing genes...\n")
our_silencing <- unique(sil_tab$gene[sil_tab$silencing == TRUE])
N_back <- length(common_genes)
N_ours <- length(our_silencing)
N_mm   <- length(mm_genes)
N_both <- length(intersect(our_silencing, mm_genes))
N_only_ours <- length(setdiff(our_silencing, mm_genes))
N_only_mm   <- length(setdiff(mm_genes, our_silencing))

cat(sprintf("  background: %d common-tested genes\n", N_back))
cat(sprintf("  our silencing: %d  | MethylMix: %d  | both: %d\n",
            N_ours, N_mm, N_both))
cat(sprintf("  only ours: %d  | only MethylMix: %d\n",
            N_only_ours, N_only_mm))

jaccard <- N_both / length(union(our_silencing, mm_genes))
cat(sprintf("  Jaccard: %.3f\n", jaccard))

# Hypergeometric on overlap
hyp_p <- phyper(N_both - 1, N_mm, N_back - N_mm, N_ours, lower.tail = FALSE)
cat(sprintf("  hypergeometric p (overlap > random): %.2e\n", hyp_p))

# Recovery of pre-registered TSGs by MethylMix
TSG18 <- c("VHL","RASSF1","SFRP1","SFRP2","SFRP4","SFRP5","DKK3","GATA5",
           "BNC1","COL14A1","PCDH17","SLIT2","CDKN2A","APAF1",
           "TIMP3","NEFH","UCHL1","PITX2")
mm_tsg <- intersect(TSG18, mm_genes)
ours_tsg <- intersect(TSG18, our_silencing)
cat(sprintf("\n  18-TSG recovery — ours: %d/18, MethylMix: %d/18, both: %d/18\n",
            length(ours_tsg), length(mm_tsg), length(intersect(mm_tsg, ours_tsg))))
cat(sprintf("    ours: %s\n", paste(ours_tsg, collapse = ", ")))
cat(sprintf("    mm:   %s\n", paste(mm_tsg,   collapse = ", ")))

# ---- 7. write summary ------------------------------------------------------
cat("[7/7] writing summary...\n")
summary_dt <- data.frame(
  metric = c("background_genes",
             "ours_n", "mm_n", "overlap_n",
             "ours_only_n", "mm_only_n",
             "jaccard", "hypergeometric_p",
             "ours_TSG_n", "mm_TSG_n", "both_TSG_n",
             "mm_runtime_sec"),
  value  = c(N_back,
             N_ours, N_mm, N_both,
             N_only_ours, N_only_mm,
             round(jaccard, 4), formatC(hyp_p, format = "e", digits = 2),
             length(ours_tsg), length(mm_tsg),
             length(intersect(mm_tsg, ours_tsg)),
             round(runtime, 1)))
fwrite(summary_dt, file.path(OUT_DIR, "benchmark_summary.csv"))

t_end <- Sys.time()
cat(sprintf("\nDone. Total: %.1f min\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))
cat(sprintf("Outputs in %s\n", OUT_DIR))
