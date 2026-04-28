# P8139 Final Project — Results Log

> Living document. Every analysis output lands here with enough detail to reproduce and interpret. `plan.md` has the plan + decisions; this file has the findings.
>
> Organized by workflow step, newest subsections appended. When a finding has an associated figure or CSV, reference the path under `project_workspace/`.

---

## Step 1 — EDA of TCGA-KIRC 450K methylation

*Run 2026-04-22. Script: `project_workspace/02_eda_kirc.R`. Outputs: `project_workspace/eda_outputs_kirc/`.*

### Cohort assembly

- Downloaded **485** open-access 450K methylation beta-value files (`project_workspace/data_kirc/GDCdata`, 5.9 GB).
- Assembled SummarizedExperiment: **485,577 CpG probes × 485 samples**.
- Sample-type breakdown:
  | Sample type | N |
  |---|---|
  | Primary Tumor (`01`) | 324 |
  | Solid Tissue Normal (`11`) | 160 |
  | Additional New Primary (`05`) | 1 |
- **Paired patients (both tumor and normal on 450K): 160** — matches GDC API pre-query exactly.
  - Paired IDs saved to `eda_outputs_kirc/paired_patient_ids.txt`.

### Missingness

- **Per-sample** (`eda_outputs_kirc/sample_missingness.csv`):
  - n=485, min 13.3%, max 24.8%, mean 14.8%.
  - Uniformly ≥13% across all samples → consistent with GDC's processed beta already having detection-p-value QC applied (failing calls set to NA).
- **Per-probe** across 485 samples:
  - Probes with zero missing: **330,174 / 485,577 (68%)** → sufficient complete-case set for differential analysis.
  - Probes with >5% samples missing: 95,091.
  - Probes with >20% samples missing: 76,722 — will be filtered in preprocessing.

### Clinical metadata quality

- 537 patient-level clinical records in `data_kirc/clinical.rds`; 319 have methylation.
- Per-field missingness (`eda_outputs_kirc/clinical_missingness.csv`):
  | Field | Missing |
  |---|---|
  | gender, race, ethnicity, vital_status, days_to_last_follow_up, primary_diagnosis, prior_malignancy, prior_treatment, ajcc_pathologic_t, ajcc_pathologic_n | 0% |
  | age_at_diagnosis | 0.2% |
  | ajcc_pathologic_m | 0.4% |
  | ajcc_pathologic_stage, tumor_grade | 0.6% |
  | days_to_death | 67% — expected (most patients alive) |
- Verdict: clinical metadata is clean enough to use directly as covariates.

### Unsupervised structure — PCA on top-10k-variance CpGs

- Computed on complete-case subset (330,174 probes), then reduced to top 10,000 by variance.
- Scree (first 5 PCs): **PC1 32.5%, PC2 19.6%, PC3 5.5%, PC4 3.8%, PC5 2.7%**.
- PC1–PC2 scatter: `eda_outputs_kirc/pca_top10k_PC12.png`
- Interpretation:
  - **PC1 cleanly separates tumor vs normal** → expected primary signal.
  - **PC2 is orthogonal to tumor/normal** — splits *both* sample types into an upper cluster (PC2 ≈ 10–15) and lower cluster (PC2 ≈ −5 to −7).
  - PC2 identity diagnosed in Step 1b below: **gender** (R² = 0.94). Addressed by dropping chr X / Y probes in preprocessing (Step 3a).

### Open technical notes from this step

1. **1 "Additional New Primary" sample** — decide in Step 3 whether to drop or merge with primary tumor.
2. **Probe filtering rules** for Step 3 preprocessing need to be fixed: complete-case requirement + cross-reactive probe list (Chen 2013) + SNP-overlap probes + sex-chromosome probes (if the question is sex-neutral).
3. **β → M transformation** for statistical testing to be applied before regression.

## Step 1b — PC1–PC10 diagnostic

*Run 2026-04-22. Script: `project_workspace/03_pc_diagnostic.R`.*

### Scree (top 10k variance CpGs, complete rows; n samples = 485)

| PC | Var |
|---|---|
| PC1 | 32.5% |
| PC2 | 19.6% |
| PC3 | 5.5% |
| PC4 | 3.8% |
| PC5 | 2.7% |
| PC6 | 1.4% |
| PC7 | 1.3% |
| PC8 | 0.9% |
| PC9 | 0.8% |
| PC10 | 0.7% |

- Elbow 1: PC2→PC3 (19.6 → 5.5). Elbow 2: PC5→PC6 (2.7 → 1.4).
- Major structure = **PC1–PC5 (~64% cumulative)**. PC6+ individually <1.5% and show no meaningful covariate associations → treated as noise.
- Figure: `eda_outputs_kirc/scree_top20.png`.

### PC ↔ covariate associations

For each PC score, associated every candidate covariate via one-way ANOVA (categorical) or Pearson correlation (numeric). R² + p-value in `eda_outputs_kirc/pc_covariate_association.csv`, heatmap in `eda_outputs_kirc/pc_covariate_heatmap.png`.

**Top driver per PC:**

| PC | Driver | R² | Type | Interpretation |
|---|---|---|---|---|
| PC1 | sample_type (tumor vs normal) | **0.78** | biological | Expected primary contrast |
| PC2 | **gender** | **0.94** | technical/biological | Sex-chromosome CpG dominance → must filter X/Y probes or adjust for gender |
| PC3 | tumor_grade (0.12), ajcc_T (0.11), vital_status (0.10), ajcc_stage (0.09) | ~0.11 | biological | Tumor aggressiveness / progression axis |
| PC4 | TSS (tissue source site) | 0.14 | technical (batch) | Must include as covariate |
| PC5 | TSS (0.11), plate (0.07) | ~0.11 | technical (batch) | Include TSS/plate as covariate |

### Implications for Step 2 / Step 3

1. **Sex handling**: drop chr X/Y probes (simplest) *and* include `gender` as a safety covariate. Otherwise PC2 noise dominates the tumor-vs-normal contrast.
2. **Batch adjustment**: include `TSS` (hospital of origin) as a covariate in every model. `plate` is likely collinear with TSS; add only if not too confounded.
3. **PC3 as secondary signal**: `tumor_grade` / `stage` / `vital_status` are mutually correlated and all load on PC3 — they point at a "tumor aggressiveness" axis. Candidate **secondary research question**: within tumor samples, which CpGs track tumor grade or T-stage?
4. The 1 "Additional New Primary" sample is near the primary-tumor cluster on PC1/PC2 — can be grouped with primary tumors or excluded without changing results.

## Step 3a — Preprocessing

*Script: `04_preprocess.R`. Outputs: `preprocess_out/`.*

Pipeline applied, in order:
| Step | Probes before | Probes after | Notes |
|---|---|---|---|
| Intersect SE probes with 450K annotation | 485,577 | 472,651 | 13K probes in SE not in official annotation — dropped |
| Drop chr X / chr Y | 472,651 | 461,295 | — |
| DMRcate `rmSNPandCH` (cross-reactive + SNP MAF>0.05) | 461,295 | 433,595 | Chen 2013 cross-reactive list applied via DMRcate |
| Drop probes with >20% sample missingness | 433,595 | 395,325 | — |
| Require complete-case across 484 samples | 395,325 | **320,953** | Final |

Additionally, the 1 `05A` "Additional New Primary" sample (`TCGA-DV-A4W0-05A-11D-A264-05`) was dropped per pre-registration → 484 samples.

β → M transform: M = log₂(β / (1−β)) with β clipped to [0.001, 0.999]. Final M matrix **320,953 × 484**, no missing cells, M range [−8.01, 8.19].

## Step 3b — Q1 primary: paired tumor vs normal differential methylation

*Script: `05_q1_paired.R`. Outputs: `q1_out/`.*

### Design amendment from plan

The pre-registered plan called for `limma + duplicateCorrelation` with covariates. In practice, this was infeasible for 320K probes × 320 samples (killed after 20 min). Switched to the canonical approach for balanced paired data: **within-pair differences D_i = M_tumor_i − M_normal_i, then moderated one-sample t-test on D via limma**. Mathematically, this design absorbs any patient-level covariate (age, gender, race, TSS) because those cancel in D_i; aliquot-level covariates (plate) are addressed in a sensitivity analysis. This change is a minor methodological amendment and is documented in the final Methods section.

### Primary results

- n = 160 paired patients, 320 samples analyzed.
- 320,953 CpGs tested.
- **FDR < 0.05: 222,386 probes (69.3%)**.
- FDR < 0.01: 190,232 probes.
- Bonferroni significant (α/n): 61,349 probes.
- Direction: among FDR<0.05 hits — 107,277 hypermethylated (tumor > normal), 115,109 hypomethylated → near-balanced.
- Genomic control λ_GC = **36.77** (interpretation below).
- Figures: `qq_q1.png`, `volcano_q1.png`, `manhattan_q1.png`.

### λ_GC interpretation (key methodological point)

A λ_GC of 36.77 would normally raise a red flag for GWAS. For **tumor-vs-normal methylation** it is expected:
1. **Permutation check** (`10_permutation_check.R`): randomly flipping each pair's tumor-normal sign and re-running on a 30K subset gives mean λ = **1.009** (sd = 0.39, range [0.61, 2.22]) across 20 permutations. So the pipeline is perfectly calibrated under the null; the observed λ = 37 reflects real pervasive signal, not confounding.
2. **Q-Q plot shape** (`qq_q1.png`): observed p-values track y = x at the origin and smoothly diverge upward — the signature of pervasive real signal, not uniform inflation.
3. **Biological expectation**: cancer methylomes are globally disordered, with large fractions of CpGs differentially methylated between tumor and normal. A λ close to 1 would actually be suspicious here (it would mean we missed the main biology).

Upshot: λ_GC should NOT be treated as a QC metric for contrasts where pervasive real change is expected. This is a point to explain carefully in the final report.

### Biologically significant effect sizes (beta-scale)

P-value is not effect size. Converting logFC (M-value) to β-value effect sizes for interpretability (`09_effectsize.R`):

| Criterion | N CpGs | % of tested |
|---|---|---|
| FDR < 0.05 | 222,386 | 69.3% |
| FDR < 0.05 **and** \|Δβ\| > 0.1 | **46,435** | 14.5% |
| FDR < 0.05 **and** \|Δβ\| > 0.2 | 9,312 | 2.9% |
| Bonferroni **and** \|Δβ\| > 0.1 | 46,156 | 14.4% |

Among "strong" hits (FDR<0.05 & |Δβ|>0.1): 22,942 hyper / 23,493 hypo — both directions well represented. This is the number to lead with in the report as the "reportable finding".

### Top 20 biologically strong hits (Q1)

From `effectsize_out/q1_top30_strong.csv`. First 10 by |t|:

| CpG | Chr | Gene | Feature | Δβ | p | FDR |
|---|---|---|---|---|---|---|
| cg13294602 | chr8 | EIF2C2 (AGO2) | Body | +0.55 (hyper) | 5×10⁻¹²¹ | 2×10⁻¹¹⁵ |
| cg12691620 | chr1 | C1orf163 | 3'UTR | −0.33 | 9×10⁻¹¹⁹ | 1×10⁻¹¹³ |
| cg13324103 | chr10 | SVIL | 5'UTR | −0.39 | 2×10⁻¹¹⁷ | 2×10⁻¹¹² |
| cg14601621 | chr9 | C9orf3 | 3'UTR | +0.52 | 2×10⁻¹¹⁷ | 2×10⁻¹¹² |
| cg17774001 | chr2 | DTNB | 3'UTR | +0.51 | 4×10⁻¹¹⁵ | 2×10⁻¹¹⁰ |
| cg08141142 | chr14 | MTA1 | 3'UTR | +0.50 | 9×10⁻¹¹² | 4×10⁻¹⁰⁷ |
| cg22274117 | chr6 | ATXN1 | 5'UTR | −0.46 | 1×10⁻¹¹¹ | 5×10⁻¹⁰⁷ |
| cg25247520 | chr8 | PVT1 | TSS200/Body | −0.48 | 2×10⁻¹⁰⁷ | 5×10⁻¹⁰³ |
| cg11201447 | chr8 | PVT1 | TSS200/Body | −0.57 | 2×10⁻¹⁰⁷ | 5×10⁻¹⁰³ |
| cg17469978 | chr7 | CAV1 | TSS200 | −0.27 | 2×10⁻¹⁰² | 2×10⁻⁹⁸ |

Biological sanity: **CAV1**, **PVT1**, **HDAC4**, **MTA1**, **ATXN1** are all well-known cancer genes; several specifically implicated in kidney cancer biology.

### Sensitivity analyses

- **S1 (unpaired, 484 samples with covariate adjustment for gender/age/TSS)**: FDR<0.05 → identified hits that overlap **89.8%** with paired primary. λ_GC = 29.2 (similarly high, as expected).
- **S2 (same-plate pairs only)**: only 5 pairs had mismatched plates (cross-plate); the 155-pair same-plate subset gives FDR<0.05 → **100% overlap** with primary. Plate confounding is negligible.

## Step 3c — Q3 secondary: within-tumor grade association

*Script: `06_q3_tumor_grade.R`. Outputs: `q3_out/`.*

### Sample selection

- 324 primary tumor samples before grade filter.
- Samples without assigned grade (GX, NA): dropped. Final **~310 tumors** across grades G1–G4.

### Primary results

- Model: limma, M ~ ordinal_grade + gender + age + TSS.
- FDR < 0.05: **99,724 probes (31.1%)**.
- λ_GC = 6.58; permutation and null-region calibration consistent with pervasive real signal at smaller effect sizes than Q1.
- After biological filter (FDR<0.05 AND |Δβ| across G1→G4 range > 0.1): **19,280 probes**.
- Top hits include **HOXB3**, **FOXD3**, **NF1**, **APCDD1**, **KCNQ1DN** — all well-known tumor-associated genes.
- Figures: `qq_q3.png`, `volcano_q3.png`.

### Sensitivity

- Grade as categorical factor (F-test across 4 levels) vs ordinal: very high overlap with ordinal primary.

## Step 3d — Q5 interpretation: GSEA / ORA

*Script: `07_q5_enrichment.R`. Outputs: `q5_out/`.*

### Q1 GSEA (MSigDB Hallmark + KEGG)

- **43 pathways** at padj<0.25.
- Top 10 mostly IMMUNE / INFLAMMATION related, negative NES → hypomethylated in tumor (consistent with immune infiltration in tumor tissue and/or tumor-intrinsic inflammation):
  - KEGG_OLFACTORY_TRANSDUCTION, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_COMPLEMENT, HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_INTERFERON_ALPHA_RESPONSE, HALLMARK_COAGULATION, HALLMARK_APOPTOSIS, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.
- **HALLMARK_HYPOXIA**: padj = 0.055, NES = −1.47 — hypoxia-response genes are hypomethylated in tumor → consistent with HIF activation.
- KEGG_RENAL_CELL_CARCINOMA: padj = 0.804 — not enriched in the tumor-vs-normal contrast.

### Q3 GSEA

- Top hits mixed direction:
  - **KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION**: NES +2.67 (genes hypermethylate with grade).
  - **HALLMARK_MYC_TARGETS_V1**, **HALLMARK_E2F_TARGETS**, **HALLMARK_G2M_CHECKPOINT**: NES −2.8, −2.7, −2.3 → hypomethylation with higher grade → consistent with proliferation activation.
  - **KEGG_RENAL_CELL_CARCINOMA**: padj = 0.030, NES = +1.82 → RCC pathway genes accumulate methylation changes with grade, supporting "higher grade → more canonical RCC epigenome".

## Step 3e — Q9 positive control (pre-registered)

*Script: `08_q9_positive_control.R`. Outputs: `q9_out/`. Overall verdict: **PASS**.*

Seven VHL promoter CpGs identified:

| CpG | Feature | Δβ (tumor − normal, via M) | FDR |
|---|---|---|---|
| cg03619761 | TSS1500 | negative | 2×10⁻¹⁰ |
| cg01289861 | TSS1500 | negative | 10⁻⁸ |
| cg13672843 | 1stExon | **POSITIVE** (large) | 3×10⁻⁷ |
| cg23977453 | TSS1500 | **positive** | 9×10⁻⁴ |
| cg13419702 | TSS1500 | ≈0 | 0.31 (n.s.) |
| cg15267345 | 5'UTR/1stExon | ≈0 | 0.47 (n.s.) |
| cg22730772 | TSS1500 | ≈0 | 0.74 (n.s.) |

**Pre-registered assertions**:
- **A1 (any VHL promoter CpG FDR<0.05 & logFC>0)**: PASS — cg13672843 (large hypermethylation in 1stExon) and cg23977453 are significantly hypermethylated in tumor.
- **A2 (HYPOXIA or RCC enriched at padj<0.25 in Q1 GSEA)**: PASS — HYPOXIA padj=0.055.
- **A3 (median logFC across FDR-sig VHL promoter CpGs > 0)**: PASS — median = +0.025.

**Important nuance for the report**: VHL promoter shows a **mixed methylation pattern** — some CpGs hypermethylated (likely alleles undergoing silencing), some hypomethylated (possibly regulatory-element demethylation in activated HIF context). Reporting this as a uniform "VHL silenced by methylation" would be inaccurate. The paper needs to acknowledge the bidirectional pattern explicitly.

## Step 4 — Self-audit (v2)

*Script: `99_audit.R`. Output: `audit_out/audit_report.csv` and `audit_report.md`.*

### Final verdict: **19 PASS, 2 WARN (explained), 0 FAIL**.

Audit v1 had flagged λ_GC > 3 as WARN for both Q1 and Q3. After adding:
- Null-region calibration test (obs/exp correlation on the least-significant 80% of p-values)
- Permutation-null validation (mean lambda = 1.009 on 20 permutations)
- Effect-size-aware hit counts (FDR AND |Δβ| > 0.1)

Audit v2 verdict breakdown:
- All preprocessing integrity checks: PASS.
- Hit counts (both any FDR and strong-effect) consistent with cancer biology: PASS.
- Direction balance (bidirectional, not one-sided): PASS.
- Sensitivity overlaps (unpaired 90%, same-plate 100%): PASS.
- Positive-control assertions (A1/A2/A3): PASS.
- GSEA pathway-specific expectations (HYPOXIA Q1, RCC Q3): PASS.
- Reproducibility (sessionInfo at each step): PASS.
- The two WARN are on the global lambdas (Q1=36.77, Q3=6.58) being higher than conventional GWAS-style thresholds; both are documented as expected-and-not-a-concern based on permutation check (true-null lambda = 1.01) and Q-Q plot shape.

---

