# P8139 Final Project — Working Plan

> Living document. Updated every time we make a substantive decision.
> If the conversation session dies, a new session should be able to read this file + `lecture_methods.md` + `results.md` and pick up where we left off.
>
> **File convention**:
> - `plan.md` (this file) — plan, decisions, open questions, shortlists, write-up checklist. Short, navigable.
> - `results.md` — every analysis finding with numbers, tables, figure references. Append-only, grows with the project.
> - `lecture_methods.md` — in-scope method catalog.
> - `project_workspace/` — code, downloaded data, raw outputs.

---

## 0. Project ground rules

- **Course**: P8139 Statistical Genetic Modeling (Columbia Biostatistics, Spring 2026).
- **Deliverables**:
  - HW5 (due 2026-04-28): short proposal / progress report on the final-project data.
  - Final report (due 2026-05-15): ≤10 pages, single-spaced, 11pt, journal-article style.
- **Working principle**: project first, proposal last. We do the analysis, get real results, then write the HW5 proposal and the final report describing what we actually did. HW5 can be submitted as a scoped-down snapshot of this document.
- **Methods discipline**: lean on methods in `lecture_methods.md`. Going beyond is allowed when necessary but must be declared. No "let's throw a fancy method at it because it sounds good".
- **Honesty discipline**: null results are fine. Broken pipelines presented as findings are not. Every reported result must be reproducible from the code we commit.

## 1. TCGA in one page (for ourselves)

*See Section 1 of conversation 2026-04-22 for full briefing. Key facts for decision-making:*

- **Scope**: ~11,000 patients across ~33 cancer projects ("TCGA-XXXX"), closed enrollment (2006–2018). Hosted on NCI GDC (Genomic Data Commons).
- **Study design**: not uniform — some projects are rich in adjacent-normal samples (paired), most are tumor-only. This dictates whether paired t-tests are even an option.
- **Data categories commonly available**:
  - **DNA methylation**: Illumina 450K (most projects) and 27K arrays (older, smaller). β-values in [0,1]; M-values recommended for testing.
  - **Gene expression**: RNA-Seq (STAR counts, TPM, FPKM) — broadest coverage.
  - **miRNA expression**: miRNA-Seq.
  - **Copy number variation**: SNP6 array, masked/unmasked CNV.
  - **Somatic mutation (SNV/indel)**: MAF files from Mutect2, VarScan, etc.
  - **Clinical**: demographics, stage, survival, treatment.
- **Sample barcodes** encode tissue type in positions 14–15 (`01` primary tumor, `11` solid tissue normal, `06` metastatic, etc.). Pair detection = same patient (first 12 chars of barcode) with both `01` and `11` samples.
- **Access tiers**:
  - *Open access*: processed/normalized data (methylation β-values, gene-level counts) — what we will use.
  - *Controlled access*: raw sequencing, germline genotypes — requires dbGaP approval; we will not use.
- **Typical adjacent-normal-rich projects** (≥30 paired 450K methylation samples historically): BRCA, LUAD, LUSC, COAD, KIRC, KIRP, LIHC, THCA, HNSC, PRAD. We will verify these numbers in Step 1, not trust them blindly.

## 2. Methods scope

Canonical list in `lecture_methods.md`. Short version of what we expect to use:

- Per-feature differential test (paired t-test or linear regression) on methylation M-values / expression log-counts.
- Multiple-comparison correction: Bonferroni (strict), Benjamini-Hochberg FDR (primary), Storey q-value (optional).
- PCA for unsupervised structure + batch/ancestry adjustment.
- Permutation for either empirical p-values per locus or a global test statistic.
- Manhattan + Q-Q + volcano plots.
- GSEA or ORA on the ranked / significant feature list.
- (If applicable) HWE filter on any genotype data we end up including.

## 3. Workflow

Iterative, with explicit exit criteria at each step so we don't wander.

### Step 1 — Dataset triage + EDA
- Query TCGAbiolinks for N candidate (cancer × omics) combos.
- TCGA coverage is sparse — **do not assume uniform availability**. For each candidate, record:
  1. Does this project have this data category at all? (some small projects skip categories)
  2. n tumor samples with this data type
  3. n independent-normal samples with this data type
  4. n **paired** patients where *both* tumor and normal have this data type (not "has normal" — has normal *with this assay*)
  5. If multi-omics: n patients in the intersection of all required data types (often much smaller than any single type)
  6. Clinical covariate coverage for the fields the planned question needs
  7. Methylation: how many on 450K vs 27K array (do not mix)
  8. Expression: which workflow (STAR counts for DE, TPM for correlation/PCA)
- **Exit criterion**: at least one combo with (a) ≥30 paired samples *with matched assay on both sides* **or** ≥100 tumor + ≥20 independent normal, (b) metadata fields needed for the planned question with <20% missingness, (c) manageable download size.

### Step 2 — Research question + analysis plan
- 1–3 concrete research questions per surviving dataset.
- Each maps to a primary test and a covariate set.
- Power / feasibility sanity check (not a formal power calc; just "is n large enough for this effect scale?").
- **Exit criterion**: question is a sentence a reviewer would understand; plan cites a specific method from `lecture_methods.md`; multiple-comparison strategy specified up front.

### Step 3 — Execute
- Reproducible R script; seed set; `sessionInfo()` captured alongside results.
- Save intermediate artifacts: sample list used, feature-level p-value table, figures — all to disk, committable.
- Run the primary test, then diagnostics (Q-Q, volcano, PCA of residuals), then secondary/sensitivity analyses.
- **Exit criterion**: code runs end-to-end from raw TCGA download to result tables/plots without manual intervention.

### Step 4 — Audit
- Assumption checks, effect-size sanity (are the top CpGs/genes biologically plausible?), covariate balance, negative controls if available.
- Compare to at least one sanity baseline (e.g. permuted labels should give flat Q-Q).
- Decide: **accept** / **revise plan** (back to Step 2) / **drop dataset**.
- **Exit criterion**: written audit note in the decision log below.

Loop Step 1→4 over multiple (dataset, question) pairs. At the end, pick the strongest one with the user and write it up.

**Pre-registration / anti-p-hacking rules** (agreed 2026-04-22):
- Every (dataset, question) we try gets a row in Section 5 of this file *before* running the analysis.
- Every row is disclosed in the final report — even "tested, null" — so we don't selectively report the one dataset that happened to give a small p.
- Final pick criterion is **effect plausibility × story coherence × methodological cleanliness**, not "lowest p".

**Step 0 — environment setup** (to do before Step 1):
- Confirm `TCGAbiolinks`, `minfi` / `sesame`, `limma`, `DESeq2` install cleanly.
- Pin Bioconductor version in `sessionInfo()`.
- Confirm disk space for at least two cohort-level methylation downloads.

## 4. Candidate datasets — living shortlist

Populated from GDC API query 2026-04-22 (`project_workspace/query_gdc_methylation.py` → `gdc_methylation_counts.json`).

Filter: DNA methylation, Illumina Human Methylation 450 platform, open access, TCGA program. "Paired 450K" = unique patients with *both* primary tumor and solid tissue normal 450K files.

**Top candidates (≥29 paired):**

| ID | Cancer | N tumor 450K | N normal 450K | **N paired 450K** | Total size (GB) | Notes | Status |
|---|---|---|---|---|---|---|---|
| D01 | TCGA-KIRC (ccRCC) | 319 | 160 | **160** | 6.4 | Largest paired n by far; VHL/HIF biology known; single histology | candidate |
| D02 | TCGA-BRCA | 784 | 97 | **91** | 11.7 | Rich molecular subtypes (ER/PR/HER2, PAM50); largest cohort | candidate |
| D03 | TCGA-THCA | 507 | 56 | **56** | 7.5 | | candidate |
| D04 | TCGA-HNSC | 528 | 50 | **50** | 7.6 | HPV+/− subgroup possible | candidate |
| D05 | TCGA-PRAD | 498 | 50 | **50** | 7.3 | | candidate |
| D06 | TCGA-LIHC | 377 | 50 | **50** | 5.6 | HBV/HCV context | candidate |
| D07 | TCGA-KIRP (papillary) | 275 | 45 | **45** | 4.2 | Kidney sibling cohort | candidate |
| D08 | TCGA-LUSC | 370 | 42 | **40** | 5.4 | | candidate |
| D09 | TCGA-COAD | 295 | 38 | **38** | 4.6 | Classic CIMP phenotype | candidate |
| D10 | TCGA-UCEC | 431 | 46 | **33** | 6.4 | | candidate |
| D11 | TCGA-LUAD | 458 | 32 | **29** | 6.7 | Borderline; smoking covariate available | borderline |

**Excluded (paired n too small for tumor-vs-normal, or platform mismatch):**

BLCA (21), ESCA (16), PAAD (10), CHOL (9), READ (7), SARC (4), CESC (3), PCPG (3), STAD (2), THYM (2), GBM (1), LGG/LAML/TGCT/MESO/UVM/ACC/KICH/UCS/DLBC (all 0), SKCM (0 paired — mostly metastatic), OV (only 10 on 450K; its 592 cases are on older 27K array — **not pooled with 450K**).

Total 450K across TCGA: 9812 files, **128.6 GB** (full download not needed).

## 5. Research questions — pre-registered 2026-04-22

Locked-in three-layer narrative. Every question below is pre-registered *before* running the analysis. All three must appear in the final report (even if null). The analysis plan is specified in full here so Step 3 just executes.

| ID | Type | Question | Design | Primary test | Covariates | MCC |
|---|---|---|---|---|---|---|
| **Q1** | Primary | Which CpGs are differentially methylated between primary KIRC tumor and matched adjacent-normal tissue? | Paired, n = 160 patients | Paired test on M-values (paired t-test via `limma` with patient as blocking factor, or one-sample t-test on within-pair M-value differences) | gender, age, TSS (nuisance) | BH FDR q < 0.05; report Bonferroni-significant count too |
| **Q3** | Secondary | Within primary tumors (n = 324), which CpGs' methylation levels are associated with tumor grade (G1/G2/G3/G4)? | Unpaired, within-tumor | Linear regression: M ~ grade (ordinal) + gender + age + TSS | same + patient not blocked | BH FDR q < 0.05 |
| **Q5** | Interpretation | What biological pathways are enriched among Q1 hits (and separately Q3 hits)? | Ranked-list enrichment | GSEA on signed -log10(p) ranked CpG→gene mapping; ORA with hypergeometric test as a robust check | — | GSEA default perm-p; ORA BH FDR |
| **Q9** | Positive control | Are the VHL/HIF-axis methylation signals recovered? | As in Q1/Q5 | VHL promoter CpGs hypermethylation test (pre-specified); HIF/hypoxia pathway enrichment check | as Q1 | Pre-specified, *not* part of multiple-testing family |

**Excluded (documented in report)**: Q2 DMR calling (possible minor extension to Q1), Q4 stage (too correlated with Q3), Q6 sex×tissue interaction (complexity / 10-page limit), Q7 unsupervised subtyping (not in course), Q8 survival (Cox not in course).

### Pre-registered analysis plan (writeup-ready)

**Sample selection**:
- Q1 input = 160 patients with *both* primary tumor (`-01A`) and adjacent-normal (`-11A`) on 450K. Drop the 1 "Additional New Primary" sample.
- Q3 input = 324 primary-tumor samples (including all patients, not just paired).

**Probe preprocessing (common to Q1/Q3)**:
1. Start from TCGAbiolinks β-value matrix (485,577 CpGs).
2. Drop probes on chr X / chr Y (sex dominates PC2 per EDA).
3. Drop cross-reactive probes (Chen et al. 2013 list).
4. Drop SNP-overlap probes (annotation-based).
5. Drop probes with >20% sample missingness; then keep complete-case probes for testing.
6. Transform β → M = log2(β / (1−β)), with β clipped to [0.001, 0.999] to avoid ±∞.

**Statistical test (Q1)**:
- Primary: `limma`-based moderated paired t-test on M-values, with `lmFit` design matrix including `tumorVsNormal` as main coefficient and `gender`, `age`, `TSS` as adjustment. Patient ID as `block` + `duplicateCorrelation` to borrow strength across pairs.
- Alternative (sanity check): naive paired t-test on within-pair M-differences, no covariates.

**Statistical test (Q3)**:
- Primary: `limma` linear model M ~ ordinal(grade) + gender + age + TSS, coefficient on `grade`.
- Alternative: Kendall's tau per CpG (distribution-free sanity check).

**Multiple testing**:
- Primary: BH FDR with threshold q < 0.05.
- Reported alongside: Bonferroni hits; Q-Q plot vs uniform null; permutation-based genomic control λ.

**Sensitivity analyses (pre-registered)**:
- Repeat Q1 without covariate adjustment → quantify how much sex/age/TSS inflates null.
- Re-do Q1 on the 164 unpaired tumor + 160 normal samples with an independent-samples design → compare overlap of top hits.
- Exclude any CpG with even 1% missingness → confirm findings are not driven by imputation / missingness.

**Pathway / biology (Q5)**:
- GSEA on signed `-log10(p) × sign(logFC)` ranked list of CpG→gene mappings (use Illumina 450K → gene map; if multiple CpGs per gene, take max absolute statistic).
- ORA (hypergeometric) on FDR-significant gene list vs 450K-array gene universe.
- Gene sets: MSigDB Hallmark (50 sets) + KEGG renal cell carcinoma + KEGG TCA cycle + selected Hypoxia/HIF-response gene sets.

**Positive-control pre-registration (Q9)**:
Three pre-specified assertions, *not* part of the family-wide MCC, written *before* results:
1. At least one VHL promoter CpG (probes mapping to `VHL` TSS ± 1500 bp) will appear hypermethylated in tumor vs normal with FDR < 0.05.
2. `HALLMARK_HYPOXIA` and/or `KEGG_RENAL_CELL_CARCINOMA` will enrich in Q5 (GSEA FDR < 0.25).
3. VHL promoter methylation direction will be positive (tumor > normal).
If any assertion fails, pause Step 3 and audit the pipeline (Step 4 early).

## 6. Current status

- **Last updated**: 2026-04-23
- **Locked-in**: omics = **DNA methylation**; cohort = **TCGA-KIRC**.
- **Report framing**: PIVOTED AGAIN to **methodological case study (3-layer framework)** with KIRC as primary application. Tumor-vs-normal biomarker discovery turned out too trivially easy (V2 degenerate) so the contribution becomes the Layer 1/2/3 framework itself.
- **Step 3 (execute analyses)**: DONE. Pipeline 04→10 ran clean. Q9 positive control PASS (3/3 assertions). Permutation null lambda = 1.009 → method calibrated.
- **Step 4 (audit v2)**: DONE. 19 PASS, 2 WARN (explained), 0 FAIL.
- **Biomarker pipeline (Step 5)**: DONE with Option B gate (post code-audit fix: removed erroneous CDKN2B alias of CDKN2A in V1 literature matcher). **V1 lenient PASS: 12/20 literature biomarkers have ≥1 CpG in Layer 1 pool, 15/20 have ≥1 CpG in Layer 2 pool.** V2 degenerate (top-50 AUC = 0.99996 indistinguishable from random L2 baselines — "problem too easy"). V3 underpowered (26 genes in top-50).
- **Active task**: write final report under methodological framing.
- **Blocked on**: nothing.

### Step 1 EDA — headline only (full details in `results.md`)
- Cohort assembled: 485 samples × 485,577 CpG. 160 paired tumor/normal patients confirmed.
- Clinical metadata clean (<1% missing on key fields).
- PCA: **PC1 (32.5%) = tumor vs normal**, as expected. **PC2 (19.6%) is an unknown latent axis**, orthogonal to tumor/normal — splits both sample types. Must diagnose before Step 2.
- → See `results.md` § "Step 1 — EDA of TCGA-KIRC 450K methylation" for numbers, tables, figures, open technical notes.

## 5b. Biomarker discovery operational plan (2026-04-23 amendment)

After Step 3 results came in, a methodological discussion surfaced several concerns that reshape *how* we report findings:
- 222K FDR-significant CpGs is statistical power, not biological reality (cell composition, purity, tissue heterogeneity bias).
- T-test null (mean = 0) conflates "population-wide small shifts" with "subgroup-driven large changes".
- "Tissue-wide biomarker" is a clinical concept (gene/pathway level), not a statistical one (single CpG).
- Need stronger validation than a single positive control.

**Decision**: pivot primary framing from descriptive Q1/Q3/Q5/Q9 to **biomarker discovery**; keep Atlas (§5) as fallback. Original Q1/Q3/Q5/Q9 pre-registered analyses remain intact — they *feed* the new workflow.

### Operational definition of "tissue-wide KIRC methylation biomarker"

A CpG that satisfies all four:
1. Statistically reliable (FDR < 0.05).
2. Biologically meaningful effect (|mean Δβ| > 0.1).
3. Consistent direction across patients (pct_consistent > 0.8, i.e. ≥128/160 pairs same direction).
4. Strong single-marker classification performance (per-CpG unpaired AUC).

### Workflow

**Stage 1 — Gate (three filters on Q1 results, revised 2026-04-23 after empirical check of filter redundancy)**:
1. FDR < 0.05 (from Q1 paired moderated t)
2. **disc_auc > 0.85** (per-CpG unpaired AUC, discrimination threshold; scale-invariant replacement for |Δβ|>0.1)
3. pct_consistent > 0.8 (binomial null p ≈ 10⁻¹⁷; within-patient directional penetrance)

*Rationale for dropping |Δβ|>0.1 and adding AUC*: β is a bounded proportion whose biological meaning is scale-dependent (0→0.1 is categorically different from 0.5→0.6). AUC is scale-invariant and directly measures discrimination utility. Empirical check (`14_consistency_vs_auc.R`) showed r(|Δβ|>0.1 mask, AUC>0.85 mask) ≈ 1.0 in this dataset, so this is mostly a conceptual refinement rather than a large numerical change. Both disc_auc > 0.85 and pct_consistent > 0.8 kept because they measure conceptually distinct properties (clinical discrimination vs biological penetrance): correlation r=0.90 but 29K CpGs pass pct_consistent yet fail AUC, confirming they are complementary not redundant.

**Stage 2 — Rank & truncate**:
- Per-CpG unpaired AUC (tumor vs normal ROC, descending)
- Top 50 CpGs, with **one-CpG-per-gene constraint** (avoid panel dominated by neighboring CpGs of same gene)

**Stage 3 — Three-way validation**:

| ID | Type | Method | Pass criterion |
|---|---|---|---|
| V1 | Concurrent | 20-biomarker literature list (Tier A+B); test gene-level overlap with top-50 and broader Layer-1 pool; hypergeometric enrichment | ≥12/20 recovered in top-50 or Layer-1 pool |
| V2 | Predictive | Patient-level 5-fold CV (same patient's tumor+normal in same fold) of logistic regression on top-50 M values; compare vs L1-baseline (50 random from FDR<0.05, 30 reps) and L2-baseline (50 random from FDR<0.05 ∩ \|Δβ\|>0.1, 30 reps) | Top-50 AUC > L2-baseline 95th percentile |
| V3 | Biological | GSEA + ORA on top-50 mapped genes | ≥1 KIRC-related pathway (RCC/HIF/Wnt/cancer) padj < 0.05 |

**Stage 4 — Decision rules (pre-specified before running)**:

| Scenario | Criteria | Action |
|---|---|---|
| **Success** | V1 ≥ 12/20 AND V2 top-50 AUC > L2-baseline 95th pct AND V3 ≥ 1 pathway padj<0.05 | Write Biomarker paper |
| **Partial** | Any 1 validation fails | Biomarker paper with honest limitations |
| **Failure** | V1 < 8/20 AND (V2 AUC − L2-baseline mean) < 0.02 | **Fallback to Atlas** (§5 reused, narrative restructured to 3-layer classification) |

### Pre-registered literature biomarker list (V1 reference)

All 20 with expected layer (predicted from literature methylation frequency) — fixed before V1 evaluation.

| Gene | Tier | Expected layer | Reason |
|---|---|---|---|
| SFRP1 | A | L1 | Hypermethylated in >60% ccRCC (Morris 2010, Urakami 2006) |
| SFRP2 | A | L1 | Morris 2010 |
| RASSF1A | A | L1 | Frequently methylated in multiple cancers including RCC (Dreijerink 2001) |
| GATA5 | A | L1 | Peters 2014, Morris 2010 |
| DKK3 | A | L1 | Wnt antagonist, Urakami 2006 |
| BNC1 | B | L1 | Morris 2010 diagnostic panel |
| COL14A1 | B | L1 | Morris 2010 diagnostic panel |
| PCDH17 | B | L1 | Costa 2011, often >50% methylated |
| SLIT2 | B | L1 | Astuti 2011 |
| SFRP4 | B | L2 | Less consistently reported than SFRP1 |
| SFRP5 | B | L2 | Less consistently reported |
| VHL | A | L2 | Only ~15–20% of ccRCC silence VHL via methylation |
| CDKN2A | B | L2 | Dulaimi 2004; variable across cohorts |
| APAF1 | B | L2 | Christoph 2006; moderate frequency |
| TIMP3 | B | L2 | Bachman 1999; moderate |
| NEFH | B | L2 | Revelo 2005 |
| UCHL1 | B | L2 | Seliger 2010 |
| PITX2 | B | L2 | Dietrich 2013 |
| SPG20 | B | L2 | Morris 2010 panel, but lower frequency than SFRP1 |
| RPRM | C | not-detected | Single-study support, low prior |

Scoring rule: "recovered" if any CpG in that gene's promoter region (TSS1500/TSS200/1stExon/5'UTR) appears in the target pool (top-50 for L1 test, Layer-2 pool for L2 test).

### Layer taxonomy (revised 2026-04-23)

**Layer 1 (tissue-wide biomarker)**:
`FDR < 0.05 AND disc_auc > 0.85 AND pct_consistent > 0.8`

**Layer 2 (subgroup-driven)**:
`FDR < 0.05 AND NOT Layer 1 AND pct_responder_at_|Δβ|>0.2 ≥ 0.10`

- Lower bound 0.10 justified: in null, P(|Δβ|>0.2) per pair ≈ 1% (given SD_Δβ ≈ 0.08), so 10% responder rate is ~5× null expectation
- No upper bound: if pct_responder gets high enough, mean Δβ naturally crosses Layer 1 gate

**Layer 3 (weak pervasive)**:
`FDR < 0.05 AND NOT Layer 1 AND NOT Layer 2`

Captures FDR-significant CpGs that show neither strong tissue-wide signal nor identifiable responder subgroup — likely cell-composition shifts, field effects, or weak widespread changes.

### Fallback plan (Atlas framing — also the primary narrative post-pivot to methodological framing)

If strict biomarker validation fails, reuse all Q1/Q3/Q5/Q9 outputs plus layer taxonomy. Narrative restructures as:
- Layer 1 (tissue-wide): the gate-passing population-wide strong signals
- Layer 2 (subgroup-driven): the gate-miss-but-subgroup-signal CpGs
- Layer 3 (weak pervasive): FDR-sig but small/diffuse
- VHL becomes the canonical Layer-2 illustration (correct biological identification)
- GSEA results per layer
- Discussion of cell-composition confounding and tissue-level interpretation

No new data, no new models — only narrative restructuring.

### Pre-registration status

- Primary Q1 test + BH FDR: **unchanged** (plan.md §5 still holds).
- Gate thresholds (0.05 / 0.1 / 0.8), panel size (top 50), one-CpG-per-gene, V1 list, V2 CV setup, baseline sampling parameters (50 random × 30 reps), and Success/Partial/Failure thresholds: **all fixed here before running** (this section acts as pre-registration for the amendment).
- All new analyses are declared as post-hoc / exploratory in Methods.

### Implementation deliverables

- `11_hit_classification.R`: Gate + AUC ranking + top-50 selection (also outputs full Layer 1/2/3 sets for Atlas fallback).
- `12_biomarker_validation.R`: V1 + V2 + V3.
- `13_biomarker_report.R`: figures/tables for paper (to be written after V-results).
- Outputs under `hitclass_out/`, `validation_out/`.

---

## 7. Decision log

Append-only. Newest entries at top. Each entry: date, decision, one-line reason.

- **2026-04-23** — **Pivoted again** to **methodological case study framing** (3-layer decomposition of KIRC methylation). *Reason:* biomarker discovery validation (§5b V2) was degenerate — top-50 panel AUC = 0.99996 indistinguishable from random 50-CpG baselines (all ≥ 0.99), confirming cancer tumor-vs-normal classification is trivially easy at scale and cannot support strong biomarker claims. BUT V1 lenient gave 12/20 literature biomarkers with ≥1 CpG in Layer 1 pool, 15/20 with ≥1 CpG in Layer 2 pool, which is a strong concurrent validity finding for the layer framework itself. The 3-layer taxonomy becomes the contribution; KIRC is the case study. 10-page report writeup starts from this angle.
- **2026-04-23** — **Revised Gate criteria from Option A to Option B** (empirically justified). Gate changed from `FDR + |Δβ|>0.1 + pct_consistent>0.8` to `FDR + disc_auc>0.85 + pct_consistent>0.8`. Rationale: (1) β is bounded; fixed |Δβ| threshold is scale-dependent. (2) AUC is scale-invariant and directly measures discrimination. (3) Empirical check (`14_consistency_vs_auc.R`) confirmed disc_auc>0.85 and pct_consistent>0.8 are correlated (r=0.9) but not redundant — 29K CpGs pass pct_consistent but fail AUC (high within-patient consistency but tiny magnitude — clinically useless). Layer 2 lower bound changed from [20%, 50%] arbitrary range to pct_responder ≥ 0.10 (justified as ~5× the null P(|Δβ|>0.2) rate). Upper bound removed (high pct_responder naturally graduates to Layer 1 via mean Δβ). *Under Option B:* Layer 1 = 36,959 (≈ same), Layer 2 = 61,899 (larger, more inclusive), Layer 3 = 123,528.
- **2026-04-23** — **Pivoted primary framing to biomarker discovery** (§5b amendment); **Atlas (§5) retained as fallback**. Added 3-stage pipeline (Gate → AUC Rank → 3-way Validation) with pre-specified Success/Partial/Failure decision rules, 20-biomarker literature reference list with expected layers, and permutation baseline comparison. *Reason:* Step 3 discussion exposed the gap between "222K FDR hits" (statistical) and "tissue-wide biomarker" (clinical). User chose biomarker-discovery over atlas for primary narrative, with a clear fallback if validation criteria fail. All thresholds/panels/reference sets fixed before running to avoid post-hoc tuning.
- **2026-04-22** — **Locked in Step 2 research-question structure and analysis plan**. Three-layer narrative: Q1 (primary) paired tumor-vs-normal differential methylation; Q3 (secondary) within-tumor CpG associations with tumor grade; Q5 (interpretation) GSEA/ORA on Q1 and Q3 hits. Q9 (positive control) pre-registers three VHL/HIF-axis assertions as pipeline sanity checks, *not* part of MCC family. Excluded stage (too correlated with grade), sex×tissue interaction (complexity/space), unsupervised subtyping (not in course), survival analysis (Cox not in course). *Reason:* user approved this structure and scope; sex-interaction question deliberately dropped to keep report coherent. Full plan in Section 5.
- **2026-04-22** — **Picked TCGA-KIRC as the first (and likely only) cohort**, starting with the full 319 tumor + 160 normal = 479 samples → 160 paired patients. *Reason:* Largest paired n in TCGA (nearly 2× the runner-up BRCA); ccRCC single-histology keeps heterogeneity low; well-established VHL/HIF methylation biology provides sanity-check signals; 6.4 GB is a manageable download.
- **2026-04-22** — **Also excluded somatic mutation and CNV**, committing fully to methylation only (Option A). *Reason:* Neither somatic mutation's canonical methods (MutSigCV, dNdScv) nor CNV's (GISTIC) are taught in the course. Using them would just re-run the same per-feature regression + MCC + GSEA framework on different features — no new course-method coverage, extra explanation overhead. Methylation already saturates the course methods applicable under TCGA's data structure.
- **2026-04-22** — **Locked in DNA methylation as the omics type** for the project. *Reason:* Primary logic is elimination among DNA-level data. (a) Germline SNP (the course's natural GWAS substrate) is TCGA-controlled-access → infeasible. (b) RNA-Seq / miRNA are transcriptomic and never appeared in course → off-theme. This leaves methylation, somatic mutation, CNV. Methylation chosen among the remainders because of HW5 question (c) signal, assay maturity, paired design, and 1:1 method transfer from SNP→CpG. Full rationale in Section 9 item #1.
- **2026-04-22** — Adopted workflow refinements: (i) pre-register each (dataset, question) pair before running; (ii) disclose all tried combos in final report; (iii) selection criterion is effect plausibility × story coherence × methodological cleanliness, not min-p; (iv) added Step 0 env-setup; (v) added reproducibility guard in Step 3 (seed + sessionInfo + saved intermediates). *Reason:* user accepted the anti-p-hacking and reproducibility adjustments proposed during workflow discussion.
- **2026-04-22** — Tightened Step 1 exit criteria to check *assay-matched* paired samples (tumor + normal both run on same assay), not just "has normal", plus array-version split (450K vs 27K not mixed) and multi-omics patient-level intersection. *Reason:* user asked about data coverage uniformity → TCGA is sparse at three levels (project / patient / sample), so triage has to verify intersections, not assume them.
- **2026-04-22** — Deferred discussion on "how to better formulate research questions" (Section 1a of prior draft) until after user reviews TCGA data structure. *Reason:* user wants to understand the data first before committing to question-formulation principles.
- **2026-04-22** — Adopted iterative 4-step workflow (triage → question → execute → audit), looping over multiple dataset/question pairs before committing. *Reason:* hedges against bad-dataset surprises; matches user's request to try several angles before locking one in.
- **2026-04-22** — Drafted `lecture_methods.md` from course lecture PDFs as the in-scope reference. *Reason:* final project must be anchored in course material.
- **2026-04-22** — Final project priority ordering: analysis → results → proposal. HW5 submission will be a snapshot of Steps 1–2 of this plan. *Reason:* user explicitly requested this order.

## 8. Open questions for the user

Things we need input on before we can proceed. Cleared as they are answered.

1. ~~**DNA methylation as default, or open to other omics?**~~ — **RESOLVED 2026-04-22**: methylation locked in. Rationale in Section 9 item #1.
2. **Compute budget.** Downloading 450K methylation for a whole cancer cohort is ~1–5 GB per cancer. OK to pull 2–3 cohorts? Any storage limit on this laptop?
3. **R environment.** Is TCGAbiolinks already installed? Bioconductor version? (If not, Step 1 starts with environment setup.)
4. **How autonomous?** User said "automatically try multiple datasets". Confirming: it's OK for the assistant to download data and run analyses without per-step approval, as long as (a) results are real, (b) every decision lands in this file, (c) user approves the final pick before report-writing?
5. **Question-formulation principles** — deferred. User wants to discuss *how* to formulate good research questions after reviewing TCGA data structure. Revisit before Step 2.

## 9. Paper write-up checklist — points the final report must address

*Running list. Append as we go. These are things we have explicitly decided are worth mentioning in the final report / proposal, with enough context that a future writer (us or a fresh session) knows why each is on the list.*

Each entry: **topic → one-line claim to make → where it lands in the paper** (Intro / Methods / Results / Discussion).

1. **Why DNA methylation?** → Argue by elimination, not by methylation's positive merits alone.
    - **Primary reason**: this is a *Statistical Genetic Modeling* course whose methods target DNA-level germline variation. Among TCGA's data categories:
        - Germline SNP array (the natural GWAS substrate) is **controlled access via dbGaP** — not feasible for a course project.
        - RNA-Seq and miRNA-Seq are **transcriptomic, not DNA-level**, and the course never used them — off-theme.
        - Somatic mutation and CNV are DNA-level and open access, but **not taught in the course** — their canonical methods (MutSigCV, dNdScv, GISTIC, mutual-exclusivity tests) are all absent from the syllabus. Using them would only mean re-running the same per-feature regression + MCC + GSEA on a different feature type — no additional course-method content, plus extra explanation burden.
        - Methylation is what survives: DNA-level, open access, course-aligned.
    - **Secondary reasons** supporting methylation specifically:
        - HW5 question (c) explicitly names DNA methylation preprocessing → instructor-signaled expectation.
        - 450K methylation is TCGA's largest and most standardized open-access DNA assay.
        - Paired tumor / adjacent-normal design is natural for methylation and boosts power.
    - **Method coverage claim**: methylation uses essentially all course methods that are applicable to TCGA's data structure. Directly usable: Lec 6 logistic/linear regression, Lec 6 PCA stratification correction, Lec 8 per-feature scan + Manhattan/Q-Q, Lec 9 Bonferroni/BH/permutation, Lec 10 GSEA/ORA and conditional analysis. By analogy: Lec 3 LD → CpG co-methylation; Lec 10 LD clumping → DMR calling. Unapplicable (for any TCGA choice, not just methylation): Lec 2 HWE, Lec 4–5 linkage, Lec 7 TDT — these require genotype / pedigree data that TCGA does not provide.
    - → **Intro / Methods intro / Discussion (re: which course methods were & weren't applicable to TCGA)**

2. **Why VHL/HIF as positive control?** → ccRCC is the archetypal "VHL/HIF disease" — ~60–80% of sporadic ccRCC cases inactivate VHL, and VHL promoter hypermethylation (in ~10–20% of cases) is a classic epigenetic silencing mechanism documented since the 1990s. Downstream HIF targets (VEGFA, CA9, etc.) also show methylation changes. Using VHL/HIF as a pipeline sanity check means: if our methylation pipeline cannot recover this textbook KIRC signal, the pipeline is broken. Alternatives considered (CDKN2A — not KIRC-specific; MGMT — GBM-specific; LINE-1 — sparse on 450K) were weaker. Pre-registered assertions, treated as methodological QC, not reported as "findings". → **Methods (pipeline validation) / Results (sanity check subsection)**

3. **Anti-p-hacking disclosure**: All pre-registered research questions (Q1, Q3, Q5, Q9) and excluded candidates (Q2/Q4/Q6/Q7/Q8) must be named in the report. Section §5 of `plan.md` is the authoritative pre-registration record. → **Methods (analysis-plan transparency)**

4. **Chosen preprocessing pipeline** (decision-log, not a result): TCGAbiolinks β-matrix → drop chrX/Y probes → drop cross-reactive (Chen 2013) → drop SNP-overlap → drop probes with >20% missingness → keep complete-case → β→M with [0.001, 0.999] clipping. The chrX/Y drop is **driven by** our own PC diagnostic finding that gender dominates PC2 (R²=0.94), not just by convention. → **Methods (preprocessing)**

5. **Sensitivity analyses we committed to run** (pre-registered): (a) Q1 without covariate adjustment → quantify inflation from sex/age/TSS; (b) unpaired-samples version of Q1 → hit-list overlap; (c) exclude any CpG with >0% missingness → confirm robustness. → **Results (sensitivity analyses subsection)**

6. **Methods amendment from pre-reg**: The pre-registered Q1 used `limma + duplicateCorrelation` with patient as random effect and covariates. In practice this was infeasible on 320K probes × 320 samples (computational; killed after 20 min). Switched to the canonical **within-pair differences + moderated one-sample t** approach, which is mathematically equivalent for balanced paired design and absorbs patient-level covariates implicitly. Plate-level confounding addressed in sensitivity S2 (same-plate-only subset, 155/160 pairs → 100% hit overlap with primary). Must be disclosed in Methods. → **Methods amendment statement**

7. **λ_GC interpretation**: observed λ = 36.77 (Q1) and 6.58 (Q3) are far above conventional GWAS QC thresholds but **are expected and not a concern** here. Proof: (a) Q-Q plot shape shows curve departing y=x immediately and smoothly, not a tail spike; (b) permutation null (sign-flipped tumor/normal labels) yields mean λ = 1.009 ± 0.39 — methodology is well calibrated; (c) 69% of CpGs being differentially methylated in cancer matches published expectation for cancer epigenomes. Report should explicitly explain this so reviewers don't reflexively flag it. → **Methods (QC discussion) / Results**

8. **VHL mixed direction (positive control nuance)**: 2/7 VHL promoter CpGs strongly hypermethylated in tumor (including cg13672843 at 1stExon, large effect), 2/7 significantly hypomethylated (TSS1500 CpGs), 3/7 n.s. Don't over-simplify to "VHL silenced by methylation" — report the bidirectional pattern. Pre-registered assertions still all PASS. → **Results (positive-control subsection) / Discussion**

9. **Biological-significance threshold**: report both FDR counts AND effect-size-filtered counts (FDR<0.05 AND |Δβ|>0.1). For Q1: 222K → 46K strong hits. For Q3: 99K → 19K strong hits. This avoids "statistical significance but biologically trivial" slip-up. → **Results (main tables)**

10. **Biomarker operational definition**: Report explicitly that our "tissue-wide biomarker" is defined as FDR<0.05 AND |Δβ|>0.1 AND pct_consistent>0.8, then ranked by per-CpG unpaired AUC. Explain why this is stricter than pure FDR filtering (statistical significance is not biomarker utility). → **Methods (operational definition), Results (stage-by-stage attrition)**

11. **Hard baselines for classifier AUC**: Our top-50 AUC must be compared against (L1) random 50 from FDR<0.05 and (L2) random 50 from FDR<0.05 ∩ |Δβ|>0.1, each with 30 replicates. This is an "ablation" showing which filter contributes value. Without this baseline, any AUC claim would be weak. → **Results (V2 validation), Discussion (marginal value of selection)**

12. **Concurrent validation via literature biomarker list**: 20 published KIRC methylation biomarkers with pre-specified expected layer (L1/L2/not-detected). Recovery rate serves as concurrent validity check. Hypergeometric enrichment test for overlap significance. → **Results (V1 validation), Methods (literature list with references)**

13. **Fallback plan disclosed**: If biomarker validation fails pre-specified thresholds, reframe as descriptive atlas (3-layer methylation landscape). Disclose this contingency even if not used — shows pre-registered honesty. → **Methods (framing choice), possibly Discussion**

14. **VHL as a Layer-2 illustration (not Layer-1 failure)**: VHL is clinically "tissue-wide" (most ccRCC have VHL-pathway disruption) but methylation-driven in only ~15–20% → correctly lands in Layer 2. Report this explicitly as evidence that the classification correctly identifies heterogeneous vs homogeneous epigenetic patterns, not as a failure to recover a classic marker. → **Results (Q9 positive control discussion)**

*(more items to be added as we discuss)*

---
