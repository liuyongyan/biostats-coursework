# P8139 Final Project — Working Plan

> Living document. Updated every time we make a substantive decision.
> If the conversation session dies, a new session should be able to read this file + `lecture_methods.md` and pick up where we left off.

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
- For each: total samples, tumor count, paired-normal count, metadata completeness, feature-level missingness.
- **Exit criterion**: at least one combo with (a) ≥30 paired samples *or* ≥100 tumor+ ≥20 independent normal, (b) metadata fields needed for the planned question, (c) manageable download size.

### Step 2 — Research question + analysis plan
- 1–3 concrete research questions per surviving dataset.
- Each maps to a primary test and a covariate set.
- Power / feasibility sanity check (not a formal power calc; just "is n large enough for this effect scale?").
- **Exit criterion**: question is a sentence a reviewer would understand; plan cites a specific method from `lecture_methods.md`; multiple-comparison strategy specified up front.

### Step 3 — Execute
- Reproducible R script; seed set; `sessionInfo()` captured.
- Run the primary test, then diagnostics (Q-Q, volcano, PCA of residuals), then secondary/sensitivity analyses.
- **Exit criterion**: code runs end-to-end from raw TCGA download to result tables/plots without manual intervention.

### Step 4 — Audit
- Assumption checks, effect-size sanity (are the top CpGs/genes biologically plausible?), covariate balance, negative controls if available.
- Compare to at least one sanity baseline (e.g. permuted labels should give flat Q-Q).
- Decide: **accept** / **revise plan** (back to Step 2) / **drop dataset**.
- **Exit criterion**: written audit note in the decision log below.

Loop Step 1→4 over multiple (dataset, question) pairs. At the end, pick the strongest one with the user and write it up.

## 4. Candidate datasets — living shortlist

*Filled in during Step 1. Format: one row per (cancer, omics) combo.*

| ID | Cancer | Omics | N tumor | N normal | N paired | Notes | Status |
|---|---|---|---|---|---|---|---|
| *(pending Step 1)* | | | | | | | |

## 5. Candidate research questions — living shortlist

*Filled in during Step 2. Format: one row per (dataset, question) pair.*

| ID | Dataset | Question | Primary test | Covariates | MCC | Status |
|---|---|---|---|---|---|---|
| *(pending Step 2)* | | | | | | |

## 6. Current status

- **Last updated**: 2026-04-22
- **Active task**: briefing on TCGA + aligning on workflow (Tasks #3, #4)
- **Next up**: Step 1 — pick initial (cancer × omics) shortlist for download & EDA (Task #5)
- **Blocked on**: user confirmation of workflow + initial dataset preference (broad net vs focused on DNA methylation)

## 7. Decision log

Append-only. Newest entries at top. Each entry: date, decision, one-line reason.

- **2026-04-22** — Adopted iterative 4-step workflow (triage → question → execute → audit), looping over multiple dataset/question pairs before committing. *Reason:* hedges against bad-dataset surprises; matches user's request to try several angles before locking one in.
- **2026-04-22** — Drafted `lecture_methods.md` from course lecture PDFs as the in-scope reference. *Reason:* final project must be anchored in course material.
- **2026-04-22** — Final project priority ordering: analysis → results → proposal. HW5 submission will be a snapshot of Steps 1–2 of this plan. *Reason:* user explicitly requested this order.

## 8. Open questions for the user

Things we need input on before we can proceed. Cleared as they are answered.

1. **DNA methylation as default, or open to other omics?** HW5 question (c) specifically mentions methylation preprocessing, which suggests the instructor expects methylation. Does the user want to follow that steer or consider expression / mutation / multi-omics too?
2. **Compute budget.** Downloading 450K methylation for a whole cancer cohort is ~1–5 GB per cancer. OK to pull 2–3 cohorts? Any storage limit on this laptop?
3. **R environment.** Is TCGAbiolinks already installed? Bioconductor version? (If not, Step 1 starts with environment setup.)
4. **How autonomous?** User said "automatically try multiple datasets". Confirming: it's OK for the assistant to download data and run analyses without per-step approval, as long as (a) results are real, (b) every decision lands in this file, (c) user approves the final pick before report-writing?

---
