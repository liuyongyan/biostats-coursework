# P8139 Statistical Genetic Modeling — In-Scope Methods

Reference list of methods covered in the course. The final project should lean on these; going outside is allowed but must be justified and should not dominate the analysis.

Source: `~/Desktop/Columbia - Biostatistics/Statistical Genetic Modeling/Lecture{1..10}*.pdf`

---

## 1. Population genetics foundations (Lec 1–2)
- **Allele / genotype frequency estimation** — direct counting.
- **Hardy-Weinberg equilibrium (HWE) test** — χ² or Fisher's exact on observed vs expected p², 2pq, q². Used for QC (filter SNPs/CpGs violating HWE in controls) and as a sanity check for genotyping error.
- **Genetic model specification** — additive / dominant / recessive / codominant / multiplicative encodings of a variant. Drives the design matrix in regression.
- **Penetrance modeling** — f = P(affected | genotype); relevant to parametric linkage.

## 2. Recombination & linkage disequilibrium (Lec 3)
- **Recombination fraction θ** — estimated from meioses.
- **LD measures** — D, D′, r². r² is the effective sample-size-relevant measure in association.
- **Haplotype block / LD structure** — used to prune correlated SNPs (clumping) and define independent signals.

## 3. Linkage analysis (Lec 4–5)
- **Parametric (model-based) LOD score** — LOD(θ) = log₁₀ L(θ)/L(½); threshold 3 for significance, −2 for exclusion. Two-point and multipoint variants.
- **Non-parametric linkage** — affected sib-pair (ASP) IBD sharing test; Kong–Cox nonparametric LOD.
- **Haseman-Elston regression** — regress squared sib trait difference on IBD sharing for quantitative traits.

*Note:* Linkage methods require pedigree/family data, which TCGA does not have. These are unlikely to appear in our project but are in scope conceptually.

## 4. Association testing (Lec 6)
- **χ² / Fisher's exact test on allele or genotype counts** — simplest case-control test.
- **Logistic regression** — logit P(case | G, covariates) = α + βG + γᵀZ. Encodes genetic model via G. Standard for binary outcome.
- **Linear regression** — for quantitative trait. Extends to methylation β/M-values, expression, etc.
- **Odds ratio / effect size** — reported alongside p-values.

## 5. Population stratification correction (Lec 6)
- **Principal component analysis (PCA)** on genotype matrix — include top PCs as covariates to adjust for ancestry. Directly applicable to methylation/expression matrices for adjusting latent structure.
- **Genomic control λ_GC** — rescale test statistics by median χ²/0.456; diagnostic Q-Q plot.
- **Stratified analysis** — split by subgroup or add covariate.

## 6. Family-based association (Lec 7)
- **TDT (McNemar)** — χ² = (b−c)²/(b+c) on transmitted vs untransmitted alleles in trios.
- **FBAT** — general framework for within-family transmission tests.
- **Quantitative TDT** — trait-stratified transmission.

*Note:* Requires trios/families, not in TCGA. Out of scope for our project unless we use an external dataset.

## 7. GWAS methodology (Lec 8)
- **Sample QC** — missingness, relatedness (kinship/IBD), stratification outliers.
- **Variant QC** — MAF threshold, call rate, HWE in controls.
- **Per-locus test** — logistic/linear regression on each variant (or CpG, or gene).
- **Genome-wide significance threshold** — 5×10⁻⁸ for SNPs; for methylation or expression the effective number of independent tests is smaller, but Bonferroni / FDR still apply.
- **Manhattan plot** — −log₁₀(p) vs genomic coordinate.
- **Q-Q plot** — observed vs expected p-values; diagnose inflation.

## 8. Permutation & multiple testing (Lec 9)
- **Permutation test** — label-shuffling to generate empirical null; empirical p = (#extreme + 1)/(B + 1). Can be applied per-locus or to a global test statistic.
- **Bonferroni FWER** — α/m. Conservative but distribution-free.
- **maxT / minP FWER** — permutation-based, accounts for correlation.
- **Benjamini-Hochberg FDR** — adjusted p = min over k≥i of m·p_(k)/k. Standard for high-dimensional screens.
- **Storey q-value** — FDR with estimated π₀.

## 9. Post-GWAS / downstream interpretation (Lec 10)
- **LD clumping** — keep lead variant; remove others with r² > threshold. Analog for omics: collapse correlated features within a gene/region.
- **Fine mapping / credible sets** — variants with cumulative posterior ≥ 99%.
- **Conditional / joint analysis** — test secondary signals holding lead fixed.
- **Gene set enrichment analysis (GSEA)** — rank-based permutation test across gene-set members; pairs with the Subramanian 2005 paper assigned in class.
- **Over-representation analysis (ORA)** — hypergeometric test on foreground vs background gene lists.
- **Functional annotation** — map variants/CpGs to genes, promoters, enhancers (ENCODE/GTEx).

---

## Map to TCGA / omics project

| Course method | How it applies to TCGA analysis |
|---|---|
| HWE test | QC filter on methylation probes (control samples), genotype data |
| Additive/dominant genetic models | Encoding SNPs or categorical covariates (e.g., stage) in regression |
| Logistic regression + covariates | Tumor vs normal classification by feature(s) |
| Linear regression per feature | Differential methylation / expression per CpG / gene |
| **Paired analysis** (extension of Lec 6 idea) | Paired t-test / mixed model for tumor–adjacent-normal pairs |
| PCA for stratification | Adjust for batch / ancestry / cell composition; visualize sample structure |
| Bonferroni / BH FDR | Control across 450k CpGs or ~20k genes |
| Permutation test | Empirical significance for global signatures or correlated tests |
| Manhattan + Q-Q plots | Visualize per-feature p-values; check inflation |
| GSEA / ORA | Biological interpretation of significant feature lists |
| LD clumping analog | Collapse correlated CpGs within a region/gene |

---

## Out-of-scope (but may appear in support)
- Cox proportional hazards for survival — not in lectures; may borrow from P8131. Declare if used.
- Mixed-effects models (lme4) — not in lectures; may use for paired/repeated designs. Declare if used.
- Multi-omics integration (iCluster, MOFA) — beyond scope; avoid unless auxiliary.
