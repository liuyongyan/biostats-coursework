# Final report — deferred items

These are improvements identified during the journal-quality review pass on `P8139_FinalReport_yl6107.Rmd`. They were judged worth doing for publication readiness but deferred from the current revision.

---

## TODO #1 — Redesign the literature biomarker benchmark

**Current state.** Table 4 reports recovery against an 18-gene ccRCC TSG list. Hypergeometric enrichment is 1.96× (p = 0.014) for strict silencing — statistically significant but modest in magnitude, and it is the weakest empirical claim in the paper.

**Root cause.** The 18-gene list mixes well-evidenced TSGs (VHL, RASSF1, SFRP1, SFRP2, DKK3, GATA5) with single-cohort or weakly supported entries (NEFH, APAF1, SFRP4, PITX2). The list dilutes the signal.

**Two options to consider:**

1. **Narrow benchmark** — restrict to ~6 high-confidence ccRCC methylation TSGs with multi-cohort or mechanistic support (suggested set: VHL, RASSF1, SFRP1, SFRP2, DKK3, GATA5). Recovery becomes ~5/6, hypergeometric p << 0.001. Strong but accused of cherry-picking.
2. **Tiered benchmark** — expand to ≥30 entries with explicit evidence-strength tiers (Tier A: multi-cohort; Tier B: single-cohort with mechanism; Tier C: single-cohort association only). Report recovery per tier. More defensible but more work.

**Recommendation.** Tier-stratified is publication-stronger. Source candidates from MethylMix output, ELMER pan-cancer atlas, and PubMed search "ccRCC + DNA methylation + TSG".

---

## TODO #2 — Direct cell-composition validation with EpiDISH

**Current state.** The report claims that the silencing filter substitutes for reference-based deconvolution (EpiDISH, CIBERSORT) for the cell-composition confound. The evidence is indirect: immune-pathway enrichment collapses 10 orders of magnitude under the silencing restriction.

**Limitation.** Without a direct comparison, the claim is a plausibility argument, not a demonstration. A reviewer can ask "did you actually run EpiDISH and check that the silencing-filtered list is no longer correlated with infiltrate fractions?"

**Plan.**
1. Run `EpiDISH::epidish` with the centDHSbloodDMC.m reference panel on the 484-aliquot β matrix to estimate per-sample fractions of B / NK / CD4T / CD8T / Mono / Neutro / Eos / Treg / Epi / Fib.
2. For each silencing CpG and for each raw-DM-only CpG, regress β on the inferred infiltrate fractions in tumor samples; compare the distribution of variance-explained.
3. If the silencing-restricted list shows substantially less correlation with infiltrate composition than the raw-DM list, this is direct numerical evidence for the §3.5 claim.

**Cost.** EpiDISH installation + ~10 min runtime. Adds ~1 page of methods/results. Significantly strengthens the most novel claim of the report.

---

## Lower priority (optional)

- Survival association of silencing genes (Cox PH on TCGA-KIRC OS / PFI). Out of P8139 scope (would draw on P8131); add as future-work bullet only.
- Activation pairs (n=2,214) currently mentioned but not analyzed. Either remove all references or write a 1-paragraph oncogene-activation supplementary analysis.
