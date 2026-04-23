# P8131 Homework 8 Grading Rubric

**Course:** P8131 Biostatistical Methods II (Spring 2026)
**Scope:** Problem 1, parts (a)-(c) — GEE and GLMM analysis of the HEALTH longitudinal RCT dataset.

This is a **principle-based rubric**. The TA reviewer must verify correctness of the student's numbers, code, and reasoning **independently against the source data, the R output, and the lecture slides** — this document does not encode a "correct" answer. Pattern-matching a student's answer against a pre-computed one would itself risk propagating errors.

---

## General Expectations (applies to all parts)

- **Reproducibility.** Code should run top-to-bottom against `HEALTH.xlsx` without manual intervention. Variable recodings (factor levels, ordering, baseline-covariate construction) must be visible in the code.
- **Written interpretation is required** for every part. A table of coefficients with no prose is not sufficient — the problem literally says "Interpret …" in each sub-part.
- **Factor handling.** `HEALTH`, `TRT`, `AGEGROUP` are character/categorical. Student should convert to factors with an intentional reference level, or explicitly recode to 0/1 for the binary response. The reference level they pick governs the sign/interpretation of every coefficient — the reviewer should read the interpretation *in light of* the reference level the student actually set, not in light of an imagined one.
- **Data sort.** For `gee()`, the data must be sorted by subject ID (Lecture 19, slide 13 explicitly notes: "Data need to be sorted by subject!"). If the student skips this, the working correlation and robust SEs can be wrong — flag it even if numbers happen to look reasonable.
- **Package/function expectations.** The lectures use `gee::gee()` for GEE (Lec 19 slides 13-20) and `lme4::glmer()` for GLMM (Lec 20 slides 15-17). Using `geepack::geeglm()` or `glmmTMB` is not wrong per se, but the student should be able to justify the swap; if they do swap, the TA should verify argument choices (e.g. `corstr`, `family`) line up with lecture intent.

---

## Part (a): Bivariate cross-sectional analysis at randomization

**(i) What the student must answer/compute**

- A direct, explicit answer to: *Is there an association between randomized group (TRT) and self-rated HEALTH at TIME == 1?* This is a yes/no question backed by a test statistic, p-value, and an effect-size summary (odds ratio, risk difference, or the 2×2 proportions).
- Interpretation must mention the RCT context — because randomization should (in expectation) balance groups at baseline, the student should note what the finding implies about successful randomization / baseline balance.

**(ii) Key steps that MUST be shown**

1. Subset the data to `TIME == 1` only.
2. State how many subjects appear at TIME==1 and confirm it matches n = 80 (flag if not — possible data issue worth a sentence).
3. Produce the 2×2 cross-tabulation of TRT × HEALTH.
4. Perform an appropriate test of association for a 2×2 table (e.g., Pearson chi-square, Fisher's exact, or a logistic/Wald test). Any of these is acceptable; the student should briefly say why their choice is appropriate (expected counts, sample size).
5. Report the test statistic, degrees of freedom (if applicable), and p-value.
6. Interpret in plain language at a stated significance level.

**(iii) Method/R function guidance**

- This part is *pre-longitudinal*. It is basic categorical data analysis — a `table()`, `chisq.test()`, `fisher.test()`, or a simple `glm(..., family=binomial)` are all appropriate. It does not draw on the Lec 19/20 machinery.
- The reviewer should verify that the test actually corresponds to the hypothesis being stated in the interpretation (e.g., if the student reports an OR from logistic regression they should not then cite a chi-square p-value from a different test without explanation).

**(iv) Common pitfalls to watch for**

- Using the full dataset (all four visits) instead of TIME==1 only. This is **wrong** — part (a) is explicitly cross-sectional at randomization.
- Using TIME==0 instead of TIME==1. The codebook says `1 = randomization`; TIME==0 does not exist in this dataset. (This is a different convention from the `respiratory` example in Lecture 19, where baseline is month 0 — do not let a student import that convention blindly.)
- Confusing "no significant association" with "treatment has no effect" — at baseline there has been *no treatment yet*, so significance here speaks to randomization balance, not to efficacy.
- Omitting the interpretation sentence entirely.
- Reporting only a p-value with no effect-size descriptor.

---

## Part (b): Longitudinal GEE with unstructured working correlation

**(i) What the student must answer/compute**

- A fitted GEE on the **follow-up only** dataset (TIME ∈ {2, 3, 4}) with response = HEALTH, predictors = (baseline health, TRT, month post-randomization, AGEGROUP). The problem specifies **unstructured** working correlation.
- A written interpretation of **each** regression coefficient (not just the treatment effect) on the appropriate scale (log-odds or odds-ratio, since HEALTH is binary → logistic link is expected). The interpretation must be explicitly *population-averaged / marginal*, per Lecture 19 slide 8.
- A statement about whether the treatment has a detectable effect on self-rated health over follow-up, and whether age group modifies or confounds this.

**(ii) Key steps that MUST be shown**

1. Construct a **baseline-health covariate** from TIME==1 and merge/broadcast it onto the follow-up rows. This is exactly the pattern used in Lecture 19 slide 13 (the `resp$baseline <- rep(subset(..., month == "0")$status, rep(4, 111))` step). The student must verify the ordering of the broadcast (sort by ID, confirm length matches).
2. Subset to TIME ∈ {2, 3, 4} (the problem says "across all study follow-up visits (but not at randomization)").
3. Convert HEALTH to a 0/1 numeric (or use a factor and know which level is "success"). The sign of every coefficient depends on this — check that the student's interpretation matches the coding they used.
4. Define a **month-post-randomization** numeric covariate. The codebook maps TIME 2/3/4 → 3/6/12 months. Using TIME as 2/3/4 is **not** the same as using months 3/6/12 — the problem says "month post randomization." Flag any student who feeds TIME (1-4 ordinal) in where months are asked for, unless they explicitly justify it.
5. Sort the data by subject ID before calling `gee()`.
6. Fit with logit link (binomial family) and `corstr = "unstructured"`. Refer to Lecture 19 slide 15 bottom (`resp_gee3 <- gee(..., corstr = "unstructured", scale.fix = FALSE)`) for the exact argument style.
7. Print `summary()` of the fit, including the estimated working correlation matrix.
8. Interpret using **robust (sandwich) SEs and z-values**, not the naive SEs (Lecture 19 slide 10 explains why; the `gee` output shows both columns).
9. Report at minimum: point estimate, robust SE, robust z or p-value, and either the coefficient (log-OR) or the exponentiated OR — whichever the student chooses, stay consistent.

**(iii) Method/R function guidance**

- The canonical template is **Lecture 19, slides 13-16** (respiratory binary-GEE example). The student's code should mirror that structure: load → construct baseline covariate → subset to follow-ups → sort by id → `gee(... , family="binomial", id=subject, corstr="unstructured", scale.fix=FALSE)`.
- The reviewer should compare the student's `gee()` call argument-by-argument against slide 15 and verify that any deviation (e.g. `scale.fix=TRUE`) is defensible.
- Interpretation language should explicitly invoke "population-average" / "marginal" / "averaged over subjects" consistent with Lecture 19 slide 8.

**(iv) Common pitfalls to watch for**

- **Including TIME==1 in the fit.** Doing so double-counts the baseline observation (once as outcome, once as covariate) and contradicts the prompt.
- **Omitting the baseline-health covariate**, or including only follow-up predictors without it.
- **Using naive SEs for inference.** The whole point of GEE is that the working correlation may be misspecified but the sandwich SE is still consistent. Students who cite naive p-values should be flagged.
- **Wrong time metric.** Using TIME codes 2/3/4 as if they were months. Months are 3, 6, 12 — non-equally-spaced.
- **AGEGROUP as numeric.** AGEGROUP has three ordered levels but is categorical; treating it as 1/2/3 numeric imposes a linearity constraint that was not asked for. A factor with two dummy coefficients is expected unless the student explicitly argues for a linear trend.
- **Misinterpreting a categorical coefficient** without stating the reference level.
- **Interpreting subject-specific when it is marginal.** E.g., "for a given woman, the odds of good health increase by …" — this is the GLMM interpretation, not GEE. See Lecture 20 slide 13.
- **Forgetting to sort by ID** before `gee()` (Lecture 19 slide 13 warning).
- **Failing to address the question literally asked.** The prompt says "describe the relationship … interpret your results" — a student who only prints the summary and writes "the treatment coefficient is -0.3" has not interpreted.

---

## Part (c): GLMM with subject-specific random intercept

**(i) What the student must answer/compute**

- A fitted generalized linear mixed effects model with a **random intercept per subject**, same fixed-effects structure as in part (b) (same response, same four covariates, same follow-up subset).
- Interpretation of each fixed-effect estimate on the **subject-specific / conditional** scale.
- A direct comparison with the GEE estimates from (b): do the signs match? Do the magnitudes differ in the direction expected from theory (subject-specific log-odds coefficients are typically *larger in magnitude* than marginal ones for logistic models with random intercepts)? The student must discuss **why** the two interpretations differ, not just that they differ.
- Report the estimated random-intercept variance and discuss what it tells us about between-subject heterogeneity.

**(ii) Key steps that MUST be shown**

1. Same data construction as part (b) (baseline covariate, follow-up subset, months coding). If the student re-does it, verify it is consistent with (b); if they reuse the (b) data frame, verify it is the right one.
2. Fit a binomial GLMM with `(1 | ID)` random intercept and the same fixed effects. The canonical template is Lecture 20 slide 15: `glmer(response ~ fixed_effects + (1 | subject), family = 'poisson', data = ...)` — here `family = binomial` is appropriate since HEALTH is binary. Verify `family="binomial"` (not poisson) is used.
3. Print `summary()` and extract: fixed-effects table, random-effect variance/SD, number of observations and groups.
4. State on the subject-specific scale what each coefficient means (per Lecture 20 slide 12: "change in log odds of a positive response per unit change in X_ij, **for the same subject**").
5. Contrast marginal vs. conditional interpretation explicitly, referring to the shrinkage/attenuation phenomenon (Lecture 20 slide 13: E(Y_i) ≠ g^-1(X_i^T β_*) unless the link is linear).
6. Compare estimates in a table or side-by-side narrative between (b) and (c).

**(iii) Method/R function guidance**

- Canonical template: **Lecture 20, slides 15-17** (`ep.GLMM1 <- glmer(...)` with `(1 | subject)`). Student's `glmer()` call should follow that style.
- Required library: `lme4` (loaded as `library(lme4)` on slide 15).
- The estimation method for binary GLMM via `glmer` is Laplace approximation by default — students do not need to change this, but if they use `nAGQ > 1` they should say why.
- Interpretation reference: Lecture 20, slide 12 ("β_2 measures the change in log odds … for the same subject") and slide 13 (contrast with marginal).

**(iv) Common pitfalls to watch for**

- **Using `family = "poisson"`** (copying Lecture 20's epilepsy example verbatim). HEALTH is binary → must be `family = binomial`.
- **Including random slopes the problem did not ask for.** The prompt says "subject-specific random intercepts" — keep it to `(1 | ID)`. A student may optionally compare a richer model, but the primary model must be random-intercept-only.
- **Interpreting GLMM coefficients as if they were population-averaged.** This is the central conceptual point of Lecture 20 — getting it wrong here is the single biggest conceptual error for this problem.
- **Claiming the difference between (b) and (c) is due to "different estimators."** It is not just estimation difference; the *estimands* differ (marginal vs. conditional). Lecture 20 slide 13 is the touchstone. Student should reference that the two quantities are genuinely different functionals of the data-generating process, not two estimates of the same thing.
- **Ignoring the random-effect variance.** A near-zero variance says subject-level heterogeneity is small and GEE/GLMM estimates should be close; a large variance predicts meaningful attenuation of GEE relative to GLMM. The student should connect the estimated variance to their observed magnitude difference.
- **Convergence warnings ignored.** If `glmer` emits a warning (singular fit, failed to converge), the student should acknowledge it and discuss.
- **Reusing the GEE interpretation verbatim.** If the student simply pastes the part-(b) interpretation under part (c), that is a failure to answer "How are the interpretations different from the GEE model?"

---

## Cross-cutting checklist for the TA reviewer

Before assigning credit, verify independently:

- [ ] The data subsetting for each part actually matches the prompt (a: TIME==1; b,c: TIME ∈ {2,3,4}).
- [ ] The baseline-health covariate is constructed correctly (one value per subject, broadcast to follow-up rows, merged by ID not by row position).
- [ ] The HEALTH response coding (which level = 1) is consistent with the direction of interpretation.
- [ ] For GEE: robust SEs are used; working correlation is unstructured; data sorted by ID.
- [ ] For GLMM: `family=binomial`, `(1|ID)`, random-effect variance reported.
- [ ] Written interpretations exist for every sub-part and are on the correct scale (marginal in b, conditional in c).
- [ ] Numbers reported in prose match what appears in the R output (no transcription errors).
- [ ] Code is reproducible end-to-end.

**Reviewer note:** Do not treat this rubric's wording as "the right answer." Re-run the student's code, inspect the output, and judge the interpretation against the slides and the output — not against an imagined gold-standard result.
