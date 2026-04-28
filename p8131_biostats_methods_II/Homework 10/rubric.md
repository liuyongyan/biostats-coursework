# P8131 Homework 10 — Grading Rubric (Principle-Based)

This rubric specifies WHAT a correct submission must contain and HOW the TA reviewer should verify it. The rubric does NOT contain pre-computed answers; the reviewer must independently re-derive numerical results and compare against the student's submission.

Course context: Survival analysis is taught in Lectures 21–23 by Prof. Bin Cheng. The teaching conventions to enforce are:
- Censoring is denoted with a trailing `+` (Lec 21, slide 7).
- KM/NA estimators and `survfit`/`Surv`/`summary(KM)` style code are demonstrated in Lec 21, slides 28–30.
- Log-rank test is implemented via `survdiff(Surv(time,cens) ~ group, data=...)` shown in Lec 22, slide 13. The hypothesis form `H0: S_A(t)=S_B(t) for all t` (equivalently `h_A=h_B`) is stated in Lec 22, slide 6.
- Cox PH model is fit via `coxph(Surv(time,delta) ~ ..., data=..., ties='breslow')` and stratified models with `strata(...)`, both demonstrated in Lec 23, slides 26–28. The larynx dataset and stage/age modeling are explicit examples in Lec 23, slides 22–28.

---

## Problem 1 — Log-rank test on small two-group data

**Question type:** Direct, explicit hypothesis test. The student must produce a test statistic, df, p-value, and an explicit Reject / Do-not-reject conclusion at α = 0.05.

**Expected derivation type:** Either (a) by-hand log-rank computation built from a per-event 2×2 risk table, or (b) R-based `survdiff(...)` output, or (c) both. By-hand work is encouraged because the data set is tiny — a hand table is a strong demonstration of understanding.

**Key steps that MUST be shown:**

1. **Data coding.** The student must correctly read the `+` notation as right-censoring. A correct submission constructs vectors `time`, `status` (1 = event, 0 = censored) and `group`. The reviewer should verify that:
   - All five Group-1 times are present with the correct status (the `+` entries censored, the others events).
   - All five Group-2 times are present with the correct status.
   - Censoring indicator convention is internally consistent (e.g., 1=event throughout, matching Lec 21 slide 7 convention `δ_i = 1` for event).

2. **Hypotheses.** Must restate `H0: h1(t)=h2(t) for all t` versus `H1: h1(t) ≠ h2(t) for some t` (or the equivalent statement on `S(t)`, per Lec 22 slide 6). Do not accept silent omission.

3. **Test mechanics — at least one of the following is REQUIRED:**

   a) *By-hand path.* Per Lec 22 slides 7–9, the student must:
      - List the **uncensored event times** in the pooled sample (only event times generate 2×2 tables — censored times do NOT).
      - At each event time, give the at-risk count `n_i`, the per-group at-risk counts `n_{A,i}, n_{B,i}`, the per-group death counts `d_{A,i}, d_{B,i}`, and the pooled deaths `m_{D,i}`.
      - Compute expected deaths `E(d_{A,i}) = m_{D,i} n_{A,i}/n_i` and hypergeometric variance `Var(d_{A,i}) = n_{A,i} n_{B,i} m_{D,i} m_{\bar D,i} / [n_i^2(n_i-1)]` (Lec 22 slide 8). The reviewer should verify the variance formula uses `n_i - 1` in the denominator and check correct handling of tables where the variance is 0 (when m_D = n_i or m_{\bar D} = 0).
      - Form `Z = Σ(d_{A,i} - E(d_{A,i})) / sqrt(Σ Var(d_{A,i}))`, and report `Z^2 ~ χ^2_1` (Lec 22 slide 9).
      - The reviewer should independently rebuild the at-risk table from the raw data and verify each cell in the student's table line by line.

   b) *R path.* Use `survdiff(Surv(time, status) ~ group)` exactly as in Lec 22 slide 13. The reviewer should re-run the same code on the supplied data and confirm the printed `Chisq`, `df`, and `p` match. Verify the student does NOT pass a continuous predictor or otherwise mis-specify the call.

4. **Reported quantities.** The submission must explicitly state:
   - Test statistic value (chi-square or |Z|).
   - Degrees of freedom (must be 1 for two-group log-rank).
   - p-value.
   - A decision at α = 0.05 with one sentence of substantive interpretation (which group has worse / better survival, or "insufficient evidence").

**Common mistakes to watch for:**
- Treating the `+` symbol as part of the time value or as a separate group label.
- Including censored times as rows in the log-rank 2×2 tables (only event times produce tables).
- Using the binomial variance `np(1-p)` instead of the hypergeometric variance.
- Forgetting that `n_i - 1` (not `n_i`) appears in the variance denominator.
- Reporting df ≠ 1 for a two-group test.
- Drawing a "significant" conclusion from a tiny p-value without checking sample size — if the chi-square is small, the student must state `Do not reject H0` cleanly.
- Reversing the role of `time` vs `status` in `Surv()` — a frequent typo that invalidates the analysis.

---

## Problem 2 — KM curves / log-rank by race within each sex (kidtran data)

**Question type:** Mixed: visual + inferential. Must (i) produce KM survival curves stratified appropriately and (ii) test for race differences within each sex.

**Expected derivation type:** R code, with output shown. The student must use the `KMsurv` package and the survival-analysis toolchain demonstrated in Lec 21 slides 28–30 and Lec 22 slide 13.

**Interpretation of the question:** The instruction "Compare the survival curves for different races in each sex group" implies a **sex-specific (separate) analysis**, not a single pooled stratified test. Two acceptable structures:
   - Run **two separate analyses**, one on the male subset and one on the female subset, each comparing race 1 vs race 2 (this is the most direct reading and aligns with how Prof. Cheng presents two-sample log-rank in Lec 22).
   - Alternatively, a stratified-by-sex log-rank (using `strata(sex)` syntax in `survdiff`) is also acceptable IF the student also produces separate KM curves per sex so the visual comparison is clear. A single pooled four-group test (race × sex with 4 levels) does NOT answer the question and should be marked down.

**Key steps that MUST be shown:**

1. **Data loading & inspection.** Load `KMsurv`, `data(kidtran)`, identify the time variable, the event indicator, race, and gender columns. The reviewer should run `?kidtran` to verify the student names the columns and codes them consistently with the actual dataset.

2. **Subsetting / stratification.** The student must clearly produce two analyses (males-only and females-only). The subset criterion (`gender == 1` for males, `gender == 2` for females, per the homework codebook) must match the homework's stated coding.

3. **KM curves.** Per Lec 21 slides 28–30, fit `survfit(Surv(time, delta) ~ race, data = <sex-subset>)` and plot. A correct submission will show:
   - Two KM plots (one per sex) OR a single panel showing both, with race 1 vs race 2 visually distinguished.
   - Censoring marks visible (`mark.time = TRUE` in base `plot.survfit` or default in `ggsurvplot`).
   - Axes labeled (time, survival probability), legend identifying race.

4. **Log-rank test.** Per Lec 22 slide 13, run `survdiff(Surv(time, delta) ~ race)` separately in each sex subset (or the equivalent stratified call). Report:
   - Chi-square statistic, df = 1, p-value, **for each sex**.
   - Reject / do-not-reject decision at α = 0.05 for each sex.

5. **Substantive interpretation.** The student must write at least one sentence per sex group describing what the test (and KM plot) implies about race differences in that subgroup. Pure code dump without prose interpretation should lose points.

**R code reference:** the student's code should look stylistically like Lec 21 slide 28 (KM fit + plot) and Lec 22 slide 13 (`survdiff`). Reviewer compares argument style line-by-line to those slides.

**Common mistakes to watch for:**
- Misidentifying which column in `kidtran` is the time vs. the event indicator (verify against `?kidtran`).
- Using the wrong code for race or sex (e.g., assuming 0/1 coding when the data are 1/2).
- Pooling across sex when the question explicitly asks "in each sex group".
- Performing a single 4-group log-rank instead of two race comparisons stratified or separated by sex.
- Plotting only one curve or mislabeling the legend.
- Forgetting to convert numeric race/sex codes to factors when needed for plot legends, leading to ambiguous outputs.
- Reporting only the plot without the formal log-rank p-value, or vice versa.

---

## Problem 3 — Cox PH on larynx data with Z1·Z4 interaction; relative risk question

**Question type:** Two parts. (a) Fit and interpret a Cox PH model. (b) Direct, numeric answer: the relative risk of stage II vs stage I at age 50.

**Expected derivation type:** R code (`coxph`) + analytical computation of the relative risk from estimated coefficients.

**Key steps that MUST be shown:**

1. **Data setup.** Load `KMsurv`, `data(larynx)`. Construct the dummies Z1, Z2, Z3 from the `stage` variable (1=ref, 2→Z1, 3→Z2, 4→Z3) and use `age` as Z4. Either explicit dummy creation or `factor(stage)` with stage 1 as the reference level is acceptable, but the student MUST be clear that **stage 1 is the reference category** — otherwise the relative-risk question cannot be answered correctly.

2. **Model fit.** Per Lec 23 slides 26–27 (the larynx example), the student fits a Cox model of the form
   `coxph(Surv(time, delta) ~ Z1 + Z2 + Z3 + Z4 + Z1:Z4, data = larynx, ties = 'breslow')`
   or the equivalent using `factor(stage)` plus an explicit `Z1*age` interaction. The reviewer should verify:
   - Correct outcome construction `Surv(time, delta)` matching the larynx column names (Lec 23 slide 26).
   - All five required terms in the linear predictor: Z1, Z2, Z3, Z4, and the Z1:Z4 interaction (do NOT also include Z2:Z4 or Z3:Z4 unless the student justifies it; the homework specifies only Z1·Z4).
   - Use of `ties='breslow'` or `'efron'` declared (Efron is the default and is fine; the lecture example uses Breslow). Either choice is acceptable as long as it is reported.

3. **Output to display.** Student must print the coefficient table from `summary(fit)`, including:
   - Coefficient estimates (log hazard ratios).
   - Standard errors.
   - Wald z and p-values.
   - `exp(coef)` (hazard ratios).
   - 95% CIs for `exp(coef)`.
   - The omnibus tests (likelihood-ratio / Wald / score) reported by `summary` are useful but not required.

4. **Interpretation paragraph.** Each of Z1, Z2, Z3, Z4, and Z1:Z4 must be interpreted, with two CRITICAL conceptual points:
   - **For Z2 and Z3 (no interaction with age):** `exp(β_{Z2})` and `exp(β_{Z3})` are the hazard ratios of stage III vs I and stage IV vs I, respectively, **at any age** (because no interaction with age is fitted for those stages).
   - **For Z1 (which has an age interaction):** `exp(β_{Z1})` alone is NOT the hazard ratio of stage II vs I in general — it is the HR at **age = 0**. Because Z1·Z4 is in the model, the stage II vs stage I log hazard ratio at age `a` is `β_{Z1} + a·β_{Z1:Z4}`. This is the exact pitfall warned about in lecture; the student MUST explicitly state this. A submission that interprets `exp(β_{Z1})` as "the HR for stage II vs I" without qualification should lose points.
   - **For Z4 (age main effect):** because Z1·Z4 is in the model, `β_{Z4}` is the per-year log hazard ratio **for stage I (and stages III, IV) patients only**, not the average effect of age. The stage-II per-year log HR is `β_{Z4} + β_{Z1:Z4}`. The student should note this asymmetry.
   - **For Z1:Z4:** sign and magnitude interpretation — does the stage-II/stage-I gap widen or narrow with increasing age?
   - The student should state that hazard ratios are constant in time (PH assumption, Lec 23 slide 6) but here vary across covariate combinations.

5. **Relative risk: stage II (age 50) vs stage I (age 50).** This is a specific numeric question. The required derivation:
   - Both subjects share Z2=Z3=0 and Z4=50.
   - Stage I covariate vector: (Z1=0, Z2=0, Z3=0, Z4=50, Z1·Z4=0).
   - Stage II covariate vector: (Z1=1, Z2=0, Z3=0, Z4=50, Z1·Z4=50).
   - Log HR = `β_{Z1}·(1-0) + β_{Z1:Z4}·(50-0) = β_{Z1} + 50·β_{Z1:Z4}`.
   - Relative risk = `exp(β_{Z1} + 50·β_{Z1:Z4})`.
   - The student MUST plug in their estimated coefficients and report a single numeric RR. A delta-method or contrast-based 95% CI (e.g., via `multcomp::glht` or a hand-computed `Var(β_{Z1}) + 2500·Var(β_{Z1:Z4}) + 2·50·Cov(β_{Z1},β_{Z1:Z4})`) is NOT required by the prompt but is a nice extra; do not penalize for absence.
   - The reviewer should independently re-fit the Cox model on `larynx` and verify the student's numeric RR matches their own coefficients to ~3 sig figs.

**Common mistakes to watch for:**
- Failing to set stage 1 as the reference level (e.g., using `factor(stage)` and letting R default-order it differently because of how the variable was coded), which will silently make all stage HRs wrong.
- Reporting `exp(β_{Z1})` directly as the answer to the RR question — i.e., ignoring the interaction term entirely. This is the most common failure mode for this problem.
- Reporting `exp(β_{Z1} + β_{Z1:Z4})` — i.e., using `Z1·Z4 = 1` instead of `50`.
- Using `1·Z1` somewhere instead of plugging the actual age 50 into the interaction.
- Including Z2:Z4 or Z3:Z4 in the model (the homework specifies only Z1·Z4).
- Mis-coding the event indicator (`delta` in larynx is 0=alive, 1=dead per Lec 23 slide 26 — verify the student does not invert it).
- Treating stage as continuous rather than as dummy variables (only valid if treated as `factor(stage)` with reference level 1, which mechanically equals the dummy formulation).
- Interpreting Z4 as "the age effect" without noting that the age effect differs between stage I/III/IV vs stage II under this model.
- Not reporting whether ties were handled (Breslow vs Efron) — minor point but lecture notes raise it.
- Reporting only the coefficient table without prose interpretation.

---

## General submission requirements (apply to all problems)

- The student must include both code and prose. A `.Rmd` knit to PDF/HTML is the expected format.
- Numerical results must be reproducible — i.e., the reviewer running the included code on the named datasets should obtain the reported values.
- Hypotheses, test statistics, df, p-values, and α-level decisions must be reported in plain language for every formal test (Problems 1 and 2).
- For Problem 3, the relative-risk number must be a single, clearly stated value, not buried inside output.
- All R code should follow the `survival` + `KMsurv` toolchain demonstrated in Lectures 21–23 (specifically: `Surv`, `survfit`, `survdiff`, `coxph`). Use of alternative packages (e.g., `flexsurv`, `rms::cph`) is acceptable only if the student demonstrates equivalent output and explicitly justifies the choice.

## Reviewer checklist

- [ ] Re-derive the Problem 1 log-rank test from raw data; compare to student's test statistic, df, p-value, and decision.
- [ ] Re-run the Problem 2 KM and log-rank for each sex on `kidtran`; verify both p-values and KM plot shapes.
- [ ] Re-fit the Problem 3 Cox model on `larynx`; verify each coefficient, the interpretation of the Z1:Z4 term, and the numeric `exp(β_{Z1} + 50·β_{Z1:Z4})` matches.
- [ ] Confirm Problem 3 uses stage 1 as reference.
- [ ] Confirm censoring is handled correctly throughout (event indicator orientation, `+` notation interpretation in Problem 1).
