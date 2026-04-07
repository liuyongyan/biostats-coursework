# P8131 Homework 7 -- Grading Rubric

## General Information

- **Topic:** Linear Mixed Effects Models (random intercept models, BLUPs, likelihood ratio tests, crossed random effects)
- **Primary lecture references:** Lecture 17 (LME theory, slides 1-22), Lecture 18 (LME code examples, slides 1-11), Lecture 16 (marginal models / covariance structures, for background)
- **Data:** `HW7_politeness_data.csv` -- 6 subjects (F1, F2, F3, M3, M4, M7), 7 scenarios each, 2 attitude conditions (pol/inf), response = frequency (pitch). Each subject has 14 observations (7 scenarios x 2 attitudes).
- **R packages expected:** `nlme` (for `lme`) and/or `lme4` (for `lmer`). Part (d) explicitly requires `lmer` from `lme4`.

---

## Part (a): Exploratory Analysis -- Boxplots

**What is asked:** Provide boxplots showing the relationship between gender/attitude and pitch, ignoring different scenarios.

### Checking Criteria

1. **Boxplots present (3 pts)**
   - The student must produce at least one boxplot (or a panel of boxplots) that displays the distribution of pitch (`frequency`) broken down by gender and/or attitude.
   - Acceptable approaches: (i) a single boxplot with four groups on the x-axis (e.g., Female-pol, Female-inf, Male-pol, Male-inf), (ii) two side-by-side boxplots (one by gender, one by attitude), or (iii) a faceted boxplot using `ggplot2` with one factor on the x-axis and the other as fill/facet.
   - The student must use `frequency` as the y-axis variable.

2. **Correct grouping variables (2 pts)**
   - Both `gender` and `attitude` must appear as grouping variables in the plot(s). A plot that only shows one of the two factors is incomplete.
   - The problem says "ignoring different scenarios," so scenario should NOT be a grouping variable in the boxplot. If the student facets by scenario, that contradicts the instructions.

3. **Brief interpretation (1 pt)**
   - The student should note the visible patterns: the TA reviewer should check whether the student's verbal description is consistent with what the boxplots actually show (e.g., which group tends to have higher pitch, whether attitude appears to affect pitch, etc.).

### Common Pitfalls
- Using scenario as a grouping variable despite the instruction to ignore it.
- Plotting raw data points without a boxplot (a scatterplot alone does not satisfy the requirement).
- Mislabeling axes or factors.

---

## Part (b): Mixed Effects Model with Random Intercept for Subject

**What is asked:** Fit a mixed effects model with random intercepts for subjects, with gender and attitude as fixed effects. Then answer four sub-questions: (i) covariance matrix for Y_i, (ii) covariance matrix for fixed effect estimates (3x3), (iii) BLUPs for subject-specific intercepts, (iv) residuals.

### Checking Criteria

#### (b.1) Model fitting (3 pts)

1. The student must fit a linear mixed effects model using `lme()` from the `nlme` package (or `lmer()` from `lme4`). The lecture sample code in Lecture 17 slides 19-20 and Lecture 18 slides 6-7 uses `lme()` with the `nlme` package -- the TA should compare the student's function call style against those examples.
2. The model must have:
   - **Response:** `frequency` (pitch)
   - **Fixed effects:** `gender` and `attitude` (no interaction, no scenario)
   - **Random effect:** random intercept grouped by `subject`
3. The student should use REML estimation (the default for `lme`). Check that the student did not switch to `method='ML'` for this step (ML is needed later in part (c) for the LRT, but the model fit here should use REML for variance estimation).
4. The student should display the model summary.

#### (b.2) Covariance matrix for Y_i (4 pts)

1. The student must derive or state the marginal covariance matrix for the response vector of a single subject, Cov(Y_i).
2. The derivation must follow the random intercept model theory from Lecture 17 slide 6: for a random intercept model, the covariance structure is compound symmetry induced by the random effect.
   - Diagonal entries: sum of the between-subject variance (random intercept variance) and within-subject variance (residual variance).
   - Off-diagonal entries: the between-subject variance (random intercept variance) only.
3. The student must state the general form of this matrix AND plug in the estimated variance components from their fitted model to produce a numerical covariance matrix.
4. The dimension of the matrix must match the number of observations per subject. The TA should verify this against the data (each subject has 14 observations: 7 scenarios x 2 attitudes), so the covariance matrix should be 14 x 14. The student may write the general pattern and state the dimension rather than writing out the full 14x14 matrix, which is acceptable.

#### (b.3) Covariance matrix for fixed effects estimates (3 pts)

1. The problem provides a hint: this is a 3x3 matrix (intercept, gender, attitude).
2. The student must extract this from the fitted model. In `nlme`, this is obtained via `vcov()` applied to the model object. In `lme4`, it is similarly `vcov()`. The TA should compare against the extraction methods shown in Lecture 17 slide 20 and Lecture 18 slide 6 (where `vcov(LMM1)` is used).
3. The student must present the actual 3x3 numeric matrix. Verify that it is labeled or described as the covariance matrix of the fixed effect estimators (not the covariance of Y_i, which is a different quantity).
4. The TA should verify that the matrix is symmetric and positive definite (diagonal entries positive, and it "looks reasonable" relative to the standard errors in the model summary -- the diagonal entries should be the squares of the standard errors for each fixed effect).

#### (b.4) BLUPs for subject-specific intercepts (3 pts)

1. The student must extract the Best Linear Unbiased Predictors (BLUPs) for the random intercepts. In `nlme`, this is done via `random.effects()` or `ranef()`. In `lme4`, it is `ranef()`. See Lecture 17 slide 20 and Lecture 18 slide 6 for the sample code.
2. There should be one BLUP value per subject (6 subjects total: F1, F2, F3, M3, M4, M7).
3. The student should understand what BLUPs represent: the predicted deviation of each subject's intercept from the population mean intercept. The BLUP formula is given in Lecture 17 slides 13 and 16. The student does not need to re-derive the formula but should extract the values correctly.
4. The BLUPs should sum to approximately zero (a property of BLUPs in balanced/near-balanced designs). The TA can use this as a quick sanity check.

#### (b.5) Residuals (2 pts)

1. The student must extract the residuals from the model. In `nlme`, this is `residuals()` or accessing `LMM1$residuals`. See Lecture 17 slide 20 where `LMM1$residuals` is used.
2. The residuals should be the within-subject residuals (i.e., observed - fitted, where fitted includes both fixed effects and the random intercept). There should be one residual per observation (84 total, though some subjects may have slightly different counts if data is unbalanced -- the TA should check against the actual data).
3. Presenting a summary or plot of residuals is acceptable; printing all values is also acceptable.

### Common Pitfalls
- Confusing Cov(Y_i) (the marginal covariance of the response) with Cov(beta-hat) (the covariance of the fixed effect estimates). These are completely different quantities and the problem asks for both.
- Forgetting to include the residual variance on the diagonal of Cov(Y_i).
- Extracting fitted values instead of BLUPs (fitted values = X*beta-hat + Z*b-hat, whereas BLUPs are just the b-hat).
- Using `coef()` instead of `ranef()` -- `coef()` gives the subject-specific intercept (population intercept + BLUP), not the BLUP alone.

---

## Part (c): Model with Interaction Term and Likelihood Ratio Test

**What is asked:** Fit a model that adds the gender-by-attitude interaction to the model in part (b). Use a likelihood ratio test to compare the two models.

### Checking Criteria

#### (c.1) Model fitting (3 pts)

1. The student must fit a new mixed effects model identical to part (b) but with the addition of a `gender:attitude` interaction term in the fixed effects.
2. The random effect structure must remain the same (random intercept for subject).
3. **Critical:** Both models (from part (b) and this part) must be re-fit using **maximum likelihood (ML)**, NOT REML, for the likelihood ratio test. This is because the two models differ in their fixed effects, and REML likelihoods are not comparable when fixed effects differ. This principle is stated clearly in Lecture 17 slide 20 ("do NOT use REML for likelihood ratio") and demonstrated in Lecture 18 slide 10 where `method='ML'` is used for the LRT.
4. If the student uses `anova()` to compare the two models, the TA should verify that both models were fit with ML. In `lme4`, if using `anova()` on `lmer` objects fit with REML, some implementations automatically refit with ML, but the student should be explicit about this.

#### (c.2) Likelihood ratio test (3 pts)

1. The student must perform the LRT using `anova()` to compare the two nested models. See Lecture 17 slide 20 and Lecture 18 slide 10 for the pattern.
2. The test has one degree of freedom (the interaction term adds one parameter to the fixed effects).
3. The student must state the null and alternative hypotheses: H0 is that the interaction coefficient is zero; H1 is that it is not zero.
4. The student must report the test statistic, p-value, and draw a conclusion about whether the interaction is significant at a standard significance level (typically alpha = 0.05).

#### (c.3) Conclusion (2 pts)

1. The student must state a clear conclusion: whether or not the interaction between gender and attitude is significantly associated with pitch, based on the p-value from the LRT.
2. The conclusion must be consistent with the reported p-value.

### Common Pitfalls
- Using REML for the LRT when the fixed effects differ between models. This is the single most important thing to check in this part.
- Comparing models that differ in random effects structure (they should not -- both models should have the same random intercept for subject).
- Forgetting to state the conclusion in words.

---

## Part (d): Mixed Effects Model with Crossed Random Intercepts (Subject and Scenario)

**What is asked:** (i) Write out the mixed effects model with random intercepts for both subjects AND scenarios. (ii) Fit using `lmer` from `lme4`. (iii) Write out Cov(Y_i). (iv) Interpret the coefficient for attitude.

### Checking Criteria

#### (d.1) Model specification in mathematical notation (4 pts)

1. The student must write the model equation explicitly. The model should have:
   - Fixed effects: intercept, gender, attitude
   - Random intercept for subject: b_i ~ N(0, sigma_b^2)
   - Random intercept for scenario: c_j ~ N(0, sigma_c^2) (or similar notation)
   - Residual error: epsilon_ij ~ N(0, sigma^2)
2. The model equation should clearly show both random intercepts as additive terms.
3. The student must state the distributional assumptions: the two random effects and the residual error are mutually independent, each normally distributed with mean zero.
4. The notation should follow the general LME framework from Lecture 17 slides 11-12 (Y_i = X_i * beta + Z_i * b_i + epsilon_i), adapted to include crossed random effects.

#### (d.2) Model fitting with lmer (3 pts)

1. The problem explicitly requires the use of `lmer()` from the `lme4` package (not `lme()` from `nlme`). This is because `lme()` in `nlme` handles nested random effects naturally but crossed random effects require special syntax, whereas `lmer()` handles crossed random effects straightforwardly. The TA should verify the student uses `lmer`.
2. The model formula must include two separate random intercept terms: one for `subject` and one for `scenario`. In `lmer` syntax, this means the formula should have two `(1|...)` terms.
3. The student should display the model summary showing the variance components for both random effects and the residual, plus the fixed effects estimates.

#### (d.3) Covariance matrix for Y_i (4 pts)

1. The student must derive or write out the covariance matrix for a subject's response vector Y_i.
2. With two crossed random effects (subject and scenario), the covariance structure is more complex than part (b). The TA should verify that:
   - The student correctly identifies the dimension of Y_i for a given subject (14 observations: 7 scenarios x 2 attitudes per scenario).
   - The diagonal entries include the subject variance, the scenario variance, and the residual variance.
   - Off-diagonal entries for observations from the **same scenario** (but different attitudes) include both the subject variance and the scenario variance.
   - Off-diagonal entries for observations from **different scenarios** (within the same subject) include only the subject variance.
3. The student must show the general structure of this matrix and explain why different off-diagonal elements differ. The matrix should exhibit a block-like structure reflecting the crossed nature of the random effects.
4. The student should plug in the estimated variance components from the `lmer` output to produce numerical values.

#### (d.4) Interpretation of the attitude coefficient (3 pts)

1. The student must provide a clear, contextually meaningful interpretation of the fixed effect coefficient for attitude.
2. The interpretation must reflect what this coefficient means in a mixed effects model context. Per Lecture 17 slide 12, the fixed effect coefficients in an LME model have a **population-average (marginal) interpretation**: E(Y_i) = X_i * beta.
3. The interpretation should specify:
   - The direction of the effect (which level of attitude is the reference category -- likely "inf" or "pol" depending on alphabetical ordering or factor level setting).
   - The magnitude: the coefficient represents the average difference in pitch between the two attitude levels, holding gender constant, across all subjects and scenarios.
   - The student should note that this is a population-level (marginal) interpretation, not a subject-specific one.
4. The interpretation must be stated in the context of the problem (pitch/frequency, politeness).

### Common Pitfalls
- Using `lme()` instead of `lmer()` when the problem explicitly asks for `lmer`.
- Writing only one random intercept term (for subject) and forgetting the scenario random intercept.
- In the covariance matrix, treating all off-diagonal entries as equal (that would only be correct if there were a single random intercept). With crossed random effects, off-diagonal entries differ depending on whether two observations share the same scenario.
- Giving a subject-specific (conditional) interpretation of the attitude coefficient instead of a population-average (marginal) interpretation, or failing to be clear about which interpretation is intended.
- Forgetting that in an LME model, fixed effect coefficients have the same marginal interpretation as in a standard linear model (Lecture 17 slide 12).

---

## Overall Grading Summary

| Part | Component | Points |
|------|-----------|--------|
| (a)  | Boxplots present and correctly structured | 3 |
| (a)  | Both gender and attitude used; scenarios ignored | 2 |
| (a)  | Brief interpretation | 1 |
| (b)  | Model correctly specified and fit | 3 |
| (b)  | Covariance matrix for Y_i (structure + numerics) | 4 |
| (b)  | Covariance matrix for fixed effects (3x3) | 3 |
| (b)  | BLUPs extracted and presented | 3 |
| (b)  | Residuals extracted | 2 |
| (c)  | Interaction model fit correctly | 3 |
| (c)  | LRT performed with ML (not REML) | 3 |
| (c)  | Clear conclusion stated | 2 |
| (d)  | Model written in mathematical notation | 4 |
| (d)  | Model fit with lmer, correct formula | 3 |
| (d)  | Covariance matrix for Y_i with crossed random effects | 4 |
| (d)  | Interpretation of attitude coefficient | 3 |
| **Total** | | **43** |

*Note: The point allocation above is a suggested framework. The instructor may adjust weights as appropriate. The TA should scale to whatever total the course uses.*

---

## Key Lecture References for the TA Reviewer

- **Random intercept model theory and Cov(Y_i) derivation:** Lecture 17, slides 3-6
- **General LME framework (Y = X*beta + Z*b + epsilon):** Lecture 17, slides 11-12
- **BLUP formula:** Lecture 17, slides 13, 16
- **Covariance of fixed effect estimates (vcov):** Lecture 17, slide 17 (theory), slide 20 (code)
- **R code for lme(): fixed.effects, random.effects, residuals, vcov:** Lecture 17, slides 19-20
- **LRT requirement to use ML not REML:** Lecture 17, slide 20; Lecture 18, slide 10
- **R code for lme() with real data (TURBT/GCase examples):** Lecture 18, slides 6-10
- **lmer for crossed random effects:** The problem directs students to use `lmer` from `lme4`; Lecture 17 slide 18 notes that `lmer` is more suitable for multiple non-nested random effects
