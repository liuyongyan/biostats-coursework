# P8131 Homework 6 -- Grading Rubric

Total: 40 points

---

## Problem 1: Variance, Covariance, Correlation, and Covariance Pattern Identification (16 points)

**Model:** Y_ij = mu + b_i + e_ij, where b_i ~ N(0, sigma_b^2), e_ij ~ N(0, sigma_e^2), all mutually independent.

### 1a. Variance of Y_ij (4 points)

- **Correct answer:** Var(Y_ij) = Var(b_i) + Var(e_ij) = sigma_b^2 + sigma_e^2
- **Key steps that MUST be shown:**
  - Recognize that mu is a constant and contributes zero variance (1 pt)
  - Use the independence of b_i and e_ij to split Var(b_i + e_ij) = Var(b_i) + Var(e_ij) (1 pt)
  - Substitute the given distributions to get sigma_b^2 + sigma_e^2 (2 pts)
- **R code:** Not required (pure derivation)
- **Common mistakes:**
  - Forgetting to justify the split of variance by citing independence
  - Including mu in the variance calculation

### 1b. Covariance between Y_ij and Y_ik, j != k (4 points)

- **Correct answer:** Cov(Y_ij, Y_ik) = Cov(b_i + e_ij, b_i + e_ik) = Var(b_i) = sigma_b^2
- **Key steps that MUST be shown:**
  - Expand using bilinearity of covariance (1 pt)
  - Use independence of e_ij and e_ik (given j != k) so Cov(e_ij, e_ik) = 0 (1 pt)
  - Use independence of b_i and e_ij so cross-terms Cov(b_i, e_ij) = 0 and Cov(b_i, e_ik) = 0 (1 pt)
  - Conclude Cov(Y_ij, Y_ik) = Var(b_i) = sigma_b^2 (1 pt)
- **R code:** Not required
- **Common mistakes:**
  - Not showing which cross-covariance terms vanish and why

### 1c. Correlation between Y_ij and Y_ik, j != k (4 points)

- **Correct answer:** Corr(Y_ij, Y_ik) = sigma_b^2 / (sigma_b^2 + sigma_e^2)
- **Key steps that MUST be shown:**
  - Write Corr = Cov(Y_ij, Y_ik) / sqrt(Var(Y_ij) * Var(Y_ik)) (1 pt)
  - Note Var(Y_ij) = Var(Y_ik) = sigma_b^2 + sigma_e^2 (1 pt)
  - Simplify to sigma_b^2 / (sigma_b^2 + sigma_e^2) (2 pts)
- **R code:** Not required
- **Common mistakes:**
  - Algebraic errors in simplification
  - Not recognizing that variances are equal for all j

### 1d. Identify the covariance pattern (4 points)

- **Correct answer:** This is the **compound symmetry** (also called exchangeable) covariance pattern.
- **This question requires a direct, explicit answer** -- the student must NAME the pattern.
- **Key steps that MUST be shown:**
  - State that the variance is constant across all observations (sigma_b^2 + sigma_e^2) (1 pt)
  - State that the covariance (sigma_b^2) is constant for all pairs j != k (1 pt)
  - Equivalently, the correlation rho = sigma_b^2 / (sigma_b^2 + sigma_e^2) is constant for all pairs (1 pt)
  - Explicitly name the pattern as compound symmetry / exchangeable (1 pt)
- **Reference:** Lecture 16, slide 8 -- compound symmetry has constant variance sigma^2 and constant correlation rho for all pairs. The covariance matrix has the form Cov(Y_i) = sigma^2 * (matrix with 1 on diagonal and rho on off-diagonal). Lecture 17, slide 6 shows the random intercept model induces exactly this structure.
- **R code:** Not required
- **Common mistakes:**
  - Deriving the variance and covariance correctly but forgetting to explicitly name the pattern
  - Confusing compound symmetry with AR(1) or Toeplitz

---

## Problem 2a: Spaghetti Plot (4 points)

- **Correct approach:** Create a spaghetti plot (line plot) with age on the x-axis, distance on the y-axis, one line per child, colored/grouped by gender.
- **Key elements that MUST be present:**
  - Lines connecting the 4 time points for each child (1 pt)
  - Visual distinction between boys and girls (color, linetype, or facet) (1 pt)
  - Proper axis labels and a legend (1 pt)
  - Brief verbal description of what the plot shows -- e.g., distance generally increases with age; boys tend to have larger distances or steeper slopes than girls (1 pt)
- **R code:** Required. Acceptable approaches include ggplot2 with `geom_line(aes(group = Child, color = factor(Gender)))` or base R equivalent.
- **Common mistakes:**
  - Plotting points without connecting lines (scatter plot, not spaghetti)
  - Not differentiating by gender
  - No interpretation of the plot at all

---

## Problem 2b: Write the Mixed Model in Marginal Form (10 points)

**Given model:** Y_ij = beta_0 + a_i + b_0 * I(sex_i=0) + b_1 * I(sex_i=1) + beta_1 * age_ij + e_ij

where a_i ~ N(0, sigma_a^2), b_k ~ N(0, sigma_b^2) for k=0,1, e_ij ~ N(0, sigma_e^2), all mutually independent.

### 2b-i. Marginal mean E(Y_ij) (3 points)

- **Correct answer:** E(Y_ij) = beta_0 + beta_1 * age_ij
- **Key steps that MUST be shown:**
  - Take expectation of both sides (1 pt)
  - E(a_i) = 0, E(b_0) = 0, E(b_1) = 0 since they are zero-mean normal random variables (1 pt)
  - Simplify to beta_0 + beta_1 * age_ij (1 pt)
- **Common mistakes:**
  - Including the gender indicator terms in the marginal mean. NOTE: In this model, b_0 and b_1 are RANDOM coefficients (not fixed), so their expected values are zero. The gender effect enters only through the variance structure, not the mean. This is a subtle but critical point.
  - Confusing the model in 2(b) (where sex effects are random) with the model in 2(c) (where sex is a fixed effect in the mean)

### 2b-ii. Marginal variance Var(Y_ij) (3 points)

- **Correct answer:**
  - For a girl (sex_i = 0): Var(Y_ij) = sigma_a^2 + sigma_b^2 + sigma_e^2
  - For a boy (sex_i = 1): Var(Y_ij) = sigma_a^2 + sigma_b^2 + sigma_e^2
  - (Same formula for both, since b_0 and b_1 have the same variance sigma_b^2)
- **Key steps that MUST be shown:**
  - Identify which random terms contribute to variance for a given subject (1 pt)
  - Use mutual independence of a_i, b_k, e_ij to sum their variances (1 pt)
  - Arrive at sigma_a^2 + sigma_b^2 + sigma_e^2 (1 pt)
- **Common mistakes:**
  - Double-counting by including both b_0 and b_1 for a single subject (each subject is either male or female, so only one of b_0, b_1 applies)

### 2b-iii. Marginal covariance Cov(Y_ij, Y_ik) for j != k (same subject i) (4 points)

- **Correct answer:**
  - Cov(Y_ij, Y_ik) = sigma_a^2 + sigma_b^2
- **Key steps that MUST be shown:**
  - Expand Cov(Y_ij, Y_ik) using the model (1 pt)
  - The shared random effects for subject i are a_i and one of b_0 or b_1 (depending on sex); e_ij and e_ik are independent (1 pt)
  - All cross-terms between independent random effects vanish (1 pt)
  - Result: Cov(Y_ij, Y_ik) = Var(a_i) + Var(b_{sex_i}) = sigma_a^2 + sigma_b^2 (1 pt)
- **Bonus/optional:** Noting that this is again a compound symmetry structure with rho = (sigma_a^2 + sigma_b^2) / (sigma_a^2 + sigma_b^2 + sigma_e^2).
- **R code:** Not required (pure derivation)
- **Common mistakes:**
  - Forgetting the b_k contribution to the covariance
  - Not clearly stating the marginal covariance matrix Cov(Y_i) in matrix form

---

## Problem 2c: Fit Marginal Models with Three Covariance Structures (10 points)

**Marginal model:** E(Y_ij) = beta_0 + beta_1 * sex_i + beta_2 * age_ij, Var(Y_i) = Sigma.

Fit with: (a) compound symmetry, (b) exponential, (c) autoregressive (AR(1)) covariance.

### 2c-i. Correct R implementation (4 points)

- **R code is REQUIRED.**
- **Correct approach:** Use `gls()` from the `nlme` package with appropriate `correlation` argument:
  - Compound symmetry: `correlation = corCompSymm(form = ~ 1 | Child)` (1 pt)
  - Exponential: `correlation = corExp(form = ~ 1 | Child)` (1.5 pts) — must match lecture sample code style
  - Autoregressive AR(1): `correlation = corAR1(form = ~ 1 | Child)` (1.5 pts)
- **Key details:**
  - The grouping variable in the formula must be `Child` (subject ID), NOT `Index`
  - Method should be `"REML"` (default for gls, but acceptable to state explicitly)
  - The mean model should be `Distance ~ Gender + Age`
- **Common mistakes:**
  - Using `lm()` or `glm()` instead of `gls()` -- these ignore the correlation structure
  - Incorrect grouping variable (e.g., using `Index` instead of `Child`)
  - Forgetting to convert Gender to a factor or treating it incorrectly
  - Using `corSymm` (unstructured) when compound symmetry is asked for
  - Using `lme()` instead of `gls()` -- `lme()` is for mixed effects models, not marginal models. While results may be similar, the question explicitly asks for a marginal model.

### 2c-ii. Report and compare coefficient estimates (3 points)

- **This question requires a direct comparison** -- the student must present estimates from all three models side by side and comment.
- **Key elements that MUST be present:**
  - Report beta_0 (intercept), beta_1 (sex/gender effect), beta_2 (age effect) for all three models (1.5 pts)
  - Compare: Fixed effect estimates should be similar across models since the mean structure is the same; differences arise due to different weighting from the covariance structure (1.5 pts)
- **Expected results (approximate):**
  - Intercept ~ 17; Gender (boy) effect ~ 2; Age effect ~ 0.66
  - These should be broadly similar across the three covariance structures
- **Common mistakes:**
  - Only reporting output without any comparison or commentary
  - Only reporting one model and not the others

### 2c-iii. Report and compare covariance estimates (3 points)

- **This question requires a direct comparison.**
- **Key elements that MUST be present:**
  - For compound symmetry: report sigma^2 and rho (1 pt)
  - For exponential: report sigma^2 and range parameter (1 pt)
  - For AR(1): report sigma^2 and phi (1 pt)
- **Expected discussion points:**
  - The implied correlation matrices differ: CS has constant off-diagonal correlation; AR(1) has correlation decaying as rho^|j-k|; exponential has correlation decaying as rho^|t_j - t_k| (generalizes AR to unequal spacing)
  - For equally spaced data (ages 8, 10, 12, 14 with spacing 2), AR(1) and exponential can produce very similar (or identical) results
  - Note which model might be most appropriate for this data
- **Common mistakes:**
  - Not extracting or displaying the covariance/correlation matrices
  - No substantive comparison, just dumping R output

---

## General Deductions

- **Notation inconsistency or sloppiness:** Up to -1 per problem
- **Missing interpretation when explicitly asked:** -1 to -2 depending on severity
- **Code runs but produces errors/warnings that are not addressed:** -1
- **Excessive filler or irrelevant content that obscures the answer:** -1

## Summary Table

| Component | Points |
|-----------|--------|
| 1: Var(Y_ij) | 4 |
| 1: Cov(Y_ij, Y_ik) | 4 |
| 1: Corr(Y_ij, Y_ik) | 4 |
| 1: Name the covariance pattern | 4 |
| 2a: Spaghetti plot | 4 |
| 2b: Marginal form derivation | 10 |
| 2c: Fit and compare three marginal models | 10 |
| **Total** | **40** |
