# P8131 Spring 2026 — Homework #9 Grading Rubric

This rubric specifies **principles** the TA reviewer should check the student's work against. It does **not** contain pre-computed answers; the TA reviewer is responsible for verifying correctness independently using the lecture notes.

Lecture references:
- **Lecture 21** (intro to survival analysis): definitions of S(t), F(t), f(t), h(t), H(t); the identity relating them; KM estimator; Nelson–Aalen estimator; Fleming–Harrington estimator; Greenwood/log-CI; R code with `survival::survfit`. The acute leukemia (Gehan / 6-MP) example on slides 17–24 and the R code block on slides 28–30 are the canonical templates.
- **Lecture 22** (comparing survival functions): general background only; HW9 does not require log-rank tests.

---

## Problem 1 — Determine S(x) and f(x) from h(x) = 2x / (1 + x²)

### Question type
Analytical derivation. The problem asks for two specific functional forms: the survival function S(x) and the density function f(x). Both should be presented as closed-form expressions of x, valid for x ≥ 0.

### Required form
- Pure analytical derivation (no R code required).
- Final answer should be two clean closed-form expressions.
- The derivation must use the identities introduced in Lecture 21 slides 13–14.

### Key steps that MUST be shown (checklist of principles)
- [ ] Must explicitly invoke the relationship **S(t) = exp(−H(t))** with **H(t) = ∫₀ᵗ h(u) du**, OR equivalently the relationship **h(t) = −d/dt log S(t)**. (Lecture 21 slide 14.)
- [ ] Must compute the cumulative hazard H(x) by integrating h(u) from 0 to x. The hint about the derivative of log(1 + x²) should be used (or an equivalent substitution u² = v) — the student should recognize that h is the derivative of a log-of-quadratic.
- [ ] Must show at least the integration step (substitution or antiderivative recognition); a one-line "by inspection" claim with no justification is insufficient.
- [ ] Must obtain S(x) by exponentiating −H(x).
- [ ] Must obtain f(x) using **f(x) = h(x) · S(x)** OR **f(x) = −dS/dx**. Whichever route is used must be stated.
- [ ] Must specify the support (x ≥ 0). Bonus if the student remarks that S(0) = 1 and S(∞) = 0 are satisfied (sanity check).

### Method/code consistency
- Compare against Lecture 21 slides 10, 13, 14 (definitions of S, h, H and the identities among f, h, S).

### Common mistakes / pitfalls (conceptual)
- Forgetting the negative sign in S(t) = exp(−H(t)).
- Treating h(x) itself as the density f(x) (these are equal only when S = 1, i.e., not in general).
- Producing an S(x) that does not satisfy S(0) = 1 or that is not monotonically decreasing on [0, ∞) — either is a red flag.
- Computing the integral by an incorrect antiderivative (the hint is meant to prevent this; a student who ignores the hint and uses an unrelated substitution should still arrive at the same closed form, but check the algebra carefully).
- Stopping at H(x) without producing S(x), or stopping at S(x) without producing f(x).

---

## Problem 2 — Hand calculation on data 1, 2, 2, 4+, 5+, 6, 7+, 8+, 9+, 10+

### Question type
Hand calculation. The problem asks the student to **write out the data table** and then derive three quantities from it:
- (a) Kaplan–Meier estimate of S(t)
- (b) Nelson–Aalen estimate of H(t)
- (c) Fleming–Harrington estimate of S(t)

Each part requires a numerical estimate at every distinct event time (not merely at one chosen t).

### Required form
- A single hand-built table is the natural format. R or calculator arithmetic is acceptable for individual fractions, but the **table layout itself must be shown by hand** (i.e., not delegated to `survfit`).
- All three estimators should appear together in the table for direct comparability, OR in three parallel tables. Either layout is acceptable as long as columns are clearly labeled.

### Key steps that MUST be shown (checklist of principles)

**Common to (a)–(c): the data table**
- [ ] Must list the **distinct event times** t_i (uncensored times only — censored times are not rows of the estimator tables, though they affect n_i).
- [ ] Columns required: **t_i, n_i (number at risk just before t_i), d_i (number of events at t_i), c_i (number censored in [t_i, t_{i+1}))**. (Same column convention as Lecture 21 slide 24.) A column for the conditional hazard estimate λ̂_i = d_i/n_i is also expected.
- [ ] Must correctly identify ties at t = 2 (two events at the same time) and incorporate them as d_i = 2.
- [ ] Must correctly handle the convention that a censored observation at time c is at-risk at all event times t ≤ c but not at any t > c. Reviewer should verify the n_i column is consistent with this convention by walking the data top to bottom.

**(a) Kaplan–Meier**
- [ ] Must apply the product formula **Ŝ(t) = ∏_{t_i ≤ t} (1 − d_i/n_i)** (Lecture 21 slide 23).
- [ ] Must produce Ŝ(t) values at each event time as a running product.
- [ ] Should make explicit (in words or by step-function notation) that Ŝ is constant between event times.

**(b) Nelson–Aalen**
- [ ] Must apply the formula **H̃(t) = ∑_{t_i ≤ t} d_i/n_i** (Lecture 21 slide 27).
- [ ] Must produce H̃(t) values at each event time as a running sum.
- [ ] H̃(t) = 0 for t < t_1 should be stated.

**(c) Fleming–Harrington**
- [ ] Must use **Ŝ_FH(t) = exp(−H̃(t))** with **H̃ taken from part (b)** (Lecture 21 slide 27 — exp(−H̃(t)) is the Fleming–Harrington estimator). The student must NOT instead use H_KM(t) = −log Ŝ_KM(t) here; the whole point of part (c) is the Nelson–Aalen-based version.
- [ ] Final values for each event time should appear in the table.

### Method/code consistency
- Compare formulas and table layout against Lecture 21 slides 22–24 (KM table) and slide 27 (Nelson–Aalen and Fleming–Harrington).

### Common mistakes / pitfalls (conceptual)
- **Wrong n_i**: confusing "number at risk just before t_i" with "number surviving after t_i". The student must use the *pre-event* count.
- **Mishandling ties at t = 2**: using d = 1 twice instead of d = 2 once, or listing t = 2 as two separate rows.
- **Mishandling censored observations** (the four "+" times 4, 5, 7, 8, 9, 10): a censored observation at c contributes to n_i for all event times t_i ≤ c. A common error is removing censored subjects from the risk set too early or too late.
- **(c) using KM-derived H instead of Nelson–Aalen H**: the lecture explicitly distinguishes Ĥ(t) = −log Ŝ_KM(t) from H̃(t) (Nelson–Aalen). Fleming–Harrington = exp(−H̃), not exp(+log Ŝ_KM). Watch for this confusion.
- Reporting only a final value at the last event time rather than the full step function.
- Off-by-one in the running product/sum.

---

## Problem 3 — `tongue` data from `KMsurv`: KM curves with log-transform 95% CI, plus 1-year survival

### Question type
**Two parts, each producing direct, explicit answers:**
1. **Plots** (one per tumor type, or overlaid): KM curve of Ŝ(t) with **pointwise 95% CI using the log transformation**.
2. **Direct numerical answers**: the **estimated 1-year survival rate AND its 95% CI** for **each** tumor type (aneuploid and diploid). This is **four numbers per group** (point estimate, lower CI, upper CI for each of two groups) that must be explicitly stated, not merely buried in code output.

### Required form
- R code (must be visible / runnable) using the `KMsurv::tongue` dataset and the `survival` package.
- One or two plots showing both the KM step function and its 95% pointwise confidence band.
- An explicit written statement of the 1-year survival estimate and CI for each group (e.g., "For aneuploid tumors, the estimated 1-year survival is __ with 95% CI (__, __)." and similarly for diploid).

### Key steps that MUST be shown (checklist of principles)
- [ ] Must load `KMsurv` and access the `tongue` dataset; must understand the variable encoding (the dataset has a tumor `type` variable, a `time` variable, and an event/censoring indicator `delta`). The student must subset/stratify by tumor type correctly (one group = aneuploid, other group = diploid).
- [ ] Must construct a `Surv()` object with the correct event indicator.
- [ ] Must call `survfit()` with the **log transformation for the CI** — i.e., the argument that requests CIs on the log scale (this is the `conf.type = "log"` style shown in Lecture 21 slide 28). The reviewer should confirm the student passed the log-transform option explicitly rather than relying on default behavior, since the problem asks for it specifically.
- [ ] Must produce a plot showing the step function AND the confidence intervals (i.e., `conf.int = TRUE` or equivalent). A plot without CIs does not satisfy the question.
- [ ] Must extract / read off the 1-year (i.e., t = 365 days, since `tongue` records time in **weeks** — reviewer must check the unit convention; if `tongue` time is in weeks, "1 year" = 52 weeks; if in days, 365 days). The reviewer should verify the student is consistent: whatever unit the dataset uses, "1 year" must be expressed in that same unit.
- [ ] Must report 1-year survival and 95% CI separately for each group. Reading these from `summary(survfit_obj, times = <1 year>)` is the standard route (Lecture 21 slide 28).

### Method/code consistency
- Compare against Lecture 21 slides 28–30, where the canonical workflow is shown:
  - `library(survival)` + `Surv(time, cens, type = "right")`
  - `survfit(Surv(...) ~ 1, data = subset(...), conf.type = "log")`
  - `plot(KM, conf.int = ..., mark.time = TRUE, ...)`
  - `summary(KM, times = c(...))` to read off survival rate and CI at chosen time(s)
- The student's code does not have to match line-by-line, but the **functions used** (`Surv`, `survfit`, `plot`, `summary`) and the **conf.type argument** should match.

### Direct-answer items (must be explicit in the writeup)
- Aneuploid 1-year survival point estimate: ________
- Aneuploid 1-year survival 95% CI: (________, ________)
- Diploid 1-year survival point estimate: ________
- Diploid 1-year survival 95% CI: (________, ________)

The reviewer should look for these four numbers as written-out statements (not only as raw R console output buried in the document).

### Common mistakes / pitfalls (conceptual)
- **Wrong time unit**: assuming `tongue$time` is in days when it is in weeks (or vice versa). The reviewer should verify the student picked the correct numerical value of "1 year" for the dataset's units before judging the answer.
- **Wrong CI type**: using the default Greenwood (plain) CI instead of the log-transform CI. The default in `survfit` is not always `"log"`; the student must request it explicitly. (Note: there are several related options — `"log"`, `"log-log"`, `"plain"` — and the homework asks for the log transformation specifically; the student should not silently use `"log-log"`.)
- **Pooling tumor types**: fitting a single KM curve on the whole dataset rather than stratifying by tumor type.
- **Mis-subsetting**: confusing which level of the `type` variable corresponds to aneuploid vs. diploid (the dataset codes them numerically — the student should confirm the coding from the `KMsurv` documentation and label the output accordingly).
- **Plot without CIs**, or **CIs without the requested log transformation**.
- **Reporting only the point estimate** without the CI, or only at an event time near 1 year rather than at exactly 1 year (use `summary(..., times = ...)` which interpolates correctly to the requested time by carrying forward the last KM value).
- **Picking the wrong row** when reading off the survival at 1 year (e.g., reading `n.event` instead of `survival`, or reading the row immediately after 1 year instead of the row representing the value at 1 year).
- Failing to clearly **label which group is which** in the final reported numbers.

---

## General notes for the TA reviewer
- For Problem 1, independently verify the integral and the resulting S(x) and f(x) from scratch — do not rely on any cached "standard answer".
- For Problem 2, independently rebuild the table row by row, paying attention to ties at t = 2 and to the correct decrement of n_i at each censored time.
- For Problem 3, independently run the canonical Lecture 21 workflow on `KMsurv::tongue` to confirm the four reported numbers, and confirm the time-unit convention of the dataset before judging "1 year".
- Penalize missing derivation steps even when the final answer is correct: this homework is graded on the work, not just the result.
- Award partial credit for correct method with arithmetic errors (Problem 2) or correct setup with one wrong argument (Problem 3).
