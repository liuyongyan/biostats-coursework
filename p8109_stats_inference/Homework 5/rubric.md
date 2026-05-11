# P8109 Homework 5 — Grading Rubric (Principles-Based)

This rubric specifies what each problem **requires** of the student and how a TA reviewer should **verify** the work. It does **not** contain pre-computed answers. The reviewer must independently confirm every numerical/symbolic result against the principles below.

General review notes:
- For every problem, locate the corresponding chunk in the student's RMD/PDF by line number and confirm the equation chain reads "prior × likelihood → kernel of known distribution → posterior" before checking any numbers.
- If a numeric answer is requested but the derivation skeleton is missing or unjustified, deduct for missing derivation regardless of whether the final number happens to be right.
- Conjugacy claims (e.g., "Beta + Binomial → Beta") must be either cited from Lecture 10 (slide 19, "Conjugacy" remark; Example 10.4) or derived; they may not be asserted without one of those.

---

## Problem 1

**What the question asks.** A specific posterior probability is required: $P(p \le 0.1 \mid X = 2, n = 30)$ under a Uniform$(0,1)$ prior and Binomial$(30, p)$ likelihood. The student must produce a single numerical answer.

**Type of work required.** Both analytical (identify the posterior as a Beta) and R (compute the Beta CDF tail).

**Key steps that MUST be shown.**
1. Likelihood written as $L(p \mid x) \propto p^{x}(1-p)^{n-x}$ with $n=30, x=2$.
2. Prior written as Uniform$(0,1) = \mathrm{Beta}(1,1)$, i.e., $\pi(p) \propto 1$.
3. Posterior kernel obtained by "posterior $\propto$ prior $\times$ likelihood":
   $$\pi(p \mid x) \propto p^{x+1-1}(1-p)^{n-x+1-1}.$$
4. Recognition of the kernel as a Beta distribution and explicit statement of its two shape parameters as functions of $x, n$. Reviewer must verify the parameter arithmetic independently.
5. The required tail probability written as a Beta CDF, e.g., $P(p \le 0.1 \mid X) = F_{\mathrm{Beta}}(0.1; \alpha^*, \beta^*)$, evaluated via R's `pbeta(0.1, ...)`. The R code chunk must appear and its numeric output must match the analytical claim.

**Common mistakes to watch for.**
- Treating the Uniform prior as a separate object from Beta$(1,1)$ and then forgetting the "$+1$" shifts in shape parameters.
- Off-by-one in shape parameters (writing $\mathrm{Beta}(x, n-x)$ instead of the conjugate-update form).
- Computing $P(p > 0.1)$ or using `1 - pbeta(...)` when the question asks for $\le$.
- Reporting the prior probability $P(p \le 0.1) = 0.1$ (under Uniform) instead of the **posterior** probability.

**Lecture reference.** Lecture 10 slide 10 (Beta posterior for Binomial under uniform prior); Example 10.4 (slide 19); slide 20 "Conjugacy" remark.

---

## Problem 2

**What the question asks.** Two specific quantitative claims to verify: (a) the posterior of $\theta$ is $N(110.38, 69.23)$; (b) the Bayes factor for $H_0: \theta \le 100$ vs. $H_1: \theta > 100$ is $0.12$. The student must demonstrate both, not merely state them.

**Type of work required.** Both analytical (Normal–Normal conjugate update; ratio of posterior tail probabilities under composite-vs-composite formulation) and R (`pnorm` for the BF).

**Key steps that MUST be shown.**

*(a) Posterior derivation.*
1. Likelihood and prior densities written explicitly with the given variances ($\sigma_0^2 = 100$, $\tau^2 = 225$, $\xi = 100$, $x = 115$, $n = 1$).
2. Posterior mean and variance written in the canonical Normal–Normal form (Lecture 10 Example 10.1, slide 5):
   $$A = \frac{\frac{n\bar X}{\sigma_0^2} + \frac{\xi}{\tau^2}}{\frac{n}{\sigma_0^2} + \frac{1}{\tau^2}}, \qquad B^2 = \frac{1}{\frac{n}{\sigma_0^2} + \frac{1}{\tau^2}}.$$
3. Substitution of the given numbers, with the reviewer independently confirming the arithmetic equals the claimed $(110.38, 69.23)$.

*(b) Bayes factor for composite-vs-composite.*
4. Identification that **both** $H_0$ and $H_1$ are composite, so the simple-vs-simple BF formula does not apply directly. Student must invoke the "General Composite Hypotheses" formula (Lecture 10 slide 15):
   $$\mathrm{BF} = \frac{\int_{\Theta_0} f(X \mid \theta)\rho_0(\theta)\,d\theta}{\int_{\Theta_1} f(X \mid \theta)\rho_1(\theta)\,d\theta}.$$
5. Equivalent (and easier) formulation as the ratio of posterior odds to prior odds:
   $$\mathrm{BF} = \frac{p_0/p_1}{\pi_0/\pi_1}, \qquad p_0 = P(\theta \le 100 \mid X), \quad \pi_0 = P(\theta \le 100).$$
6. Numerical evaluation of both $p_0$ (using the posterior $N(110.38, 69.23)$) and $\pi_0$ (using the prior $N(100, 225)$) via `pnorm`, and the ratio displayed.

**Common mistakes to watch for.**
- Computing the BF as a likelihood ratio at point estimates (this is the simple-vs-simple shortcut from slide 13 and is **wrong** here — both hypotheses are composite).
- Forgetting to divide by the prior odds (i.e., reporting $p_0/p_1$ as the BF).
- Confusing the variance with the standard deviation when calling `pnorm` (the posterior is parameterized as $N(\mu, \sigma^2)$ in lecture; `pnorm` takes $\sigma$, not $\sigma^2$).
- Treating $\sigma_0^2 = 100$ as a standard deviation rather than a variance (or vice versa for the prior).

**Lecture reference.** Lecture 10 Example 10.1 (slides 4–5); slide 15 "General Composite Hypotheses"; Example 10.2 revisited (slide 16) for the BF-via-tail-probability pattern.

---

## Problem 3

**What the question asks.** Specific numerical answers to two decimal places for parts (a)–(d), comparing posterior Bayes estimators (under different priors) with the MLE under two sample sizes.

**Type of work required.** Analytical (closed-form posterior under Beta–Binomial conjugacy; closed-form MLE) with arithmetic; R is optional but acceptable for arithmetic.

**Key steps that MUST be shown.**

*(a) Posterior Bayes estimator under Beta(3,5) prior, $X = 15$, $n = 20$.*
1. Likelihood $L(p \mid x) \propto p^{x}(1-p)^{n-x}$.
2. Conjugate update $p \mid X \sim \mathrm{Beta}(\alpha + x, \beta + n - x)$ stated explicitly (Lecture 10 Example 10.9, slide 38).
3. Posterior mean reported as $E[p \mid X] = \dfrac{\alpha + x}{\alpha + \beta + n}$. Reviewer verifies arithmetic.

*(b) MLE.*
4. $\hat p_{\mathrm{MLE}} = X/n$ stated and evaluated.

*(c) Bayes estimator under Uniform = Beta(1,1) prior.*
5. Same conjugate update with $\alpha = \beta = 1$. Reviewer verifies the new posterior parameters and posterior mean.

*(d) Sample size 200, 150 survivors — three-way comparison.*
6. Recompute MLE, Beta(3,5)-prior Bayes estimator, and Uniform-prior Bayes estimator using the new $(x, n)$.
7. The student must **comment** on the comparison: that as $n$ grows, the posterior mean is dominated by the data and converges toward the MLE regardless of prior. A numerical table or sentence-level comparison is sufficient; absence of any commentary is a deduction.

**Common mistakes to watch for.**
- Reporting the posterior **mode** (= MLE under Uniform) instead of the posterior **mean** as the "Bayes estimator." Under squared error loss the Bayes estimator is the posterior mean (Lecture 10 Corollary 10.1, slide 29).
- Forgetting that Uniform = Beta(1,1) and writing the posterior as Beta$(x, n-x)$ instead of Beta$(x+1, n-x+1)$.
- Mislabeling the prior parameters (writing $\alpha + n - x$ in place of $\beta + n - x$).
- For (d): asserting the priors "don't matter" without showing the numerical convergence; the point of the exercise is to *demonstrate* the shrinkage effect quantitatively.

**Lecture reference.** Lecture 10 Example 10.4 (slide 19); slide 20 Remarks (MLE vs. Bayes; uniform = Beta(1,1)); Example 10.9 (slide 38) for the Beta–Binomial posterior mean form.

---

## Problem 4

**What the question asks.** Derive a specific closed-form expression: the marginal (unconditional) density of $X$ when $X \mid \theta \sim \mathrm{Gamma}(\alpha, 1/\theta)$ and $\theta \sim \mathrm{Gamma}(\beta, 1/\gamma)$.

**Type of work required.** Purely analytical (integration; recognition of a Gamma kernel).

**Key steps that MUST be shown.**
1. Setup: $f(x) = \int_0^\infty f(x \mid \theta)\,\pi(\theta)\,d\theta$ stated explicitly.
2. The integrand written out as a product:
   $$f(x)\,=\,\int_0^\infty \frac{\theta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\theta x}\;\cdot\;\frac{\gamma^\beta}{\Gamma(\beta)}\theta^{\beta-1}e^{-\gamma\theta}\;d\theta.$$
3. Constants pulled outside the integral; $\theta$-terms collected so the integrand is recognized as the kernel of a Gamma density in $\theta$ with new shape and rate. The student must explicitly identify the new shape ($\alpha + \beta$) and new rate ($x + \gamma$) — or whatever they obtain — as a Gamma kernel.
4. Use of the Gamma normalizing constant identity
   $$\int_0^\infty \theta^{a-1} e^{-b\theta}\,d\theta = \frac{\Gamma(a)}{b^{a}}$$
   to evaluate the integral in closed form.
5. Final marginal density written as a clean function of $x$, with the support $x > 0$ stated. Reviewer independently checks that the resulting density integrates to 1 over $x > 0$ (or recognizes it as a known scaled-Beta-prime / compound Gamma form).

**Common mistakes to watch for.**
- Confusing the two Gamma parameterizations (rate vs. scale). The problem uses $\mathrm{Gamma}(\alpha, 1/\theta)$ which is rate-form $\theta$. Errors here cascade.
- Forgetting to combine the $\theta$ exponents correctly: $\theta^\alpha \cdot \theta^{\beta-1} = \theta^{\alpha+\beta-1}$, not $\theta^{\alpha+\beta}$.
- Writing the answer in unsimplified integral form without identifying the closed form (the problem says "obtain", which means "evaluate").
- Dropping the $x^{\alpha-1}$ factor that lives outside the $\theta$-integral.

**Lecture reference.** Lecture 10 slide 3 (marginal density formula $f(X) = \int f(X\mid\theta)\pi(\theta)\,d\theta$); the integration technique is the same kernel-recognition pattern used throughout Examples 10.1, 10.4, and Lecture 11 Examples 7 and 8.

---

## Problem 5

**What the question asks.** Three sub-questions: (a) prior mean (closed form); (b) posterior mean (closed form, given the posterior); (c) prove admissibility of a specific estimator $\hat\theta$ under squared error loss.

**Type of work required.** Purely analytical.

**Key steps that MUST be shown.**

*(a) Prior mean.*
1. Set up $E[\theta] = \int_\beta^\infty \theta \cdot \gamma\beta^\gamma\theta^{-(\gamma+1)}\,d\theta$.
2. Direct integration; the condition $\gamma > 2$ (vs. $\gamma > 1$) is required for the **second** moment, but the first moment exists for $\gamma > 1$. Student should note the convergence condition explicitly.

*(b) Posterior mean.*
3. Likelihood for $Y_i \sim \mathrm{Uniform}(0, \theta)$:
   $$L(\theta \mid \mathbf{y}) = \theta^{-n}\,\mathbf{1}\{\theta \ge z\}, \quad z = \max(y_1,\dots,y_n).$$
   The indicator $\theta \ge z$ is critical and must appear.
4. Posterior $\propto$ prior $\times$ likelihood:
   $$\pi(\theta \mid \mathbf{y}) \propto \theta^{-n} \cdot \theta^{-(\gamma+1)}\,\mathbf{1}\{\theta \ge \max(z, \beta)\} = \theta^{-(n+\gamma+1)}\,\mathbf{1}\{\theta \ge \max(z, \beta)\}.$$
   The lower endpoint of the support must be the **maximum** of $z$ and $\beta$ (because the prior support is $\theta \ge \beta$ and the likelihood support is $\theta \ge z$).
5. Recognition that this is a Pareto-type kernel with the same family as the prior — the student should identify the conjugacy and read off the posterior parameters. Posterior mean computed directly by integration; reviewer verifies the closed form independently.

*(c) Admissibility.*
6. The given estimator $\hat\theta = \frac{\gamma+n}{\gamma+n-1}\max\{z,\beta\}$ must be shown to coincide with the **posterior mean** computed in (b) (i.e., the Bayes estimator under squared error loss; Corollary 10.1(i), Lecture 10 slide 29).
7. Invocation of the lecture result that "a Bayes estimator, when unique, is always admissible" (Lecture 10 slide 26, Bayes Risk slide). The student should argue uniqueness of the posterior mean (e.g., posterior is non-degenerate / proper) and conclude admissibility.

**Common mistakes to watch for.**
- Forgetting the indicator function in the Uniform$(0,\theta)$ likelihood, which gives a wrong posterior support.
- Using only $z$ (not $\max\{z,\beta\}$) as the lower endpoint of the posterior support.
- For (c), trying to prove admissibility from scratch (e.g., by Cramér–Rao or risk inequality) rather than via the "unique Bayes ⇒ admissible" shortcut taught in lecture.
- Confusing admissibility with minimaxity or unbiasedness.
- Asserting admissibility without first showing $\hat\theta = E[\theta \mid \mathbf{Y}]$.

**Lecture reference.** Lecture 10 slide 26 (Bayes Risk; "Bayes estimator, when unique, is always admissible"); Corollary 10.1 (slide 29); Example 10.7 (slides 30–31) for the Uniform$(0,\theta)$ likelihood pattern.

---

## Problem 6

**What the question asks.** A specific decision: which mode of travel (plane vs. car) has smaller expected loss (Bayes risk) under the given crash probabilities and loss table. The student must produce both numerical risks and a recommendation.

**Type of work required.** Purely analytical — straightforward expected-value arithmetic.

**Key steps that MUST be shown.**
1. Identification of the prior as $P(\mathrm{crash})$ vs. $P(\mathrm{no\ crash})$ for **each** mode of travel — the priors differ by mode ($5\times 10^{-6}$ for plane, $10^{-2}$ for car).
2. Bayes risk for each decision computed as expected loss with respect to its mode-specific prior:
   $$r(d_{\mathrm{plane}}) = 10^{6}\cdot(5\times 10^{-6}) + 5\cdot(1 - 5\times 10^{-6}),$$
   $$r(d_{\mathrm{car}}) = 100\cdot(10^{-2}) + 0\cdot(1 - 10^{-2}).$$
3. Numerical comparison and an **explicit conclusion** (which decision is preferred). Reviewer verifies the arithmetic.

**Common mistakes to watch for.**
- Using a single common prior for both modes (the problem gives different crash probabilities for plane and car; these are the relevant priors over $\theta = \{\mathrm{crash}, \mathrm{no\ crash}\}$ under each decision).
- Ignoring the "no crash" loss (e.g., using only the crash row).
- Confusing loss values with probabilities or vice versa.
- Failing to state the conclusion in plain language (the problem asks "which decision...").

**Lecture reference.** Lecture 10 slide 22 (risk function as expected loss); slide 26 (Bayes risk = $\int R_T(\theta)\pi(\theta)\,d\theta$).

---

## Problem 7

**What the question asks.** Three specific deliverables: (a) the Bayes risk of each of $d_1, d_2, d_3$ under the given prior $\{1/8, 3/8, 1/4, 1/4\}$; (b) which decision is Bayes (the minimizer); (c) the minimax decision(s).

**Type of work required.** Purely analytical — finite-sum arithmetic and a max-then-min over the risk table.

**Key steps that MUST be shown.**

*(a) Bayes risk of each $d_j$.*
1. Formula stated explicitly:
   $$r(d_j) = \sum_{i=1}^{4} R(d_j, \theta_i)\,\pi(\theta_i), \qquad j = 1, 2, 3.$$
2. Three numerical values reported, one per column of the risk table. Reviewer verifies arithmetic.

*(b) Bayes decision.*
3. Identification of the column achieving the minimum Bayes risk; ties (if any) must be reported as "all minimizers."

*(c) Minimax decisions.*
4. Computation of the maximum risk for each decision: $\max_i R(d_j, \theta_i)$ — i.e., the column maximum of the risk table.
5. Selection of the decision(s) achieving the minimum of these column maxima.
6. If there are ties at the minimax step, **all** tied decisions must be reported (the question says "decisions", plural).

**Common mistakes to watch for.**
- Computing $\min_i R(d_j, \theta_i)$ (column minima) instead of $\max_i R(d_j, \theta_i)$ for the minimax step.
- Mixing up rows and columns (rows are $\theta$-states, columns are decisions).
- Applying the prior weights in the minimax step (minimax is **not** prior-weighted; it is worst-case over $\theta$).
- Reporting only one minimax decision when multiple decisions tie.

**Lecture reference.** Lecture 10 slide 26 (Bayes risk definition); slide 33 (Minimax estimator definition); Example 10.8 (slides 34–37) for the worked finite-decision-space pattern.

---

## Problem 8

**What the question asks.** Two deliverables: (a) write down the risk function $R_d(\theta)$ for an arbitrary decision rule $d$ on $\{0,1\}$; (b) **show** (i.e., derive, not assert) that the minimax decision satisfies $d_1 = 3/4, d_2 = 1/4$.

**Type of work required.** Purely analytical.

**Key steps that MUST be shown.**

*(a) Risk function.*
1. Computation of $R_d(\theta) = E_\theta[(d(X) - \theta)^2]$ over $X \in \{0, 1\}$:
   $$R_d(\theta) = (d_1 - \theta)^2\,P(X = 1 \mid \theta) + (d_2 - \theta)^2\,P(X = 0 \mid \theta) = (d_1 - \theta)^2\,\theta + (d_2 - \theta)^2\,(1-\theta).$$
   This expression must appear before any minimax argument.

*(b) Minimax derivation.*
2. The student must use **Theorem 10.2** (Lecture 10 slide 33): a Bayes estimator with **constant risk** is minimax. The argument should mirror Example 10.9 (slide 38–39):
   - Place a Beta$(\alpha, \beta)$ prior on $\theta$, derive the Bayes estimator $T(X) = \frac{X+\alpha}{\alpha+\beta+n}$ (with $n=1$ here).
   - Express its risk as a quadratic in $\theta$ and identify the values of $\alpha, \beta$ that make the coefficients of $\theta$ and $\theta^2$ vanish, leaving a constant risk.
3. The two values $d_1 = T(1)$ and $d_2 = T(0)$ must be computed from the chosen $(\alpha, \beta)$ and shown to equal $3/4$ and $1/4$ respectively. The reviewer must verify the algebra independently.
4. Final invocation of Theorem 10.2 to conclude that this Bayes estimator (because it has constant risk) is minimax.

   *Alternative acceptable approach:* Direct minimax argument by writing $R_d(\theta)$ as a polynomial in $\theta$, setting the coefficients of $\theta$ and $\theta^2$ to zero (so the risk is constant in $\theta$), and solving for $d_1, d_2$. Either route should produce the same pair.

**Common mistakes to watch for.**
- Solving $\partial R_d / \partial \theta = 0$ as if $\theta$ were the optimization variable (it is not — $\theta$ is the unknown state of nature; the optimization is over $d_1, d_2$).
- Setting up a Bayes estimator without using the Beta–Binomial conjugacy ($n=1$ gives Bernoulli, but the conjugacy still holds with sample size 1).
- Using $X+\alpha/(\alpha+\beta)$ (forgetting the $+n$ in the denominator).
- Stating "$d_1=3/4, d_2=1/4$" without the constant-risk verification — this is asserting the answer, not deriving it. The problem says "show".
- Confusing minimax with Bayes (the Bayes estimator depends on the prior; the minimax estimator is the one Bayes estimator whose risk is **flat** in $\theta$).

**Lecture reference.** Lecture 10 slide 33 (Theorem 10.2: constant-risk Bayes ⇒ minimax); Example 10.9 (slides 38–39) — this is the direct template for the problem; the student is essentially being asked to redo Example 10.9 with $n = 1$.

---

## Cross-cutting items the reviewer should also check

- **Notation consistency.** The student should use the same symbol (e.g., $\theta$ vs. $p$) consistently within each problem; mixing is acceptable across problems but not within one.
- **Support statements.** Every posterior derivation must state the support of the posterior, especially in Problem 5 where $\max\{z, \beta\}$ matters.
- **Prior/posterior labeling.** Sentences like "the prior posterior is..." are deductions; ensure each density is correctly named.
- **R code provenance.** Where R is used (Problems 1, 2, possibly 3 for arithmetic), a code chunk must appear in the RMD and its output must match the analytical claim. A bare numerical answer without a code chunk or by-hand computation is unverifiable.
- **Scope discipline.** Problem 2 explicitly says "show that..." with target values; the student must reach those targets, not merely set up the problem. Same for Problem 8.
