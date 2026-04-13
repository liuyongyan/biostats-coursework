# P8109 Homework 3 -- Grading Rubric

This rubric specifies **principle-based checking criteria** for each problem. It contains NO pre-computed answers or solutions. The grader must independently verify all derivations.

All problems on this homework are **analytical only** -- no R code is required.

---

## Problem 1: Poisson with T = X1 + 2X2

**Context:** X1, X2 independent Poisson(lambda) -- same parameter lambda for both. T = X1 + 2X2 (note: the coefficient 2 on X2 is what makes T non-sufficient).

### Part (a): Show T is NOT sufficient for lambda (10 points)

- **Requires a specific answer:** Yes -- an explicit expression for P(X1=1, X2=1 | T=3), then a conclusion.
- **Type:** Analytical derivation.
- **Level of detail:** All intermediate steps must be shown -- the joint probability, the marginal P(T=3), and the ratio.

**Relevant lecture material:** Definition 5.1 (Lecture 5, slide 2); Example 5.1 (slides 3-4) for the method of checking sufficiency via conditional distributions.

**Required steps:**

1. **(2 pts)** Identify the event {T=3} = {X1 + 2X2 = 3}. Enumerate all pairs (X1, X2) of non-negative integers satisfying X1 + 2X2 = 3. The student must find all feasible combinations. (Key observation: since the coefficient on X2 is 2, the feasible set is different from X1 + X2 = 3.)
2. **(2 pts)** Compute the numerator P(X1=1, X2=1) using the Poisson pmf for each variable and independence: P(X1=1) * P(X2=1).
3. **(3 pts)** Compute P(T=3) by summing joint probabilities P(X1=a, X2=b) over all feasible pairs from step 1. Each such probability uses the Poisson pmf and independence.
4. **(3 pts)** Form the ratio to get P(X1=1, X2=1 | T=3), simplify, and observe that the result **depends on lambda**. Conclude that by Definition 5.1, T is not sufficient.

**Common mistakes to watch for:**
- Misreading T as X1 + X2 instead of X1 + 2X2. If the student uses T = X1 + X2, the entire solution is based on the wrong statistic.
- Failing to enumerate all (X1, X2) pairs satisfying X1 + 2X2 = 3, or including invalid pairs.
- Not simplifying the conditional probability enough to see that it depends on lambda.
- Concluding non-sufficiency without explicitly showing the lambda-dependence in the final expression.
- Algebra errors when manipulating Poisson pmf products (especially with factorial terms).

---

### Part (b): Sufficient statistic via the factorization theorem (8 points)

- **Requires a specific answer:** Yes -- must state the sufficient statistic for lambda.
- **Type:** Analytical derivation using the factorization theorem.
- **Level of detail:** The joint likelihood must be written out, factored explicitly into g_lambda(T(X)) * h(X), with both parts identified.

**Relevant lecture material:** Theorem 5.1 (Lecture 5, slide 5); Example 5.2 (slide 7) -- the Poisson factorization example, which is directly applicable here.

**Required steps:**

1. **(2 pts)** Write the joint likelihood L(lambda | X1, X2) as the product of two independent Poisson(lambda) pmfs.
2. **(3 pts)** Combine and factor the likelihood. The student must separate terms that depend on lambda (potentially through T(X)) from terms that depend only on the data but not lambda. The student should identify the natural sufficient statistic that emerges from the exponential family structure of the Poisson.
3. **(2 pts)** Explicitly identify g_lambda(T) and h(X), cite the factorization theorem, and state the sufficient statistic.
4. **(1 pt)** Cross-check consistency: the sufficient statistic found here should NOT be T = X1 + 2X2 (which was just shown to be non-sufficient in part (a)).

**Common mistakes to watch for:**
- Proposing T = X1 + 2X2 as sufficient after just showing it is not.
- Not clearly labeling the g_lambda and h components in the factorization.
- Proposing a statistic that is not consistent with the factored form.
- Confusion about why the natural sufficient statistic for two iid Poisson variables works but X1 + 2X2 does not (the weights matter).

---

## Problem 2: Discrete Distribution and Sufficiency

**Context:** Discrete X with pmf table (x=1,2,3; theta = theta1, theta2, theta3); T = 0 if x odd, 1 if x even.

### Part (a): Conditional distribution tables (10 points)

- **Requires a specific answer:** Yes -- two tables of conditional probabilities, one for T=0 and one for T=1.
- **Type:** Direct computation from the given pmf table.
- **Level of detail:** Intermediate steps (marginal probabilities P(T=0; theta) and P(T=1; theta)) should be shown before forming the ratios.

**Relevant lecture material:** Definition 5.1 (Lecture 5, slide 2); Example 5.1 (slides 3-4) for the conditional distribution approach.

**Required steps:**

1. **(2 pts)** Identify which x values correspond to T=0 (odd x: x=1, x=3) and T=1 (even x: x=2).
2. **(4 pts)** For T=0: compute P(T=0; theta_i) = p(1; theta_i) + p(3; theta_i) for each theta. Then compute p(x | T=0; theta_i) = p(x; theta_i) / P(T=0; theta_i) for x = 1 and x = 3, for all three theta values. Present as a table.
3. **(2 pts)** For T=1: since x=2 is the only even value, p(x=2 | T=1; theta) = 1 for all theta. The student should recognize and state this. Present as a (trivial) table.
4. **(2 pts)** Verify internal consistency: conditional probabilities within each table (for a given theta) must sum to 1.

**Common mistakes to watch for:**
- Conditioning on the wrong set of x values for a given T.
- Arithmetic errors when summing probabilities from the table.
- Not presenting results for all three theta values in each table.
- Omitting the T=1 table because it seems trivial.

---

### Part (b): Comment on the nature of T (5 points)

- **Requires a specific answer:** Yes -- must state whether T is sufficient or not, with justification.
- **Type:** Conceptual interpretation based on part (a).

**Relevant lecture material:** Definition 5.1 (Lecture 5, slide 2).

**Required steps:**

1. **(3 pts)** Examine whether the conditional distribution p(x | T; theta) depends on theta. The student must check BOTH the T=0 and T=1 cases. For T=1 the conditional is degenerate and free of theta. For T=0, the student must inspect their table to see whether the conditional probabilities for x=1 and x=3 change across the three theta values.
2. **(2 pts)** State the conclusion about sufficiency. T is sufficient if and only if the conditional distribution given T does not depend on theta for ANY value of T. If even one conditional table depends on theta, T is not sufficient. The student must connect their observation to Definition 5.1.

**Common mistakes to watch for:**
- Drawing a conclusion that contradicts their own computed tables.
- Checking only the T=1 table and claiming sufficiency, while ignoring that the T=0 table might depend on theta.
- Not referencing the definition of sufficiency.

---

## Problem 3: Sufficient Statistic for f(x; theta) = theta^{-1} on (|theta|, 2*theta)

**(12 points)**

- **Requires a specific answer:** Yes -- must state the sufficient statistic.
- **Type:** Analytical derivation using the factorization theorem.
- **Level of detail:** Joint likelihood with indicator functions, conversion to order-statistic conditions, explicit factorization.

**Relevant lecture material:** Theorem 5.1 (Lecture 5, slide 5); Example 5.7 (slide 18) for Uniform(0, theta) where the support depends on theta; Example 6.1 (Lecture 6, slides 3-4) for Uniform(theta, theta+1) where (X_(1), X_(n)) is sufficient.

**Required steps:**

1. **(2 pts)** Write the joint likelihood as a product of individual densities, each with indicator function enforcing |theta| < x_i < 2*theta.
2. **(2 pts)** Address the parameter space. Since the support (|theta|, 2*theta) must be non-empty, we need |theta| < 2*theta, which forces theta > 0. The student should note this (or at least work consistently with theta > 0).
3. **(3 pts)** Convert the indicator constraints to order statistics. The constraints |theta| < x_i for all i and x_i < 2*theta for all i translate to conditions on the sample minimum X_(1) and sample maximum X_(n). The student must correctly determine which order statistic enforces which bound.
4. **(3 pts)** Factor the joint likelihood into g_theta(T) * h(X), identifying the theta-dependent part (involving theta^{-n} and the indicator conditions) and the theta-free part.
5. **(2 pts)** State the sufficient statistic, citing the factorization theorem. Since both bounds depend on theta, the sufficient statistic should capture information from both tails of the sample.

**Common mistakes to watch for:**
- Ignoring the indicator functions entirely and treating this as if it were an exponential family.
- Ignoring the absolute value |theta| in the lower bound.
- Getting the direction of inequalities wrong when converting to order statistics (e.g., confusing min vs. max).
- Proposing the sample mean or sample sum as sufficient (these are NOT sufficient when the support depends on theta and the density is flat).
- Not recognizing that both bounds of the support depend on theta (unlike Uniform(0, theta) where only the upper bound depends on theta).

---

## Problem 4: Rao-Blackwell with N(theta, 1)

**Context:** X1, X2 iid N(theta, 1); d = (X1+X2)/2; d* = E(d | X1).

### Part (a): Expectations of d and d* (5 points)

- **Requires a specific answer:** Yes -- explicit values for E(d) and E(d*).
- **Type:** Analytical computation.

**Relevant lecture material:** Theorem 5.3 part (b) (Lecture 5, slide 19); proof on slide 21 showing E[d*] = E[E(d|T)] = E[d].

**Required steps:**

1. **(2 pts)** Compute E(d) = E[(X1+X2)/2] using linearity of expectation and E(Xi) = theta.
2. **(3 pts)** Compute E(d*). Two valid approaches:
   - **Tower property:** E[d*] = E[E(d | X1)] = E[d] by the law of iterated expectations. The student must invoke this property.
   - **Direct approach:** First compute d* = E[(X1+X2)/2 | X1] explicitly (using independence of X1 and X2), then take E[d*].

**Common mistakes to watch for:**
- Not using or citing the tower property for E(d*).
- When computing E(d | X1), forgetting that X1 is fixed in the conditional expectation while X2 retains its marginal distribution (by independence).

---

### Part (b): Show Var(d*) < Var(d) (7 points)

- **Requires a specific answer:** Yes -- must prove the strict inequality.
- **Type:** Analytical derivation.

**Relevant lecture material:** Theorem 5.3 part (c) (Lecture 5, slides 19, 21-22) and its proof using the variance decomposition (Eve's law).

**Required steps (two valid approaches):**

**Option 1 -- Direct computation:**
1. **(2 pts)** Compute d* = E(d | X1) explicitly. This requires E[(X1+X2)/2 | X1], which uses independence (E(X2 | X1) = E(X2) = theta).
2. **(3 pts)** Compute Var(d) and Var(d*) separately from the explicit expressions.
3. **(2 pts)** Compare and show strict inequality.

**Option 2 -- Eve's law / variance decomposition:**
1. **(3 pts)** Apply Var(d) = Var(E(d|X1)) + E(Var(d|X1)) = Var(d*) + E(Var(d|X1)).
2. **(2 pts)** Show E(Var(d|X1)) > 0 strictly. Compute Var(d|X1) = Var((X1+X2)/2 | X1) = Var(X2/2 | X1) = Var(X2)/4 (by independence), which is strictly positive.
3. **(2 pts)** Conclude Var(d*) = Var(d) - E(Var(d|X1)) < Var(d).

**Common mistakes to watch for:**
- Only showing Var(d*) <= Var(d) without establishing **strict** inequality.
- When computing E(X2 | X1), not using independence correctly.
- Incorrectly applying the Rao-Blackwell theorem here (the theorem requires a sufficient statistic -- X1 is NOT sufficient, but the variance inequality still holds by Eve's law; the issue comes up in part (c)).

---

### Part (c): Why d* cannot be used as improved estimator of theta (6 points)

- **Requires a specific answer:** Yes -- must give a clear conceptual explanation.
- **Type:** Conceptual reasoning.

**Relevant lecture material:** Theorem 5.3 and Remarks (Lecture 5, slides 19-20). The Rao-Blackwell theorem specifically requires T to be a SUFFICIENT statistic.

**Required content:**

1. **(3 pts)** Identify the key issue: X1 alone is NOT a sufficient statistic for theta in the N(theta, 1) model with two observations. The sufficient statistic is the sum (or mean) of both observations (by the factorization theorem / exponential family, see Lecture 5, Example 5.5). Since d* = E(d | X1) conditions on a non-sufficient statistic, the Rao-Blackwell theorem does not apply.
2. **(3 pts)** Explain the practical consequence: d* is a function of X1 alone, effectively discarding X2. It ignores information about theta contained in X2. While Var(d*) < Var(d), this does not make d* a superior estimator because there exist other estimators (based on the sufficient statistic) that have even lower variance. The sample mean d itself, when viewed as a function of the sufficient statistic, already uses all the data optimally in a different sense.

**Common mistakes to watch for:**
- Failing to identify that X1 is not sufficient as the core reason.
- Incorrectly claiming d* is biased (it IS unbiased, as shown in part (a)).
- Claiming the problem is that Var(d*) is not actually smaller (it IS smaller, as shown in part (b) -- that is not the issue).
- Vague answers that do not specifically mention sufficiency or the Rao-Blackwell conditions.

---

## Problem 5: Show Uniform(theta, theta+1) Family Is Not Complete

**(15 points)**

- **Requires a specific answer:** Yes -- must exhibit a specific non-zero function h with E[h] = 0 for all theta.
- **Type:** Proof by construction (counterexample to completeness).
- **Level of detail:** The function must be exhibited, E[h] = 0 verified for all theta, and h shown to be non-zero.

**Relevant lecture material:** Definition 6.2 (Lecture 6, slide 9); Example 6.4 part (i) (slide 10) for the strategy of exhibiting a counterexample. Also Example 6.1 (slides 3-4) for the Uniform(theta, theta+1) distribution and its sufficient statistics.

**Required steps:**

1. **(3 pts)** State the definition of completeness (Definition 6.2): T is complete for theta if E[h(T)] = 0 for all theta implies h(T) = 0 a.s. To show the family is NOT complete, one must find a function h that is not identically zero but has zero expectation for all theta.
2. **(3 pts)** Exhibit a specific non-zero function h(x). The student has freedom in choice. Valid approaches include:
   - Periodic functions with period 1 and zero integral over one period (e.g., h(x) = sin(2*pi*x), or h(x) = x - floor(x) - 1/2).
   - Any function satisfying: integral from theta to theta+1 of h(x) dx = 0 for all theta. Differentiating this in theta gives h(theta+1) - h(theta) = 0, so h must be periodic with period 1. Thus any non-constant function with period 1 and zero mean over one period works.
3. **(5 pts)** Verify that E[h(X)] = 0 for ALL theta. This requires computing the integral from theta to theta+1 of h(x) * 1 dx and showing it equals zero for every theta in the parameter space.
4. **(2 pts)** Verify that h is not identically zero (exhibit specific values where h is nonzero, or argue from the form of h).
5. **(2 pts)** Conclude: since a non-zero function h exists with E[h(X)] = 0 for all theta, the Uniform(theta, theta+1) family is not complete by the definition of completeness.

**Common mistakes to watch for:**
- Proposing a function that works for some but not all theta.
- Choosing a function that is actually identically zero.
- Not verifying E[h] = 0 for ALL theta (only checking finitely many values is insufficient).
- Attempting to use Theorem 6.2 (exponential family completeness) -- Uniform(theta, theta+1) is NOT an exponential family, so that theorem does not apply. However, the non-applicability of Theorem 6.2 does not prove non-completeness; a constructive counterexample is required.
- Confusing "not complete" with "not sufficient" -- these are different concepts.
- Working with n > 1 observations and a sufficient statistic like (X_(1), X_(n)) is acceptable, but the student must handle the more complex joint distribution.

---

## Problem 6: UMVUE of p^r * q^s for Binomial

**Context:** X1, ..., Xn iid Binomial(m, p); T = X1 + ... + Xn; estimate p^r * q^s where q = 1-p, r and s are known positive integers with r + s < mn.

### Part (a): Distribution of T (6 points)

- **Requires a specific answer:** Yes -- must state the distribution of T with correct parameters.
- **Type:** Derivation using MGFs (as suggested) or another valid method.
- **Level of detail:** The derivation must be shown, not just the result.

**Relevant lecture material:** Lecture 5, slide 10 (exponential family sufficient statistics); general probability theory on sums of independent Binomials.

**Required steps:**

1. **(2 pts)** Write the MGF of a single Xi ~ Binomial(m, p).
2. **(2 pts)** Use independence to write the MGF of T = X1 + ... + Xn as the product of n individual MGFs.
3. **(2 pts)** Recognize the resulting MGF as that of a known distribution and state the distribution of T clearly with both parameters. (Alternatively, argue via the interpretation of each Binomial(m, p) as a sum of m independent Bernoulli(p) trials, so T is the sum of nm Bernoulli(p) trials.)

**Common mistakes to watch for:**
- Incorrect MGF formula for the Binomial.
- Incorrectly computing the parameters of the resulting distribution (e.g., getting Binomial(nm, np) instead of Binomial(nm, p) -- the p parameter does NOT get multiplied by n).
- Claiming the result without showing any derivation or justification.
- If using the Bernoulli decomposition argument, not justifying why the sum of independent Bernoullis is Binomial.

---

### Part (b): Show U(T) is unbiased for p^r * q^s (10 points)

- **Requires a specific answer:** Yes -- must derive the given formula for U(T) by solving E[U(T)] = p^r * q^s.
- **Type:** Algebraic derivation involving combinatorial identities.
- **Level of detail:** All key combinatorial manipulations must be shown explicitly.

**Relevant lecture material:** Example 6.6 (Lecture 6, slide 18) for a similar technique of finding an unbiased estimator as a function of a complete sufficient statistic.

**Required steps:**

1. **(2 pts)** Write E[U(T)] using the pmf of T (from part (a)): E[U(T)] = sum over t of U(t) * C(mn, t) * p^t * q^{mn-t}, and set this equal to p^r * q^s.
2. **(3 pts)** Manipulate the right-hand side p^r * q^s to express it in a form compatible with the left-hand side. A standard approach: write p^r * q^s = p^r * q^s * 1 = p^r * q^s * (p + q)^{mn-r-s}, then expand (p + q)^{mn-r-s} using the binomial theorem.
3. **(3 pts)** Match the resulting expansion with the left-hand side term by term. This involves re-indexing the sum and comparing coefficients of p^t * q^{mn-t} on both sides. The key combinatorial identity step is to relate C(mn, t) * U(t) to C(mn-r-s, t-r) (or equivalent manipulations).
4. **(2 pts)** Arrive at the formula U(t) = C(mn-r-s, t-r) / C(mn, t) for t = r, r+1, ..., mn-s, and U(t) = 0 otherwise. Justify the range: the binomial coefficient C(mn-r-s, t-r) requires 0 <= t-r <= mn-r-s, i.e., r <= t <= mn-s.

**Common mistakes to watch for:**
- Errors in binomial coefficient algebra or re-indexing sums.
- Not correctly determining the range of T for which U(T) is nonzero.
- Skipping the key combinatorial step without justification.
- Circular reasoning: simply plugging in the given U(T) and verifying is acceptable (the problem says "show that"), but the student must still carry out the verification rigorously.

---

### Part (c): Argue U(T) is the UMVUE (6 points)

- **Requires a specific answer:** Yes -- must provide a brief but complete argument.
- **Type:** Conceptual argument citing appropriate theorems.

**Relevant lecture material:** Theorem 6.2 (Lecture 6, slide 13) on completeness for exponential families; Theorem 6.3 / Lehmann-Scheffe (slide 14); Example 6.4(ii) (slides 10-12) for completeness of the Binomial sum.

**Required steps:**

1. **(2 pts)** Establish that T is a sufficient statistic for p. This follows from the factorization theorem applied to the joint Binomial likelihood, or from the exponential family structure (Lecture 5, slide 10, Example 5.3).
2. **(2 pts)** Establish that T is a complete statistic. The student should argue that the Binomial(mn, p) family (with p in (0,1)) is a one-parameter exponential family of full rank, so by Theorem 6.2, T is complete. Alternatively, the student may cite or reproduce the polynomial argument from Example 6.4(ii) (Lecture 6, slides 10-12): E[h(T)] = 0 for all p implies a polynomial in p/(1-p) is identically zero, forcing all coefficients (and hence h) to be zero.
3. **(2 pts)** Apply the Lehmann-Scheffe Theorem (Theorem 6.3): U(T) is an unbiased estimator of p^r * q^s (shown in part (b)) that is a function of the complete sufficient statistic T, hence U(T) is the UMVUE.

**Common mistakes to watch for:**
- Claiming UMVUE from sufficiency alone, without establishing completeness.
- Invoking only the Rao-Blackwell theorem without Lehmann-Scheffe -- Rao-Blackwell gives variance reduction but does not by itself guarantee UMVUE or uniqueness.
- Confusing sufficiency with completeness.
- Not mentioning the Lehmann-Scheffe theorem (by name or statement).
- Failing to identify the exponential family structure or otherwise justify completeness.

---

## Point Allocation Summary

| Problem | Points | Key Concept |
|---------|--------|-------------|
| 1(a) | 10 | Conditional distribution, sufficiency definition |
| 1(b) | 8 | Factorization theorem |
| 2(a) | 10 | Conditional distribution (discrete) |
| 2(b) | 5 | Sufficiency assessment |
| 3 | 12 | Factorization with indicator functions, order statistics |
| 4(a) | 5 | Iterated expectation / tower property |
| 4(b) | 7 | Variance comparison, Eve's law |
| 4(c) | 6 | Rao-Blackwell requires sufficient statistic |
| 5 | 15 | Completeness definition, constructive counterexample |
| 6(a) | 6 | MGF / reproductive property of Binomial |
| 6(b) | 10 | Unbiased estimator derivation via combinatorics |
| 6(c) | 6 | Lehmann-Scheffe theorem |
| **Total** | **100** | |

---

## General Grading Principles

1. **Notation and rigor:** Deduct 1 point per problem for sloppy or undefined notation (e.g., using the same symbol for different quantities, not defining new variables).
2. **Theorem citations:** Students should cite relevant definitions and theorems (by name or number) when applying them. Deduct 1 point if a major theorem is applied without identification.
3. **Logical flow:** Each derivation should proceed in a clear, logical sequence. Gaps in reasoning (skipped steps that are not obvious) should result in partial credit rather than full credit.
4. **Correct but different approaches:** Accept any mathematically valid approach, even if it differs from the lecture method, provided the reasoning is complete and correct.
5. **Partial credit:** Award partial credit for correct intermediate steps even if the final conclusion is wrong, provided the student demonstrates understanding of the relevant concepts.
6. **Cross-part consistency:** Within multi-part problems, answers across sub-parts must be mutually consistent. Flag and deduct for contradictions (e.g., finding a statistic is not sufficient in one part but claiming it is in another).
7. **All derivations are analytical** -- no R code is required or expected for any problem on this homework.
