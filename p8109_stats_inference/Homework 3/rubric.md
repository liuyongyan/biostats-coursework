# P8109 Homework 3 -- Grading Rubric

This rubric specifies **principle-based checking criteria** for each problem. The grader must independently verify the student's work against the stated principles; no pre-computed answers are provided.

---

## Problem 1: Poisson Sufficiency

**Context:** X1, X2 iid Poisson(lambda), T = X1 + X2.

### Part (a): Show T is NOT sufficient for lambda

- **Type:** Analytical derivation leading to a specific conditional probability expression, then a conclusion.
- **Required approach:** The student must compute P(X1=1, X2=2 | T=3) using the definition of conditional probability:
  - Numerator: joint probability P(X1=1, X2=2) using the Poisson pmf.
  - Denominator: P(T=3), recognizing that T = X1+X2 follows Poisson(2*lambda) (sum of independent Poissons).
- **Key criterion for the sufficiency conclusion:** The student must evaluate the resulting conditional probability and observe whether it depends on lambda. According to Definition 5.1 (Lecture 5, slide 2), T is sufficient only if the conditional distribution f(x | T) does not depend on the parameter. The student must clearly state whether their computed expression depends on lambda and draw the correct conclusion about sufficiency.
- **Relevant lecture material:** Definition 5.1 (Lecture 5, slide 2); Example 5.1 (slides 3-4) for the method of checking sufficiency via conditional distributions.
- **Common mistakes to watch for:**
  - Incorrectly computing P(T=3) -- e.g., using Poisson(lambda) instead of Poisson(2*lambda).
  - Failing to simplify the conditional probability fully before concluding.
  - Stating the wrong conclusion about dependence on lambda.
  - Confusing "not sufficient" with "sufficient" (the problem says "show that T is NOT sufficient" -- but the grader should verify the student's algebra independently to confirm whether the conclusion follows from the computation).

### Part (b): Find a sufficient statistic via the factorization theorem

- **Type:** Analytical derivation using the factorization theorem.
- **Required approach:** The student must write the joint likelihood L(lambda | X1, X2) = product of two Poisson pmfs, then factor it into g_lambda(T(X)) * h(X) per Theorem 5.1 (Lecture 5, slide 5).
- **Key criteria:**
  - The joint likelihood must be written out correctly as a product of two Poisson pmfs.
  - The factorization must explicitly identify which part is g_theta(T(X)) (depends on lambda only through the statistic) and which part is h(X) (does not depend on lambda).
  - The student must clearly state what the sufficient statistic is, based on the factorization.
- **Relevant lecture material:** Theorem 5.1 (Lecture 5, slide 5); Example 5.2 (slide 7) is directly analogous (Poisson case with n observations).
- **Common mistakes to watch for:**
  - Not explicitly identifying the g and h components.
  - Proposing a statistic that is not consistent with their own factorization.
  - Note: the student's answer to (a) and (b) should be logically consistent -- the sufficient statistic found in (b) should differ from T = X1 + X2 if (a) correctly shows T is not sufficient, OR the student should reconcile any apparent contradiction. The grader must check this consistency independently.

---

## Problem 2: Discrete Distribution Sufficiency

**Context:** Discrete X with given pmf table; T = 0 if x odd, 1 if x even.

### Part (a): Conditional distribution tables

- **Type:** Direct computation; the answer is a table of conditional probabilities.
- **Required approach:** For each value of T (0 and 1), compute p(x | T; theta) = P(X=x, T=t) / P(T=t) for each x and each theta.
  - When T=0 (x is odd): the relevant x values are 1 and 3. Sum the probabilities for x=1 and x=3 to get P(T=0; theta), then divide each by this sum.
  - When T=1 (x is even): the relevant x value is 2. Since x=2 is the only even value, P(X=2 | T=1; theta) should be computed.
- **Key criteria:**
  - The conditional probabilities must be computed correctly from the given pmf table using the definition of conditional probability.
  - Two separate tables should be presented (one for T=0, one for T=1), each showing the conditional distribution across all relevant x values for all three theta values.
  - The conditional probabilities within each table (for a given theta) must sum to 1.
- **Relevant lecture material:** Definition 5.1 (Lecture 5, slide 2) for the concept of conditioning.
- **Common mistakes to watch for:**
  - Conditioning on the wrong set of x values for a given T.
  - Arithmetic errors in dividing probabilities.
  - Not presenting results for all three theta values.

### Part (b): Comment on the nature of T

- **Type:** Conceptual interpretation based on the result of part (a).
- **Required approach:** The student must examine the conditional distributions from part (a) and determine whether they depend on theta. Per Definition 5.1, if the conditional distribution p(x | T; theta) does not depend on theta, then T is sufficient.
- **Key criteria:**
  - The student must explicitly state whether the conditional distributions depend on theta or not.
  - The conclusion about sufficiency (or lack thereof) must be logically consistent with the tables computed in part (a).
  - The student should reference the definition of sufficiency to justify their conclusion.
- **Common mistakes to watch for:**
  - Drawing a conclusion that contradicts their own computed tables.
  - Not clearly linking the observation to the definition of sufficiency.

---

## Problem 3: Sufficient Statistic for Uniform-like Distribution

**Context:** X1,...,Xn iid from f(x; theta) = theta^{-1} for |theta| < x < 2*theta, 0 elsewhere.

### Requirements

- **Type:** Analytical derivation using the factorization theorem.
- **Required approach:**
  - Write the joint likelihood as a product of the individual densities.
  - Crucially, handle the support constraint |theta| < x_i < 2*theta for all i. This translates into indicator functions that constrain theta based on the data.
  - The student must recognize that the support depends on theta, so the indicator functions are part of g_theta and contribute to identifying the sufficient statistic.
  - Apply the factorization theorem (Theorem 5.1) to identify g_theta(T(X)) and h(X).
- **Key criteria:**
  - The student must correctly express the support constraints in terms of order statistics or min/max of the sample.
  - The constraint |theta| < x_i < 2*theta must be translated into conditions on theta involving the sample minimum X_(1) and/or maximum X_(n).
  - The factorization must correctly separate the theta-dependent part (including the indicator/support constraints) from the theta-free part.
  - The proposed sufficient statistic must be consistent with the factorization.
- **Relevant lecture material:** Theorem 5.1 (Lecture 5, slide 5); Example 5.7 (slide 18) for a case where the support depends on theta and the sufficient statistic involves an order statistic. Also Example 6.1 (Lecture 6, slides 3-4) for the Uniform case with theta-dependent support.
- **Common mistakes to watch for:**
  - Ignoring the indicator functions / support constraints entirely.
  - Incorrectly converting the constraint |theta| < x < 2*theta into conditions on order statistics. Note the absolute value: |theta| means the constraint depends on the sign of theta. The student should carefully handle the case theta > 0 (where |theta| = theta, so theta < x < 2*theta). Grader should verify whether the problem implicitly restricts theta > 0 (since f is a valid density only when theta > 0, because the support length is 2*theta - |theta| = theta when theta > 0, matching the normalizing constant theta^{-1}).
  - Proposing only a single order statistic when the sufficient statistic may need to be a pair (or vice versa) -- the grader should verify independently based on the factorization.
  - Not recognizing that this distribution's support depends on theta, making it non-exponential-family-like.

---

## Problem 4: Rao-Blackwell Theorem Application

**Context:** X1, X2 iid N(theta, 1); d = (X1+X2)/2; d* = E(d | X1).

### Part (a): Expectations of d and d*

- **Type:** Direct computation of expectations.
- **Required approach:**
  - E(d): Compute E[(X1+X2)/2] using linearity of expectation and the fact that E(Xi) = theta.
  - E(d*): Apply the tower property (law of iterated expectations): E[d*] = E[E(d | X1)] = E(d).
- **Key criteria:**
  - Both expectations must be computed and stated explicitly.
  - The tower property must be used (at least implicitly) for E(d*).
  - Both results should be clearly stated.
- **Relevant lecture material:** Theorem 5.3 part (b) (Lecture 5, slide 19) -- that d* = E(d|T) is unbiased if d is unbiased; Proof on slide 21.
- **Common mistakes to watch for:**
  - Not recognizing or applying the tower property for E(d*).
  - Unnecessary computation when the tower property gives the answer directly.

### Part (b): Show Var(d*) < Var(d)

- **Type:** Analytical proof of strict variance inequality.
- **Required approach:** The student should either:
  - (Option 1) Compute d* = E(d | X1) explicitly by evaluating E[(X1+X2)/2 | X1], then compute Var(d*) and Var(d) separately and compare.
  - (Option 2) Use the variance decomposition formula: Var(d) = Var(E(d|X1)) + E(Var(d|X1)) = Var(d*) + E(Var(d|X1)), and show that E(Var(d|X1)) > 0, giving strict inequality.
- **Key criteria:**
  - The student must show strict inequality (not just <=).
  - If using direct computation, both variances must be correctly derived.
  - If using the variance decomposition (as in the proof of Theorem 5.3, Lecture 5, slides 21-22), the student must verify that the remainder term E[Var(d|X1)] is strictly positive (not just non-negative).
- **Relevant lecture material:** Theorem 5.3 part (c) and its proof (Lecture 5, slides 19, 21-22). The variance decomposition Var(d) = Var(d*) + E[Var(d|T)] is the key identity.
- **Common mistakes to watch for:**
  - Only showing Var(d*) <= Var(d) without establishing strictness.
  - Incorrectly computing E(d | X1) -- this requires computing E(X2 | X1), and the student must use the independence of X1 and X2.
  - Confusing conditioning on X1 (not sufficient) with conditioning on a sufficient statistic.

### Part (c): Why d* cannot be used as improved estimator

- **Type:** Conceptual explanation.
- **Required approach:** The student must explain why, despite having lower variance, d* is not a useful "improved" estimator of theta. The key issue is that X1 alone is NOT a sufficient statistic for theta (the sufficient statistic for N(theta,1) with n=2 is X1+X2, not X1 alone). The Rao-Blackwell theorem (Theorem 5.3) guarantees improvement only when conditioning on a SUFFICIENT statistic. Conditioning on a non-sufficient statistic (X1) loses information and produces an estimator (d*) that depends only on X1, discarding X2 entirely. The resulting d* is simply a function of X1 alone and thus ignores half the data.
- **Key criteria:**
  - The student must identify that X1 is NOT a sufficient statistic for theta.
  - The student must connect this to the Rao-Blackwell theorem's requirement of conditioning on a sufficient statistic.
  - The student should note that d* is effectively just a function of X1, meaning it throws away the information in X2.
  - The student should recognize that while Var(d*) < Var(d), the estimator d* is inferior to d in the sense that d uses all the data while d* does not (and there exist other estimators based on all the data that dominate d*).
- **Relevant lecture material:** Theorem 5.3 and Remarks (Lecture 5, slides 19-20). The Rao-Blackwell theorem specifically requires T to be sufficient.
- **Common mistakes to watch for:**
  - Failing to mention sufficiency as the key issue.
  - Incorrectly claiming that d* is biased (it is unbiased by the tower property).
  - Vague answers that do not precisely identify why the Rao-Blackwell framework does not apply here.

---

## Problem 5: Incompleteness of Uniform(theta, theta+1)

### Requirements

- **Type:** Proof by exhibiting a counterexample -- find a non-zero function h(T) of a sufficient statistic T such that E[h(T)] = 0 for all theta.
- **Required approach:** Per Definition 6.2 (Lecture 6, slide 9), to show incompleteness one must find a function h(T) that is not identically zero but has E[h(T)] = 0 for all theta. The student must:
  1. Identify or propose a sufficient statistic T for this family (see Example 6.1 in Lecture 6 for the Uniform(theta, theta+1) case -- the sufficient statistic is (X_(1), X_(n))).
  2. Construct a non-trivial function h of T (or of X) such that its expectation is zero for all theta.
  3. Verify that E[h(T)] = 0 for all theta by computing the expectation under the Uniform(theta, theta+1) distribution.
  4. Verify that h is not identically zero.
- **Key criteria:**
  - The student must produce an explicit function h and show it satisfies both conditions (expectation zero for all theta; function not identically zero).
  - The computation of E[h(T)] = 0 must be valid for ALL theta (not just specific values).
  - The function h must be a function of a sufficient statistic (or equivalently, the student can work with any statistic and note that completeness of the family implies completeness of any sufficient statistic).
- **Relevant lecture material:** Definition 6.2 (Lecture 6, slide 9); Example 6.4 part (i) (slide 10) for a demonstration of showing a statistic is NOT complete by exhibiting a counterexample. Example 6.1 (slides 3-4) for the sufficient statistic of the Uniform(theta, theta+1) family.
- **Common mistakes to watch for:**
  - Failing to verify that E[h(T)] = 0 for ALL theta (not just one value).
  - Choosing a function h that is actually identically zero (trivial).
  - Not working with a sufficient statistic (though one can also show the family itself is not complete, which implies no complete sufficient statistic exists).
  - Confusing "not complete" with "not sufficient."
  - A common and valid approach involves a function of X_(n) - X_(1) (the range), since for Uniform(theta, theta+1) the range has a distribution that does not depend on theta in a way that "fills up" the parameter space enough for completeness.

---

## Problem 6: UMVUE of p^r * q^s for Binomial

**Context:** X1,...,Xn iid Binomial(m, p); estimate p^r * q^s where q = 1-p, r+s < mn.

### Part (a): Distribution of T = X1 + ... + Xn

- **Type:** Derivation using moment generating functions (or convolution).
- **Required approach:** The student must derive the distribution of the sum of n independent Binomial(m, p) random variables. The standard approach is:
  - Write the mgf of a single Binomial(m, p) variable.
  - Use independence to write the mgf of T as the product of n such mgfs.
  - Recognize the resulting mgf as that of a known distribution.
  - Alternatively, argue via the reproductive property of the Binomial: the sum of independent Binomial(m, p) variables is Binomial(nm, p).
- **Key criteria:**
  - The mgf (or other method) must be correctly computed.
  - The final distribution of T must be clearly identified with its parameters.
  - The derivation should be complete, not just stated as a fact.
- **Relevant lecture material:** This draws on probability prerequisites (mgf techniques). The result is used in the subsequent parts.
- **Common mistakes to watch for:**
  - Incorrect mgf formula for the Binomial distribution.
  - Claiming the sum is Binomial without justification (the problem says "by using mgfs or otherwise," so some derivation is expected).

### Part (b): Show U(T) is unbiased for p^r * q^s

- **Type:** Algebraic verification that E[U(T)] = p^r * q^s.
- **Required approach:** The student must compute E[U(T)] using the distribution of T found in part (a) and show it equals p^r * q^s. This involves:
  - Writing out E[U(T)] = sum over T from r to mn-s of U(T) * P(T=t).
  - Substituting the given formula for U(T) = C(mn-s-r, T-r) / C(mn, T).
  - Substituting the Binomial(mn, p) pmf for P(T=t).
  - Performing algebraic manipulations involving binomial coefficients (such as Vandermonde's identity or the technique of splitting binomial coefficients) to show the sum equals p^r * q^s.
- **Key criteria:**
  - The expectation must be written out explicitly as a sum.
  - The algebraic manipulation of binomial coefficients must be clearly shown.
  - The key identity or technique used to evaluate the sum must be stated or derived (e.g., recognizing a hypergeometric-type sum, using the identity C(mn, T) = C(mn, r) * ... or the Vandermonde convolution, or splitting the Binomial(mn, p) pmf into components).
  - The final result E[U(T)] = p^r * q^s must be clearly established.
- **Relevant lecture material:** Example 6.6 (Lecture 6, slide 18) for a similar computation of conditioning to find an unbiased estimator based on a sufficient statistic.
- **Common mistakes to watch for:**
  - Errors in binomial coefficient algebra.
  - Not accounting for the range of summation (T = r to mn-s) correctly.
  - Skipping critical steps in the combinatorial argument.
  - Not verifying that U(T) is indeed a function of the sufficient statistic T only.

### Part (c): Argue U(T) is the UMVUE

- **Type:** Brief conceptual argument invoking the Lehmann-Scheffe theorem.
- **Required approach:** The student must argue that:
  1. T = X1 + ... + Xn is a sufficient statistic for p (from the factorization theorem applied to the Binomial family).
  2. T is a complete statistic (because the Binomial(mn, p) family is a one-parameter exponential family with full rank -- invoke Theorem 6.2, Lecture 6, slide 13).
  3. U(T) is an unbiased estimator of p^r * q^s that is a function of the complete sufficient statistic T.
  4. By the Lehmann-Scheffe theorem (Theorem 6.3, Lecture 6, slide 14), U(T) is the UMVUE.
- **Key criteria:**
  - Sufficiency of T must be established or referenced.
  - Completeness of T must be established, citing the exponential family result (Theorem 6.2) or proving it directly (as in Example 6.4(ii), Lecture 6, slides 11-12).
  - The Lehmann-Scheffe theorem must be explicitly invoked by name or theorem number.
  - The argument must clearly connect all three ingredients: complete, sufficient, and unbiased.
- **Relevant lecture material:** Theorem 6.2 (Lecture 6, slide 13) for completeness of exponential family; Theorem 6.3 (Lecture 6, slide 14) for Lehmann-Scheffe; Example 6.4(ii) (slides 11-12) for completeness of the Binomial sum.
- **Common mistakes to watch for:**
  - Claiming UMVUE without establishing completeness.
  - Not mentioning sufficiency.
  - Invoking Rao-Blackwell alone (which gives minimum variance among Rao-Blackwellized estimators, but not UMVUE without completeness).
  - Confusing sufficiency with completeness.
  - Failing to verify that the Binomial family is an exponential family (or failing to establish completeness by another route).

---

## General Grading Notes

- **Notation and rigor:** Students should use proper notation and clearly define any symbols introduced. Steps should be logically connected.
- **Referencing theorems:** Students should cite relevant definitions and theorems (by name or number) when applying them. Simply writing a formula without identifying the theorem is insufficient for full credit.
- **Logical consistency:** Within multi-part problems, answers across sub-parts should be mutually consistent. Flag contradictions.
- **All derivations are analytical** unless otherwise noted -- no R code is required for any problem on this homework.
