# P8106 Data Science II — Homework 4 Grading Rubric

**Topic:** Tree-based methods and ensemble methods (regression tree, random forest, boosting; classification tree, boosting)

**Source materials consulted:**
- Homework PDF: `/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8106_data_science_II/Homework 4/DSII_HW4_2026.pdf`
- Lecture 10 (trees): `/Users/yliu/Desktop/Columbia - Biostatistics/Data Science II/Lecture 10/DSII_L10_2026.pdf` — covers regression trees (pp. 2–20), pruning / cost-complexity / CV selection of alpha (pp. 15–20), classification trees with Gini / cross-entropy / classification error (pp. 21–26)
- Lecture 8 R code reference: `/Users/yliu/Desktop/Columbia - Biostatistics/Data Science II/Lecture 8/L8_2026.Rmd` — `caret::train` workflow, `trainControl`, CV, variable importance via `vip`, resampling comparisons
- Lecture 9 R code reference: `/Users/yliu/Desktop/Columbia - Biostatistics/Data Science II/Lecture 9/L9_2026.Rmd` — classification `trainControl` with `twoClassSummary`, `classProbs = TRUE`, `pROC`, `confusionMatrix`

**Important note for the reviewer:** As of this rubric's preparation, only Lecture 10 (on single trees) is available in the lecture archive; no lecture R code file exists yet for trees/RF/boosting. The homework assumes students have been exposed to bagging / random forest / boosting via a later lecture (likely Lecture 11). Because that lecture's R code file is not available here, the reviewer should verify that the student uses function / argument conventions consistent with the course style demonstrated in L8/L9 Rmd files (caret-based workflow, same `trainControl` idiom, `vip` / `varImp` for importance, `set.seed` before tuning, etc.). If the class-specific R code for trees/ensembles becomes available, compare line by line; otherwise accept any standard implementation that follows the course's `caret` / `tidymodels` conventions.

---

## General expectations (all problems)

- The student must submit a reproducible document (Rmd + PDF/HTML) with clearly labeled code chunks and narrative.
- A `set.seed(...)` call must appear before any data split and before any CV/tuning. The same seed should be reused consistently within a problem so results are reproducible.
- Data splitting must use the exact proportions specified (80/20 in Problem 1, 70/30 in Problem 2). The split should be done once per dataset and reused across sub-parts. Acceptable split idioms include `rsample::initial_split()` (as used in L8/L9 Rmd) or `caret::createDataPartition()`.
- Categorical predictors in `College.csv` / `auto.csv` must be handled appropriately (factor conversion where needed). For `auto.csv`, `mpg_cat` must be treated as a factor for classification.
- When building models with `caret::train`, `trainControl` should specify CV / repeatedcv with a sensible fold count (5 or 10). For classification, `classProbs = TRUE` and `summaryFunction = twoClassSummary` with `metric = "ROC"` follows the course convention from L9 (optional but recommended).
- Reported "test error" must be computed on the held-out test set only — NOT CV error, NOT training error. The reviewer should verify the student clearly distinguishes the two.
- For regression problems, test error = MSE or RMSE on the test set. Either is acceptable as long as clearly labeled.
- For classification problems, test performance should include at least one of: test classification error / accuracy, confusion matrix, test AUC.
- Any random-forest / boosting tuning must specify a tuning grid (not a single default) and the grid should be chosen reasonably (not degenerate, not absurdly large). Tuning should be by CV.
- **Do NOT penalize minor differences** in numerical answers that arise from differing seeds, CV fold counts, or tuning grids — the reviewer must independently verify that the methodology is correct, not that the numbers match a pre-computed answer key.

---

## Problem 1 — Regression on the College data (response: `Outstate`)

The problem explicitly requires an 80/20 split. Verify that the student states and implements this split.

### 1(a) Regression tree

**Question type:** Asks for a built model AND a plot. No specific numerical answer is demanded, but interpretation is implied by "create a plot of the tree."

**Expected deliverable:**
- R code that fits a regression tree on the training set.
- A readable tree plot.

**Accepted approaches:**
- `rpart::rpart()` + `rpart.plot::rpart.plot()` or `partykit::as.party()` plot, OR
- `caret::train(..., method = "rpart")` with a `cp` tuning grid via CV, then plotting the `finalModel`, OR
- `tree::tree()` + `plot()` + `text()`.

Either full-grown-then-pruned (with cost-complexity CV) or `caret` CV tuning of `cp` is acceptable. The latter is more consistent with the L8/L9 caret workflow.

**Key steps the reviewer MUST check line by line:**
1. Model is fit ONLY on the training data (not full data).
2. A tuning / pruning step is present — either cost-complexity pruning chosen by CV (`printcp`, `plotcp`, `prune`) or caret `train` with a `cp` (or `maxdepth`) grid and CV.
3. A tree plot is produced and is readable (labeled splits, terminal values).
4. Brief narrative identifies at least the root split or the most important splitting variables (light interpretation).

**Common mistakes to watch for:**
- Fitting on the entire dataset instead of training data.
- Reporting only a textual tree summary with no plot.
- Failing to prune / tune — just accepting `rpart` defaults without discussing why.
- Confusing classification vs regression (`method = "anova"` is the regression setting for `rpart`; `method = "class"` would be wrong here).
- Treating this as a classification problem by mistake — `Outstate` is continuous.

**NOT required:**
- A numerical test error is not asked for in (a).
- A variable importance plot is not asked for in (a).

---

### 1(b) Random forest

**Question type:** Asks for two concrete deliverables: (i) variable importance, (ii) test error. Both must be present as explicit outputs.

**Expected deliverable:**
- A random forest fit on the training data with tuning.
- Variable importance output (plot or ranked table).
- Test set error on the held-out 20%.

**Accepted approaches:**
- `randomForest::randomForest()` with manual tuning, OR
- `caret::train(..., method = "rf" or "ranger")` with a tuning grid over `mtry` (and `min.node.size` / `splitrule` if using `ranger`), OR
- `ranger::ranger()` directly with CV-based tuning.

The `caret` + `ranger` route matches the course style most closely (follow L8/L9 caret pattern).

**Key steps the reviewer MUST check line by line:**
1. `set.seed` before fitting.
2. Training done only on training data.
3. `mtry` (and optionally other RF hyperparameters) is tuned — a grid is specified, not defaults-only. Tuning is by CV or OOB.
4. Variable importance is reported — either via `importance()`, `varImpPlot()`, `vip::vip()`, or `caret::varImp()`. A numeric ranking or a plot is acceptable.
5. Predictions are made on the test set and a test error metric (MSE or RMSE) is computed and reported with units / clear labeling.
6. A brief comment identifies the top few important variables.

**Common mistakes to watch for:**
- Reporting OOB error as "test error" — OOB is an internal estimate, NOT test error. This is the single most common error on this problem.
- Not tuning `mtry` at all.
- Computing test error on the training set by accident.
- Using classification importance measures without noticing this is a regression task (RF will default to `%IncMSE` / `IncNodePurity` for regression — verify the student's output is the regression variant).
- Reporting importance ONLY as a plot with no interpretation, or ONLY as a table with no plot — either alone is fine, but there should be at least one clearly presented form.

---

### 1(c) Boosting

**Question type:** Same as (b) — asks for variable importance and test error as explicit deliverables.

**Expected deliverable:**
- Gradient boosting model fit on training data with tuning.
- Variable importance.
- Test set error.

**Accepted approaches:**
- `gbm::gbm()` with manual CV / held-out monitoring, OR
- `caret::train(..., method = "gbm")` with a grid over `n.trees`, `interaction.depth`, `shrinkage`, `n.minobsinnode`, OR
- `xgboost` / `caret::train(method = "xgbTree")` with an appropriate tuning grid.

Any of the above is acceptable if CV-tuned properly.

**Key steps the reviewer MUST check line by line:**
1. `set.seed` before tuning.
2. Tuning grid is specified, covers at least two of: number of trees, tree depth / interaction depth, learning rate / shrinkage. A single-point grid is NOT acceptable as "tuning."
3. CV (or OOB for gbm) is used to pick the best hyperparameters.
4. Distribution / loss is appropriate for regression (`distribution = "gaussian"` for `gbm`, or `objective = "reg:squarederror"` for xgboost). Using `"bernoulli"` or a classification objective is wrong here.
5. Variable importance is reported (via `summary.gbm`, `vip::vip`, or `caret::varImp`).
6. Test error computed on held-out 20% and clearly labeled.
7. Light narrative comparing to the RF result from (b) is a plus, not required.

**Common mistakes to watch for:**
- Wrong `distribution` argument (e.g., `"bernoulli"` for a regression response).
- Not tuning — just using `gbm` defaults.
- Reporting training error or CV error as "test error."
- Overfitting: using an enormous `n.trees` without CV / early stopping.
- Forgetting to pass `n.trees = best_iter` to `predict.gbm` — a frequent silent bug.

---

## Problem 2 — Classification on the `auto.csv` data (response: `mpg_cat`)

70/30 split required. Verify.

### 2(a) Classification tree — CV optimal size vs 1-SE rule

**Question type:** The problem asks TWO explicit questions:
1. Which tree size corresponds to the lowest CV error?
2. Is that the same size as the one picked by the 1-SE rule?

Both questions require a direct, explicit answer in the narrative (a sentence stating the two sizes and whether they agree). A number-free answer is not acceptable.

**Expected deliverable:**
- Classification tree grown on training data.
- CV error as a function of tree size / complexity parameter.
- Explicit identification of (a) the minimum-CV size and (b) the 1-SE size.
- A sentence stating whether they agree.

**Accepted approaches:**
- `rpart::rpart(..., method = "class")` followed by `printcp()` and `plotcp()`. The minimum-CV `cp` is selected via `which.min(cptable[,"xerror"])`. The 1-SE rule selects the smallest tree whose `xerror` is within one SE of the minimum (`xerror + xstd` of the minimum row). Then use `prune()` with each selected `cp` to get the two tree sizes.
- `caret::train(..., method = "rpart")` with a `cp` grid and a CV `trainControl`, inspecting `model$results` to identify the best `cp` and the 1-SE `cp` (`selectionFunction = "oneSE"` in `trainControl`).
- `tree::tree()` + `cv.tree()` + `prune.misclass()` — this is a direct ISL-style approach that maps cleanly onto the lecture (L10 pp. 15–18).

**Key steps the reviewer MUST check line by line:**
1. Classification-mode fit (`method = "class"` for `rpart`, or factor response for `tree()`).
2. CV is actually performed (either via `rpart`'s built-in `cp` table / `cv.tree()`, or via `caret::train` with `trainControl(method = "cv", ...)`).
3. The CV error curve is shown — either a `plotcp` plot, a `cv.tree` plot, or a `caret::train` plot.
4. A specific number is given for the minimum-CV tree size (number of terminal nodes, NOT `cp` value alone).
5. A specific number is given for the 1-SE tree size.
6. An explicit verbal answer: "yes, they are the same" OR "no, they differ (min-CV picks X, 1-SE picks Y)."

**Common mistakes to watch for:**
- Reporting `cp` instead of tree size in terminal nodes — the question asks about tree size.
- Confusing `xerror` (CV error) with `rel error` (training error) in the `rpart` cptable.
- Applying the 1-SE rule incorrectly: the rule picks the smallest (simplest) tree with CV error within 1 SE of the minimum, not the largest.
- Silently picking 1-SE via `trainControl(selectionFunction = "oneSE")` without saying so, then never articulating the 1-SE tree size.
- Running `cv.tree` without setting the seed → nondeterministic answer (acceptable if the student still answers the question with whatever they got, but they should fix the seed).
- Failing to answer the comparison question in words.

---

### 2(b) Boosting for classification

**Question type:** Asks for two deliverables:
1. Variable importance.
2. Test data performance.

Both must be explicitly present.

**Expected deliverable:**
- Boosted classifier fit on training data with CV tuning.
- Variable importance (plot or ranked table).
- Test set performance, reported with a clearly labeled metric.

**Accepted approaches:**
- `gbm::gbm(..., distribution = "bernoulli")` (requires converting factor to 0/1), OR
- `caret::train(..., method = "gbm")` with binary factor response and `trainControl` using `classProbs = TRUE`, `summaryFunction = twoClassSummary`, `metric = "ROC"` (this matches L9 style), OR
- `xgboost` / `caret::train(method = "xgbTree")` with binary objective / `"logLoss"` or `"ROC"` metric.

**Key steps the reviewer MUST check line by line:**
1. `set.seed` before tuning.
2. Response is treated as a factor (classification), and the appropriate classification loss is used (`distribution = "bernoulli"` for `gbm`, NOT `"gaussian"`).
3. A real tuning grid is specified and CV-tuned.
4. Variable importance is reported (via `summary.gbm`, `vip::vip`, or `varImp`).
5. Test set predictions are made on the held-out 30%.
6. Test performance is reported. Acceptable forms (at least one, preferably more): test classification error / accuracy, confusion matrix (`caret::confusionMatrix` as in L8/L9), test AUC with an ROC curve (`pROC::roc` as in L8/L9).
7. The reported metric is computed on the TEST set, not CV.

**Common mistakes to watch for:**
- `distribution = "gaussian"` on a binary classification problem — a common silent bug.
- Reporting CV ROC as "test ROC."
- Reporting only accuracy without noting class balance — for `mpg_cat`, classes are roughly balanced, but the student should still be clear whether they are reporting accuracy, error, or AUC.
- Forgetting `type = "prob"` when predicting on the test set for AUC computation.
- For `gbm` with `predict`, forgetting `n.trees = best_iter` or `type = "response"`.
- Label encoding mismatch: confusing the positive class in `confusionMatrix` (verify `positive = ...` argument matches the student's reported direction).
- Not reporting variable importance at all and only giving test performance (or vice versa).

---

## Global deduction guidance

Reviewer should deduct points for:
- Missing reproducibility (no `set.seed`) — small deduction.
- Wrong split proportion or split done on the wrong dataset — moderate deduction.
- Reporting training/CV error as test error — moderate deduction (this is a conceptual error the course emphasizes).
- Missing required deliverable (e.g., no variable importance plot in 1(b) or 1(c)) — moderate deduction per missing item.
- Wrong loss / objective (regression vs classification mixup) — large deduction.
- No tuning where tuning is clearly expected (RF, boosting) — moderate deduction.
- Not answering the explicit comparison question in 2(a) in words — moderate deduction.

Reviewer should NOT deduct for:
- Using `caret` vs `tidymodels` vs direct package calls — all are acceptable as long as methodology is sound.
- Differing numerical answers that arise from seed / fold / grid differences, provided the methodology is correct.
- Choice of RMSE vs MSE for regression test error, provided the metric is clearly labeled.
- Minor stylistic choices in plots.

---

## Reviewer workflow checklist

For each submission, the TA should:

1. Open the student's Rmd and verify it knits cleanly (or verify the PDF matches the Rmd).
2. Check the data import and split code against the stated proportions.
3. For each sub-problem, locate the code chunk(s) and walk through the "Key steps MUST check" list above line by line.
4. Independently sanity-check numerical test errors against the student's own output (no answer key comparison needed; just verify that the reported number matches what their code produces and is on the correct dataset).
5. Verify every explicit question in the homework text is answered in narrative form, not only in code output.
6. Note any conceptual errors from the "Common mistakes" list.
