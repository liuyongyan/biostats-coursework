# P8106 Data Science II — Homework 5 Grading Rubric

> **Reviewer instructions.** This rubric is **principle-based**: it specifies what the student's submission must contain and how it must be checked, not the numerical answers. Do **not** mark a student wrong because their numbers differ from a TA-computed reference; mark based on whether the methodology is sound, the questions are answered explicitly, and the code style follows the lecture R scripts. When checking code consistency, open the relevant lecture `.Rmd` file side-by-side and compare the student's pipeline to the canonical lecture pipeline (function names, tuning-grid construction, seeding, prediction, error reporting).
>
> Lecture R reference files:
> - Lecture 11 (ensembles, esp. boosting via `gbm` / `caret` `method="gbm"`): `Data Science II/Lecture 11/L11_2026.Rmd`
> - Lecture 13 (SVM/SVC via `e1071::tune.svm` and `caret` `method="svmLinear"|"svmLinear2"|"svmRadial"|"svmRadialSigma"|"svmRadialCost"`): `Data Science II/Lecture 13/L13_2026.Rmd`
> - Lecture 14 (hierarchical clustering with `hclust`, `dist`, `cutree`, `fviz_dend`, scaling): `Data Science II/Lecture 14/L14_2026.Rmd`

---

## Global expectations (apply to the whole submission)

### Key things that MUST be present
- A reproducible R / R Markdown document (or `.R` script + knitted PDF/HTML) with all code visible enough to be reviewed.
- A single, explicitly chosen random seed set **before** the train/test split, and seeds set before each tuning step that involves resampling (mirroring how every chunk in the lecture `.Rmd`s starts with `set.seed(...)`).
- For Problem 1: an explicit 70 / 30 train–test split using a documented method (`rsample::initial_split(..., prop = 0.7)`, `caret::createDataPartition`, or equivalent). The same split must be reused across (a), (b), and (c).
- The response variable `mpg_cat` must be a factor with explicit, consistent level ordering before any model is fitted (the lecture code does `factor(diabetes, c("pos","neg"))` precisely to fix the positive class for `confusionMatrix`).
- Predictors used must be exactly the seven specified in the prompt: `cylinders`, `displacement`, `horsepower`, `weight`, `acceleration`, `year`, `origin`. Penalize use of `name` or any other extra column; penalize dropping any of the seven. `origin` must be handled as either a factor or numeric consistently — flag inconsistencies.
- Numerical results that the question asks for (error rates, cluster memberships, scaling-effect verdict) must be reported **in prose or a table inside the write-up**, not only buried in raw R output.

### Method-consistency checks
- Library calls match the lecture stack: `caret`, `tidymodels`/`rsample`, `e1071`, `kernlab`, `gbm`, `ranger`/`randomForest`, `factoextra`, `ISLR`. Heavy use of an off-syllabus framework (e.g., `mlr3`, `tidymodels` with custom recipes the lectures didn't show) is acceptable only if the submission clearly mirrors the same workflow.
- Tuning is done by cross-validation on **training data only**; test data must never enter a tuning grid or `tune.svm` call.

### Common mistakes to flag
- Setting the seed only once at the top, after the data has already been read and split via something stochastic.
- Refitting a fresh split for each subpart, so (a)/(b)/(c) are not comparable.
- Using `accuracy` interchangeably with `error rate` without converting (`error = 1 - accuracy`); the question asks for **error rates**.
- Reporting only the test error and silently omitting the training error (or vice-versa) in 1(a) and 1(b).
- Hard-coded numerical answers in the prose that disagree with the displayed R output.

---

## Problem 1 — Support Vector Machines on `auto.csv`

The student must use the same `auto.csv` data that was used in HW3, restrict to the seven listed predictors, and use a single train/test split (70/30) shared across (a), (b), (c).

### 1(a) Support Vector Classifier (linear kernel)

**Direct question to be answered:** *"What are the training and test error rates?"* — **both** numbers must appear explicitly.

#### Key things that MUST be present
- A linear-kernel SVM is fit on the training set with `cost` tuned over a grid via cross-validation. Acceptable canonical pipelines (Lecture 13):
  - `e1071::tune.svm(..., kernel = "linear", cost = <grid>)` followed by `best.model` extraction, **or**
  - `caret::train(..., method = "svmLinear", tuneGrid = data.frame(C = <grid>), trControl = trainControl(method = "cv"))`, **or**
  - `caret::train(..., method = "svmLinear2", tuneGrid = data.frame(cost = <grid>), ...)`, **or**
  - `tidymodels` with `svm_linear(cost = tune())` and `tune_grid()`.
- Some plot or table showing the CV tuning curve over `cost` (lecture uses `plot(linear.tune)` or `plot(svml.fit, highlight = TRUE, xTrans = log)`).
- The selected `cost` (best parameters) must be reported.
- **Both** training and test misclassification error rates must be computed and reported as numbers in the write-up. Acceptable computation paths: predict on training/test data and compute `mean(pred != truth)`, or compute `1 - confusionMatrix(...)$overall["Accuracy"]`.

#### Method-consistency checks
- Tuning grid for `cost` should span several orders of magnitude on a log scale (the lecture uses `exp(seq(-5, 2, len = 50))` for the linear kernel). A single value or a very narrow grid should be flagged.
- Scaling: `tune.svm`/`svm` defaults to `scale = TRUE`; `caret`'s svm methods scale internally. If the student fits via a non-scaling pathway, they should justify it.
- Compare the student's chunk against the **"Linear boundary"** section of `L13_2026.Rmd` line-by-line for argument style.

#### Common mistakes to flag
- Reporting only test error, or only the confusion matrix without converting to an error rate.
- Tuning `cost` on the test data, or fitting one model with a fixed `cost` and never tuning at all.
- Treating `mpg_cat` as numeric so the model fits a regression silently.
- Using `kernel = "radial"` here by mistake.

---

### 1(b) Support Vector Machine with radial kernel

**Direct question to be answered:** *"What are the training and test error rates?"* — both numbers explicitly.

#### Key things that MUST be present
- A radial-kernel SVM fit on the training set with **both** tuning parameters tuned: `cost` **and** `gamma` (or, equivalently in `caret`, `C` and `sigma`). Canonical pipelines (Lecture 13):
  - `e1071::tune.svm(..., kernel = "radial", cost = <grid>, gamma = <grid>)`, **or**
  - `caret::train(..., method = "svmRadialSigma", tuneGrid = expand.grid(C = ..., sigma = ...))`, **or**
  - `caret::train(..., method = "svmRadialCost", ...)` (single `sigma` from `kernlab::sigest`), **or**
  - `tidymodels` `svm_rbf(cost = tune(), rbf_sigma = tune())`.
- Best tuning parameters are reported.
- Training **and** test error rates are reported numerically.

#### Method-consistency checks
- Both `cost` and the kernel-width parameter must have been varied on a non-trivial grid (the lecture uses `cost = exp(seq(1, 7, len = 50))`, `gamma = exp(seq(-10, -2, len = 20))`). A single fixed `gamma` is acceptable only if the student explicitly chose `svmRadialCost` (which uses `kernlab::sigest`); otherwise flag.
- Compare to the **"Radial kernel"** section of `L13_2026.Rmd`.

#### Common mistakes to flag
- Tuning only `cost` while leaving `gamma`/`sigma` at the default (this defeats the purpose of the radial kernel).
- Choosing a tuning grid that is clearly bottoming/topping out at a corner of the grid without commentary or grid expansion.
- Forgetting to predict and report **training** error in addition to test error.
- Re-splitting the data here so the train/test set differs from 1(a).

---

### 1(c) Comparison with boosting

**Direct question to be answered:** *"Do SVC or SVM provide better predictive performance on this dataset?"* — the write-up must contain an explicit verdict (e.g., "SVM-radial outperforms SVC and boosting on test error" or "boosting performs best, both SVMs are worse"). A plain table without a written conclusion is not enough.

#### Key things that MUST be present
- A boosting model fit on the **same training set** as 1(a)/(b). Canonical pipelines (Lecture 11, classification section):
  - `gbm(..., distribution = "adaboost" | "bernoulli", n.trees, interaction.depth, shrinkage, cv.folds)` with `gbm.perf(..., method = "cv")` to pick `n.trees`, **or**
  - `caret::train(..., method = "gbm", tuneGrid = expand.grid(n.trees, interaction.depth, shrinkage, n.minobsinnode), trControl, distribution = "adaboost"|..., metric = "ROC"|"Accuracy", verbose = FALSE)`.
- A test error rate (or equivalently AUC / accuracy) for the boosting model, **computed on the same test set** used for SVC and SVM.
- A side-by-side comparison (table or text) of SVC vs SVM-radial vs boosting on the test set, and an explicit answer to the question of which performs better.

#### Method-consistency checks
- For `gbm` the response must be encoded correctly (for `distribution = "adaboost"` the lecture converts the factor to numeric 0/1; for `distribution = "bernoulli"` similarly). Mismatched response encoding is a real bug.
- Tuning grid for boosting should vary at least `n.trees`, `interaction.depth`, and `shrinkage` (compare `gbmA.grid` in `L11_2026.Rmd`).
- The same `set.seed` discipline should be applied as in earlier subparts.

#### Common mistakes to flag
- Comparing on the wrong metric (e.g., training accuracy for SVMs vs CV accuracy for boosting). The comparison must be on the **test set**.
- Comparing SVMs and boosting that were trained on different splits.
- Not stating a verdict in plain English.
- Using boosting's CV result as the comparison number while using SVM's test result as the comparison number for SVMs (apples-to-oranges).

---

## Problem 2 — Hierarchical clustering on `USArrests`

Data: `ISLR::USArrests` (or `datasets::USArrests`), all 50 states, four features (`Murder`, `Assault`, `UrbanPop`, `Rape`).

### 2(a) Hierarchical clustering, complete linkage, Euclidean distance, no scaling

**Direct question to be answered:** *"Which states belong to which clusters?"* — the write-up must explicitly enumerate the membership of each of the three clusters (e.g., a table or three labelled lists of state names).

#### Key things that MUST be present
- Clustering is computed on the **unscaled** data: `hclust(dist(USArrests), method = "complete")`. Compare to the `hc.complete` line of `L14_2026.Rmd`.
- A dendrogram is plotted (acceptable: `plot(hc)`, `factoextra::fviz_dend(hc, k = 3, ...)`, or any equivalent). The cut into three clusters should be visible (e.g., `k = 3` rectangles or a horizontal cut line).
- The dendrogram is cut into three groups, e.g., `cutree(hc, k = 3)` or `cutree(hc, h = <height>)`. Either form is acceptable so long as exactly three clusters result.
- The membership of each cluster is reported in the write-up as state names — not just a vector of integers. A small table mapping state to cluster number (or three lists "Cluster 1: …, Cluster 2: …, Cluster 3: …") is the expected presentation.

#### Method-consistency checks
- Distance is the default Euclidean (no `method` argument to `dist`, or explicit `method = "euclidean"`).
- Linkage is `method = "complete"` in `hclust`.
- Row names of the input are preserved so that `cutree` returns named state assignments (the lecture's Pokemon example uses `rownames(dat1) <- dat[,1]` for exactly this reason).

#### Common mistakes to flag
- Cutting at `k = 4` instead of `k = 3` (the lecture example uses `k = 4` — students often copy that without adjusting).
- Reporting cluster sizes (e.g., "17, 14, 19") without listing which states are in each cluster.
- Scaling the data here in 2(a) by accident — this part must be on **unscaled** data so that 2(b) can compare against it.
- Using `kmeans` instead of `hclust`, or using a non-complete linkage.

---

### 2(b) Hierarchical clustering after scaling to unit standard deviation

**Direct questions to be answered (BOTH must be addressed):**
1. *"Does scaling the variables change the clustering results?"* — explicit yes/no answer with supporting evidence (e.g., comparison of cluster memberships across 2(a) and 2(b), a contingency table of `table(cluster_a, cluster_b)`, or a side-by-side dendrogram).
2. *"In your opinion, should the variables be scaled before the inter-observation dissimilarities are computed?"* — explicit opinion, with at least one sentence of justification grounded in the data (e.g., the variables are on different scales: `Assault` is in the hundreds while `Murder` and `Rape` are in single/double digits and `UrbanPop` is a percentage, so unscaled Euclidean distance is dominated by `Assault`).

#### Key things that MUST be present
- Variables are scaled to unit standard deviation before clustering: `scale(USArrests)` (centered & scaled by default) — the lecture uses exactly `scale(...)` on the Pokemon stats. Scaling only by the standard deviation without centering is also defensible since the question only requires unit SD; do not penalize either choice.
- Same linkage (`complete`) and same distance (Euclidean) as in 2(a) — the only thing that changes is the scaling.
- A second dendrogram (or `cutree` output) for the scaled clustering.
- A direct comparison of memberships between 2(a) and 2(b) (e.g., `table(cutree_a, cutree_b)`, or a list of states whose cluster changed).
- An explicit yes/no on whether scaling changed results.
- An explicit opinion on whether to scale, with reasoning.

#### Method-consistency checks
- The scaling must happen **before** `dist()`. Calling `dist()` first and then scaling is wrong.
- The student should not change the linkage or distance metric between 2(a) and 2(b); doing so confounds the scaling effect.
- Compare to the `dat1 <- scale(...)` then `hclust(dist(dat1), method = "complete")` pattern in `L14_2026.Rmd`.

#### Common mistakes to flag
- Vague answers like "the results look different" without showing how.
- Refusing to give an opinion on whether to scale (the question explicitly asks for one).
- Justifying "should scale" purely on principle without referencing the actual variable scales in `USArrests` (the strongest answers note that `Assault` ranges into the hundreds while the other three are on much smaller scales, so unscaled Euclidean distance is essentially a `Assault`-only distance).
- Cutting at a different `k` between 2(a) and 2(b), making the comparison meaningless.
- Re-running on a different subset of variables.

---

## Final write-up checks

- Every question mark in the homework PDF has a corresponding direct sentence in the student's write-up. Specifically:
  - 1(a): "training error = …, test error = …"
  - 1(b): "training error = …, test error = …"
  - 1(c): explicit verdict on SVC vs SVM vs boosting.
  - 2(a): explicit listing of which states are in each of the three clusters.
  - 2(b): explicit yes/no on whether scaling changes results, **and** explicit opinion on whether to scale.
- Code is reproducible from the seed onward.
- Output figures (tuning curves, dendrograms) are present and legible.
- The submission is internally consistent: numbers in the prose match numbers in the displayed output.
