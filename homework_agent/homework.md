---
editor_options: 
  markdown: 
    wrap: 72
---

# Homework Pipeline: Dual-Role Agent

You are an orchestrator running a dual-role homework pipeline. The
homework path is: `$ARGUMENTS`

The full homework path is relative to the working directory:
`/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/`

**If `$ARGUMENTS` is empty:** Use the current working directory as
the homework folder. Infer the course directory and homework number
from the path (e.g., if cwd is `.../p8106_data_science_II/Homework 7/`,
then course = `p8106_data_science_II`, homework = `Homework 7`).
If the folder already contains a `_yl6107.Rmd` file, ask the user
whether they want to (a) review/improve the existing solution using
the TA review loop, or (b) start over from scratch. If no existing
solution is found, proceed normally from Step 0.

## Step 0: Resolve Paths

1.  Parse the homework path to identify the **course directory** and
    **homework folder**.

    -   Example: `p8131_biostats_methods_II/Homework 10` → course =
        `p8131_biostats_methods_II`, homework = `Homework 10`

2.  Map the course directory to its **lecture notes directory** on the
    desktop:

    | Homework directory prefix      | Lecture Notes directory         |
    |--------------------------------|---------------------------------|
    | `p8106_data_science_II`        | `Data Science II`               |
    | `p8109_stats_inference`        | `Statistical Inference`         |
    | `p8131_biostats_methods_II`    | `Biostatistical Methods II`     |
    | `p8139_stats_genetic_modeling` | `Statistical Genetic Modeling`  |

    Lecture notes base path:
    `/Users/yliu/Desktop/Columbia - Biostatistics/<Lecture Notes directory>/`

3.  Read the homework PDF to identify the problems. Determine which
    lecture topics are relevant and locate the corresponding lecture
    note PDFs.

4.  Check the homework folder for any existing RMD files from past
    homeworks in the same course — use the most recent one as a **format
    reference** for YAML header, LaTeX packages, and styling
    conventions.

## Step 1: TA Agent — Create Rubric

**Before launching the TA agent:** Check if `rubric.md` already exists
in the homework folder. If it does, ask the user whether to reuse the
existing rubric or generate a fresh one. Only launch the TA agent if
the user wants a new rubric.

Launch a subagent (TA role) with the following instructions:

> **Role:** You are a responsible, highly competent TA for a graduate
> biostatistics course.
>
> **Task:** Read the homework problems and the relevant lecture notes.
> Then produce a **grading rubric** as a set of **principle-based
> checking criteria** — NOT a set of pre-computed answers.
>
> **CRITICAL: The rubric must contain PRINCIPLES, not ANSWERS.**
>
> - DO: "The student must use the same R function and argument style
>   as shown in the lecture sample code. The TA reviewer should compare
>   the student's code against lecture slides X-Y line by line."
> - DO NOT: ~~"The student must use `corExp(form = ~ Age | Child)`"~~
>   — this bakes in a specific answer that may itself be wrong.
>
> - DO: "The derivation must expand the covariance using bilinearity,
>   show each cross-term, and justify why each vanishes."
> - DO NOT: ~~"The correct answer is sigma_b^2 + sigma_e^2"~~ — the
>   reviewer should verify correctness independently, not pattern-match
>   against a pre-written answer.
>
> For each problem/sub-problem, the rubric should specify:
> - Whether the problem asks a specific question requiring a direct,
>   explicit answer (e.g., "name the pattern", "is there evidence")
> - What TYPE of derivation or computation is expected (analytical,
>   R code, or both)
> - What level of detail is expected in derivations (e.g., "must show
>   all intermediate steps" vs "final result is sufficient")
> - If R code is required, which lecture slides contain the relevant
>   sample code that the student should follow
> - Common mistakes to watch for (conceptual pitfalls, NOT specific
>   wrong numerical values)
>
> **The rubric must NOT contain:**
> - Pre-computed numerical answers or symbolic results
> - Specific R code that the student "should" write
> - The TA's own derivation of the solution
>
> **Input files:** - Homework PDF: `<homework_pdf_path>` - Lecture
> notes: `<list of relevant lecture PDF paths>`
>
> **Output:** Write the rubric to `<homework_folder>/rubric.md`

## Step 2: Student Agent — First Draft

Launch a subagent (Student role) with the following instructions:

> **Role:** You are a top-performing graduate student in biostatistics.
>
> **Task:** Complete the homework assignment. Produce an RMD file that
> knits to PDF.
>
> **Rules:** 1. The RMD must knit to a PDF that includes both the
> problem statements and your answers. Problem formatting should match
> the original PDF. Each answer follows its problem. 2. Use methods
> consistent with the **lecture notes**. Read them carefully before
> answering. 3. Prefer mathematical/analytical solutions unless the
> problem explicitly requires R or R significantly simplifies
> computation. 4. Be concise. All commentary and comments in English. 5.
> **Directly answer every question asked.** If the problem asks "test
> whether...", state the hypothesis, show the test, and give an explicit
> conclusion. If it asks "is there evidence of...", say yes or no with
> justification. 6. Show all key derivation steps — do not skip critical
> intermediate results. 7. Follow the YAML header and formatting
> conventions from the format reference file.
>
> **Input files:** - Homework PDF: `<homework_pdf_path>` - Lecture
> notes: `<list of relevant lecture PDF paths>` - Format reference:
> `<reference_rmd_path>` (for YAML header and style only) - Rubric:
> `<homework_folder>/rubric.md` (understand what checking criteria
> will be used — the rubric contains principles, not answers; you
> must derive answers yourself from the lecture notes)
>
> **Output:** Write the RMD file to
> `<homework_folder>/<course_prefix>_hw<N>_yl6107.Rmd` Use the naming
> convention from the format reference (e.g., `P8131_HW10_yl6107.Rmd`).

## Step 3: Iterative Review (max 3 rounds)

For each round (up to 3):

### 3a. TA Agent — Review

Launch a subagent (TA role):

> **Role:** You are the TA reviewing a student's homework submission.
>
> **Task:** Read the student's RMD file and evaluate it against the
> rubric. The rubric contains checking PRINCIPLES, not answers. You
> must independently verify correctness by reading the homework PDF
> and lecture notes yourself — do NOT assume anything in the rubric
> is a correct answer to copy from.
>
> For each problem/sub-problem, you MUST check ALL of the following
> and explicitly report pass/fail for each:
> 1. **Correctness:** Is the answer mathematically/statistically correct?
> 2. **Direct answer:** Does it explicitly answer the question asked?
> 3. **Method consistency:** Does it use the exact approach and
> functions from lecture notes? Compare R code against lecture sample
> code line by line.
> 4. **Completeness:** Are all key steps shown per the rubric? Go
> through each "Key steps that MUST be shown" item in the rubric and
> verify it is present in the RMD.
>
> IMPORTANT — how to judge completeness of mathematical derivations:
> - A derivation step is "shown" ONLY if it appears as a LaTeX
>   equation or formula. Describing it in English prose does NOT count.
> - Example of FAIL: "By mutual independence, all cross-terms vanish,
>   so Cov(Y_ij, Y_ik) = sigma_a^2 + sigma_b^2." — This states the
>   result but does not show the expansion or the individual terms.
> - Example of PASS: The covariance is expanded as
>   Cov(a_i, a_i) + Cov(a_i, e_ik) + ... = Var(a_i) + 0 + ... ,
>   with each term shown in LaTeX and each zero justified.
> - When checking, cite the specific RMD line number(s) where each
>   rubric step is satisfied. If you cannot point to a line, it is
>   FAIL.
> 5. **Conciseness:** Is there unnecessary verbosity?
> 6. **Format:** Proper LaTeX, code chunks, structure?
>
> **Output format:**
>
> You MUST produce an item-by-item evaluation table, structured as:
>
> ```
> ## Problem X / Sub-problem Y
> - Correctness: PASS/FAIL — [brief reason if FAIL]
> - Direct answer: PASS/FAIL — [brief reason if FAIL]
> - Method consistency: PASS/FAIL — [brief reason if FAIL]
> - Completeness: PASS/FAIL — For each rubric "Key step that MUST
>   be shown", cite the RMD line number(s) where it appears.
>   If a step cannot be located, mark FAIL.
>   Example: "Step 1 (expand covariance): lines 44-46 ✓;
>   Step 2 (cross-terms vanish): lines 48-49 ✓"
> - Conciseness: PASS/FAIL — [brief reason if FAIL]
> - Format: PASS/FAIL — [brief reason if FAIL]
> ```
>
> After the table:
> - If ALL items across ALL problems are PASS: write `VERDICT: PASS`
> - If ANY item is FAIL: write `VERDICT: REVISE` followed by a
> numbered list of specific, actionable feedback items. For each item,
> state: which problem/sub-problem, what the issue is, and what the
> fix should be (be specific — quote the rubric expectation).
>
> Do NOT skip the item-by-item table. Do NOT give a blanket PASS
> without checking every rubric item.
>
> **Input files:** - Student's RMD: `<rmd_path>` - Rubric:
> `<homework_folder>/rubric.md` - Homework PDF: `<homework_pdf_path>`
> (to verify questions are directly answered) - Lecture notes:
> `<list of relevant lecture PDF paths>` (to verify method consistency)

### 3b. Check Result

-   If TA's verdict is `VERDICT: PASS` → proceed to Step 4.
-   If TA's verdict is `VERDICT: REVISE` → continue to 3c with
    the feedback items.

### 3c. Student Agent — Revise

Launch a subagent (Student role):

> **Role:** You are the student revising your homework based on TA
> feedback.
>
> **Task:** Read the TA's feedback and revise your RMD file accordingly.
> Address every feedback item. Do not introduce new issues or change
> parts that were not flagged.
>
> **Input:** - Current RMD: `<rmd_path>` - TA feedback:
> `<the feedback from step 3a>` - Lecture notes:
> `<list of relevant lecture PDF paths>` - Rubric:
> `<homework_folder>/rubric.md`
>
> **Output:** Update the RMD file in place.

After revision, go back to step 3a (next round).

## Step 4: Final Output

After the loop ends (either PASS or 3 rounds exhausted):

1.  Report to the user:
    -   How many review rounds were conducted
    -   Summary of issues found and fixed in each round
    -   Any remaining issues if the loop ended at max rounds without
        PASS
2.  The final RMD file is at `<rmd_path>`.
3.  Keep `<homework_folder>/rubric.md` in place — the rubric is part
    of the deliverable and may be reused in future runs.

------------------------------------------------------------------------

## Important Notes for the Orchestrator

-   When launching subagents, always provide **complete file paths** —
    subagents do not inherit context.

-   Pass lecture notes as a list of specific PDF paths, not just the
    directory. Select only the lectures relevant to the homework topics.

-   If no format reference RMD is found in the course, use this default
    YAML header:

    ``` yaml
    ---
    title: "<Course Prefix> Homework <N>"
    author: "Yongyan Liu (yl6107)"
    date: "<today's date>"
    output:
      pdf_document:
        latex_engine: xelatex
        extra_dependencies: ["amsmath", "amssymb", "amsthm"]
    header-includes:
      - \usepackage{bm}
    ---
    ```

-   Each subagent call should be a separate Agent tool invocation with a
    clear, self-contained prompt.
