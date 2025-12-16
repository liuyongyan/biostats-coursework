# biostats-coursework

Homework assignments and cheat sheets from my M.S. in Biostatistics program (2025–2027) at Columbia University Mailman School of Public Health.

## Courses

| Course Code | Course Name | Contents |
|-------------|-------------|----------|
| P8104 | Probability | Homework 1–11, Cheatsheets (Midterm, Final) |
| P8105 | Data Science I | Homework 1–3, 5–6, Midterm Project, Final Project |
| P8130 | Biostatistical Methods I | Homework 1–5, Cheatsheets (Midterm, Final) |
| P8149 | Human Population Genetics | Homework 1–4, Cheatsheets (Midterm, Final) |

## Directory Structure

```
├── p8104_probability/
├── p8105_data_science/
├── p8130_biostats_methods/
└── p8149_population_genetics/
```

## Naming Conventions

### Directories
- Course folders: `p[code]_[name]` (lowercase, underscores)
- Subfolders: `Homework [N]/` or `Cheatsheet - [Exam]/`

### Files
| Type | Pattern | Example |
|------|---------|---------|
| Assignment | `P[code]_HW[N].pdf` | `P8104_HW1.pdf` |
| Submission | `P[code]_HW[N]_yl6107.pdf/.Rmd` | `P8104_HW1_yl6107.Rmd` |
| Cheatsheet | `P[code]_CheatSheet_[Exam].pdf/.Rmd` | `P8104_CheatSheet_Final.pdf` |
| Solution | `P[code]_HW[N]_sol.pdf` | `P8149_HW1_sol.pdf` |

**Note:** P8105 Data Science I uses lowercase naming (`p8105_hw1_yl6107.Rmd`) following the course's GitHub-based workflow.

## Environment Setup

To knit the R Markdown files, the following software and packages are required.

### Software
- **R** (tested with R 4.5.1)
- **TinyTeX** or a full TeX distribution (for PDF output)

### R Packages

Install all required packages:

```r
install.packages(c(
  "rmarkdown", "knitr", "tidyverse", "broom", "dplyr",
  "ggplot2", "ggrepel", "glmnet", "janitor", "maps",
  "MASS", "modelr", "patchwork", "readxl", "scales"
))

# For P8105 Data Science assignments
install.packages("devtools")
devtools::install_github("p8105/p8105.datasets")

# For PDF output
install.packages("tinytex")
tinytex::install_tinytex()
```

### Output Formats
- **P8104, P8130, P8149**: PDF documents (requires LaTeX)
- **P8105**: GitHub documents (Markdown)
