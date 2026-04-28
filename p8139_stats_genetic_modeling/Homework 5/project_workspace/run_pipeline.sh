#!/bin/bash
# Runs the silencing-focused TCGA-KIRC pipeline end-to-end. Exits on first
# failure. Designed to be re-runnable: each step writes its own log.
#
# Stages
#   01–04  acquisition + preprocessing (KIRC 450K + clinical)
#   05–06  paired DM and permutation null
#   07–09  silencing screen (MEA + sensitivity grid + pathway/immune dropout)
#   10     VHL subgroup case study
#   11–12  MethylMix benchmark
#   13–15  TCGA-KIRP cross-cohort replication
#   20–21  figures for the final report

set -euo pipefail
cd "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"

run_step () {
  local name="$1"; shift
  local script="$1"; shift
  local log="$1"; shift
  echo "=== $name ==="
  if ! Rscript "$script" > "$log" 2>&1; then
    echo "!!! STEP FAILED: $name  (see $log for full output)"
    tail -40 "$log"
    exit 1
  fi
  tail -20 "$log"
}

run_step "01 download KIRC"           01_download_kirc.R         01_download_kirc.log
run_step "02 EDA KIRC"                02_eda_kirc.R              02_eda_kirc.log
run_step "03 PC diagnostic"           03_pc_diagnostic.R         03_pc_diagnostic.log
run_step "04 preprocess"              04_preprocess.R            04_preprocess.log
run_step "05 Q1 paired DM"            05_q1_paired_dm.R          05_q1_paired_dm.log
run_step "06 permutation null"        06_permutation_null.R      06_permutation_null.log
run_step "07 methylation silencing"   07_methylation_silencing.R 07_methylation_silencing.log
run_step "08 silencing FDR grid"      08_silencing_fdr_grid.R    08_silencing_fdr_grid.log
run_step "09 silencing pathway"       09_silencing_pathway.R     09_silencing_pathway.log
run_step "10 VHL subgroup"            10_vhl_subgroup.R          10_vhl_subgroup.log
run_step "11 MethylMix benchmark"     11_methylmix_benchmark.R   11_methylmix_benchmark.log
run_step "12 MethylMix compare"       12_methylmix_compare.R     12_methylmix_compare.log
run_step "13 KIRP download"           13_kirp_download.R         13_kirp_download.log
run_step "14 KIRP pipeline"           14_kirp_pipeline.R         14_kirp_pipeline.log
run_step "15 KIRP sensitivity"        15_kirp_sensitivity.R      15_kirp_sensitivity.log
run_step "20 figure main"             20_figure_main.R           20_figure_main.log
run_step "21 figure DM diag"          21_figure_dm_diag.R        21_figure_dm_diag.log

echo "=== Pipeline complete ==="
