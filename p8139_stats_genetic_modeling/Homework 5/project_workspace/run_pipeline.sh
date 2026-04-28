#!/bin/bash
# Runs the Step 3 -> Step 4 pipeline in order. Exits non-zero on first failure.
set -euo pipefail
cd "/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"

run_step () {
  local name="$1"; shift
  local script="$1"; shift
  local log="$1"; shift
  echo "=== $name ==="
  # Run Rscript; redirect to log; fail the whole pipeline on non-zero.
  if ! Rscript "$script" > "$log" 2>&1; then
    echo "!!! STEP FAILED: $name  (see $log for full output)"
    tail -40 "$log"
    exit 1
  fi
  # Show brief summary
  tail -20 "$log"
}

run_step "04 preprocess"          04_preprocess.R            04_preprocess.log
run_step "05 Q1 paired"           05_q1_paired.R             05_q1_paired.log
run_step "06 Q3 grade"            06_q3_tumor_grade.R        06_q3_tumor_grade.log
run_step "07 Q5 enrichment"       07_q5_enrichment.R         07_q5_enrichment.log
run_step "08 Q9 positive control" 08_q9_positive_control.R   08_q9_positive_control.log

echo "=== 99 audit ==="
Rscript 99_audit.R > 99_audit.log 2>&1 || echo "(audit returned non-zero — FAILures present; see 99_audit.log)"
tail -60 99_audit.log

echo "=== Pipeline complete ==="
