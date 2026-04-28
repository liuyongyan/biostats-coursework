#!/bin/bash
# Wait for KIRP download (27) to finish, then run pipeline (28),
# write a DONE marker. Robust: polls every 60s, doesn't depend on a specific PID.

WORK="/Users/yliu/Desktop/Columbia - Biostatistics/_homeworks/p8139_stats_genetic_modeling/Homework 5/project_workspace"
LOG="$WORK/kirp_chain.log"
cd "$WORK"

{
  echo "[$(date)] chain started; waiting for 27_download_kirp.R to exit..."
} >> "$LOG"

# Wait for 27 to finish
while pgrep -f '27_download_kirp' > /dev/null; do
  sleep 60
done

echo "[$(date)] download script no longer running" >> "$LOG"

# Verify expected outputs exist before running pipeline
NEEDED=( "data_kirp/kirp_450k_SE.rds" "data_kirp/kirp_rnaseq_SE.rds" "data_kirp/clinical.rds" )
MISSING=0
for f in "${NEEDED[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "[$(date)] MISSING: $f" >> "$LOG"
    MISSING=1
  else
    echo "[$(date)] OK: $f ($(du -h "$f" | cut -f1))" >> "$LOG"
  fi
done

if [[ "$MISSING" -ne 0 ]]; then
  echo "[$(date)] download incomplete — pipeline NOT run. Investigate kirp_download.log." >> "$LOG"
  touch KIRP_DOWNLOAD_FAILED
  exit 1
fi

echo "[$(date)] all KIRP data present, launching 28_kirp_pipeline.R..." >> "$LOG"
Rscript 28_kirp_pipeline.R > kirp_pipeline.log 2>&1
RC=$?
echo "[$(date)] pipeline exit code: $RC" >> "$LOG"

if [[ "$RC" -eq 0 ]]; then
  touch KIRP_DONE
  echo "[$(date)] KIRP_DONE marker written." >> "$LOG"
else
  touch KIRP_PIPELINE_FAILED
  echo "[$(date)] KIRP_PIPELINE_FAILED marker written." >> "$LOG"
fi
