#!/bin/bash
# Run every simulation whose figures go into results/, on an LNT machine.
# Strictly sequential: one MATLAB at a time, no parfor anywhere.
#
#   nohup ./run_all.sh > /dev/null 2>&1 &      # survives logout
#   JOBS="compare_fec" ./run_all.sh             # run a subset
#
# One log per script in logs/, same format as run_pas_showcase.sh.
cd ~/Probabilistic-Amplitude-Shaping || exit 1
mkdir -p logs results

JOBS="${JOBS:-showcase_pas compare_shaping compare_fec compare_blocklength}"

for job in $JOBS; do
    LOG="logs/${job}_$(date +%Y%m%d_%H%M%S).log"
    echo "=== inicio $(date) en $(hostname) ===" | tee "$LOG"
    matlab -nodisplay -nosplash -batch "$job" 2>&1 | tee -a "$LOG"
    rc=${PIPESTATUS[0]}
    echo "=== fin $(date) (exit=$rc) ===" | tee -a "$LOG"
done
