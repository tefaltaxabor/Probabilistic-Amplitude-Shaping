#!/bin/bash
cd ~/Probabilistic-Amplitude-Shaping || exit 1
mkdir -p logs results
LOG="logs/showcase_$(date +%Y%m%d_%H%M%S).log"
echo "=== inicio $(date) en $(hostname) ===" | tee "$LOG"
matlab -nodisplay -nosplash -batch showcase_pas 2>&1 | tee -a "$LOG"
echo "=== fin $(date) (exit=$?) ===" | tee -a "$LOG"
