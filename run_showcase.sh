#!/bin/bash
# Run the PAS showcase on the LNT/ICE cluster.
#
#   Smoke test (seconds, no Monte Carlo):
#       ./run_showcase.sh 1,2
#   Full run (submit to SLURM if available, otherwise run locally):
#       sbatch run_showcase.sh
#
# The panel list can be overridden with the first argument or $PARTS.

#SBATCH --job-name=pas_showcase
#SBATCH --output=logs/pas_showcase_%j.out
#SBATCH --error=logs/pas_showcase_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results

PARTS="${1:-${PARTS:-1:6}}"

# Load MATLAB if the cluster uses environment modules.
if command -v module >/dev/null 2>&1; then
    module load matlab 2>/dev/null || echo "warning: 'module load matlab' failed, trying PATH"
fi

if ! command -v matlab >/dev/null 2>&1; then
    echo "ERROR: matlab not on PATH. Check 'module avail matlab'." >&2
    exit 1
fi

echo "=== $(date) | host $(hostname) | parts [$PARTS] ==="
matlab -nodisplay -nosplash -batch "showcase_pas([$PARTS])"
echo "=== $(date) | finished ==="
