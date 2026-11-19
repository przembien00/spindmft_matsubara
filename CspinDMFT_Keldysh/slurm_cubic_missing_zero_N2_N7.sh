#!/usr/bin/env bash
# Missing zero-field sweep jobs, excluding N=7 at beta >= 4.
# Submit from CspinDMFT_Keldysh after creating logs/:
#   mkdir -p logs
#   sbatch slurm_cubic_missing_zero_N2_N7.sh
# N=2 zero-field runs use 500000 real-time samples per core in the dedicated plan.
# Preparation samples and all N=7 sample counts remain unchanged.
# IDs retain their meaning in slurm_cubic_beta_sweep.sh and the sample plan:
#   N=2 FM:  beta 1.5, 3, 3.5, 4.5       -> 72,144,168,216
#   N=2 AFM: beta 1, 1.5, 3, 3.5, 4.5, 5 -> 51,75,147,171,219,243
#   N=7 FM:  beta 3, 3.5                  -> 156,180
#   N=7 AFM: beta 2.5, 3, 3.5             -> 135,159,183
#SBATCH --job-name=cluster_keldysh_missing_zero
#SBATCH --array=51,72,75,135,144,147,156,159,168,171,180,183,216,219,243
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=med
#SBATCH --mem-per-cpu=3G
#SBATCH --time=8:00:00
#SBATCH --output=logs/slurm_%x_%A_%a.out
#SBATCH --error=logs/slurm_%x_%A_%a.err

set -euo pipefail
# Slurm runs a spool copy, so locate the original script from the submit directory.
submit_dir="${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
if [[ -f "${submit_dir}/slurm_cubic_beta_sweep.sh" ]]; then
    sweep="${submit_dir}/slurm_cubic_beta_sweep.sh"
elif [[ -f "${submit_dir}/CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh" ]]; then
    sweep="${submit_dir}/CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh"
else
    echo "Cannot find slurm_cubic_beta_sweep.sh under ${submit_dir}." >&2
    exit 1
fi
export SAMPLE_PLAN="$(dirname "$sweep")/Timing/cubic_beta_sweep/sample_plan_missing_zero_N2_N7.tsv"
exec bash "$sweep"
