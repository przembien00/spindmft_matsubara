#!/usr/bin/env bash
#SBATCH --job-name=sdmft_field_beta2.5
#SBATCH --array=0-10
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=short
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# This works whether sbatch is called from the repository root or spinDMFT/.
if [[ -x "${SLURM_SUBMIT_DIR}/spinDMFT/executable_DOUBLE.out" ]]; then
    SPINDMFT_DIR="${SLURM_SUBMIT_DIR}/spinDMFT"
elif [[ -x "${SLURM_SUBMIT_DIR}/executable_DOUBLE.out" ]]; then
    SPINDMFT_DIR="${SLURM_SUBMIT_DIR}"
else
    echo "Could not find executable_DOUBLE.out under ${SLURM_SUBMIT_DIR}" >&2
    exit 1
fi

cd "${SPINDMFT_DIR}"

# The h_z=0.1 point already exists. Tasks 9-19 use it as the initial
# condition for the additional low-field series.
fields=(
#    0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.15 0.125
    0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.005 0.001
)
field="${fields[${SLURM_ARRAY_TASK_ID}]}"

run_args=(
    --spinmodel=ISO
    --cstype=C
    --JL=-2.449489743
    --beta=2.5
    --Bname=z
    --Babs="-${field}"
    --numTimeSteps=30
    --numSamplesPerCore=500000
    --numSamplesPerSet=100
    --mhStepSize=0.9
    --mhBurnIn=100
    --numBlocks=20
    --iterlimit=30
    --critneg=100
    --project=Magnetization_field_scan
)

if (( SLURM_ARRAY_TASK_ID >= 9 )); then
    initfile="Magnetization_field_scan/spinmodel=ISO__JL=-2.4495__beta=2.5__h=z_h_abs=-0.1"
    run_args+=(--loadinit --initcorrfile="${initfile}")
fi

echo "Starting spinDMFT: beta=1.5, h_z=${field}, task=${SLURM_ARRAY_TASK_ID}"

mpirun -n 16 ./executable_DOUBLE.out "${run_args[@]}"

echo "Finished spinDMFT: beta=1.5, h_z=${field}"
