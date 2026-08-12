#!/usr/bin/env bash
#SBATCH --job-name=sdmft_hz0.001_beta
#SBATCH --array=0-40
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

betas=(
    1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95
    2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45
    2.5 2.55 2.6 2.65 2.7 2.75 2.8 2.85 2.9 2.95
    3.0 3.05 3.1 3.15 3.2 3.25 3.3 3.35 3.4 3.45
    3.5
)
beta="${betas[${SLURM_ARRAY_TASK_ID}]}"

run_args=(
    --spinmodel=ISO
    --cstype=C
    --JL=-2.449489743
    --beta="${beta}"
    --Bname=z
    --Babs=-0.001
    --numTimeSteps=30
    --numSamplesPerCore=500000
    --numSamplesPerSet=100
    --mhStepSize=0.9
    --mhBurnIn=100
    --numBlocks=20
    --iterlimit=30
    --critneg=100
    --project=Magnetization_beta_scan
)

initfile="Magnetization_beta_scan/spinmodel=ISO__JL=-2.4495__beta=${beta}__h=z_h_abs=-0.1"
run_args+=(--loadinit --initcorrfile="${initfile}")

echo "Starting spinDMFT: beta=${beta}, h_z=0.001, initialized from h_z=0.1, task=${SLURM_ARRAY_TASK_ID}"

mpirun -n 16 ./executable_DOUBLE.out "${run_args[@]}"

echo "Finished spinDMFT: beta=${beta}, h_z=0.001"
