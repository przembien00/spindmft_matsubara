#!/usr/bin/env bash
# 6 beta values x 3 discretizations = 18 jobs.
# Submit from either the repository root or spinDMFT_Keldysh/ with:
#   sbatch spinDMFT_Keldysh/slurm_beta_discretization_antithetic_scan.sh
#   sbatch slurm_beta_discretization_antithetic_scan.sh
# To submit only a subset, for example tasks 0 through 5:
#   sbatch --array=0-5 slurm_beta_discretization_antithetic_scan.sh
#SBATCH --job-name=keldysh_beta_discretization
#SBATCH --array=0-17
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=med
#SBATCH --mem-per-cpu=2G
#SBATCH --time=8:00:00
#SBATCH --output=logs/slurm_%x_%A_%a.out
#SBATCH --error=logs/slurm_%x_%A_%a.err

set -euo pipefail

if [[ -x "${SLURM_SUBMIT_DIR}/spinDMFT_Keldysh/executable_DOUBLE.out" ]]; then
    KELDYSH_DIR="${SLURM_SUBMIT_DIR}/spinDMFT_Keldysh"
elif [[ -x "${SLURM_SUBMIT_DIR}/executable_DOUBLE.out" ]]; then
    KELDYSH_DIR="${SLURM_SUBMIT_DIR}"
else
    echo "Could not find spinDMFT_Keldysh/executable_DOUBLE.out under ${SLURM_SUBMIT_DIR}" >&2
    exit 1
fi

cd "${KELDYSH_DIR}"

# Match the beta range of slurm_beta_sample_antithetic_scan.sh.  The existing
# production grid is (50, 200); each grid below is finer at fixed Tmax=15.
betas=(0.2 0.5 1 1.5 2 2.5)
num_imag_time_steps=(75 100 150)
num_real_time_steps=(300 400 600)
samples_per_core=100000
project=Keldysh_stype_C_hz0.5_antithetic_discretization

task_id="${SLURM_ARRAY_TASK_ID}"
num_betas="${#betas[@]}"
num_grids="${#num_imag_time_steps[@]}"
expected_tasks=$(( num_betas * num_grids ))
if (( task_id < 0 || task_id >= expected_tasks )); then
    echo "SLURM_ARRAY_TASK_ID=${task_id} is outside 0-$(( expected_tasks - 1 ))" >&2
    exit 1
fi

grid_index=$(( task_id % num_grids ))
beta_index=$(( task_id / num_grids ))

beta="${betas[beta_index]}"
num_imag="${num_imag_time_steps[grid_index]}"
num_real="${num_real_time_steps[grid_index]}"

run_args=(
    --spinmodel=ISO
    --beta="${beta}"
    --numImagTimeSteps="${num_imag}"
    --numRealTimeSteps="${num_real}"
    --Tmax=15
    --gaussianFactorization=fft
    --fftCrossFrequencyCutoff=10
    --numSamplesPerCore="${samples_per_core}"
    --project="${project}"
    --fileext="numImagTimeSteps=${num_imag}__numRealTimeSteps=${num_real}"
    --spinInsertionStrategy=prefix
    --samplingStrategy=independent
    --cstype=C
    --Bname=z
    --Babs=0.5
)

echo "Starting task ${task_id}: beta=${beta}, numImagTimeSteps=${num_imag}, numRealTimeSteps=${num_real}, samples/core=${samples_per_core}"
mpirun -n 16 ./executable_DOUBLE.out "${run_args[@]}"
echo "Finished task ${task_id}: beta=${beta}, numImagTimeSteps=${num_imag}, numRealTimeSteps=${num_real}"
