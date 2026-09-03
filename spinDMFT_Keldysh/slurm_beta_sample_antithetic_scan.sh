#!/usr/bin/env bash
# 6 beta values x 7 sample counts x 4 cases = 168 jobs.
# Submit from either the repository root or spinDMFT_Keldysh/ with:
#   sbatch spinDMFT_Keldysh/slurm_beta_sample_antithetic_scan.sh
#   sbatch slurm_beta_sample_antithetic_scan.sh
# To submit only a subset, for example tasks 0 through 41:
#   sbatch --array=0-41 slurm_beta_sample_antithetic_scan.sh
#SBATCH --job-name=keldysh_beta_samples
#SBATCH --array=0-167
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=short
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%A_%a.out
#SBATCH --error=slurm_%x_%A_%a.err

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

betas=(0.2 0.5 1 1.5 2 2.5)
samples=(1000 5000 10000 50000 100000 500000 1000000)

# Case order: A, antithetic A, C with h_z=0.5, antithetic C with h_z=0.5.
projects=(
    Keldysh_stype_A
    Keldysh_stype_A_antithetic
    Keldysh_stype_C_hz0.5
    Keldysh_stype_C_hz0.5_antithetic
)

task_id="${SLURM_ARRAY_TASK_ID}"
num_cases="${#projects[@]}"
num_samples="${#samples[@]}"
num_betas="${#betas[@]}"
expected_tasks=$(( num_cases * num_samples * num_betas ))
if (( task_id < 0 || task_id >= expected_tasks )); then
    echo "SLURM_ARRAY_TASK_ID=${task_id} is outside 0-$(( expected_tasks - 1 ))" >&2
    exit 1
fi

case_index=$(( task_id % num_cases ))
grid_index=$(( task_id / num_cases ))
sample_index=$(( grid_index % num_samples ))
beta_index=$(( grid_index / num_samples ))

beta="${betas[beta_index]}"
samples_per_core="${samples[sample_index]}"
project="${projects[case_index]}"

run_args=(
    --spinmodel=ISO
    --beta="${beta}"
    --numImagTimeSteps=50
    --numRealTimeSteps=200
    --Tmax=15
    --gaussianFactorization=fft
    --fftCrossFrequencyCutoff=10
    --numSamplesPerCore="${samples_per_core}"
    --project="${project}"
    --fileext="samples_per_core=${samples_per_core}"
)

case "${case_index}" in
    0)
        run_args+=(--cstype=A)
        ;;
    1)
        run_args+=(--cstype=A --antitheticPairs)
        ;;
    2)
        run_args+=(--cstype=C --Bname=z --Babs=0.5)
        ;;
    3)
        run_args+=(--cstype=C --Bname=z --Babs=0.5 --antitheticPairs)
        ;;
esac

echo "Starting task ${task_id}: beta=${beta}, samples/core=${samples_per_core}, project=${project}"
mpirun -n 16 ./executable_DOUBLE.out "${run_args[@]}"
echo "Finished task ${task_id}: beta=${beta}, samples/core=${samples_per_core}, project=${project}"
