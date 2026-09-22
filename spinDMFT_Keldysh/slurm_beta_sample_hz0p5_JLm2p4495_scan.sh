#!/usr/bin/env bash
# 6 beta values x 7 sample counts = 42 jobs at h_z=0.5 and JL=-2.4495.
# Submit from either the repository root or spinDMFT_Keldysh/ with:
#   sbatch spinDMFT_Keldysh/slurm_beta_sample_hz0p5_JLm2p4495_scan.sh
#   sbatch slurm_beta_sample_hz0p5_JLm2p4495_scan.sh
# To submit only a subset, for example tasks 0 through 20:
#   sbatch --array=0-20 slurm_beta_sample_hz0p5_JLm2p4495_scan.sh
#SBATCH --job-name=keldysh_beta_samples_JLm2p4495
#SBATCH --array=0-41
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=med
#SBATCH --mem-per-cpu=2G
#SBATCH --time=8:00:00
#SBATCH --output=logs/slurm_%x_%A_%a.out
#SBATCH --error=logs/slurm_%x_%A_%a.err

set -euo pipefail

# Each MPI rank uses its allocated CPUs for BLAS; avoid nested oversubscription.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

if [[ -x "${SLURM_SUBMIT_DIR}/spinDMFT_Keldysh/executable_DOUBLE.out" ]]; then
    KELDYSH_DIR="${SLURM_SUBMIT_DIR}/spinDMFT_Keldysh"
elif [[ -x "${SLURM_SUBMIT_DIR}/executable_DOUBLE.out" ]]; then
    KELDYSH_DIR="${SLURM_SUBMIT_DIR}"
else
    echo "Could not find spinDMFT_Keldysh/executable_DOUBLE.out under ${SLURM_SUBMIT_DIR}" >&2
    exit 1
fi

cd "${KELDYSH_DIR}"

# Keep exactly the beta and sample grids of slurm_beta_sample_scan.sh.
betas=(0.2 0.5 1 1.5 2 2.5)
samples=(1000 5000 10000 50000 100000 500000 1000000)
project=Keldysh_stype_C_hz0.5_JLm2p4495

task_id="${SLURM_ARRAY_TASK_ID}"
num_samples="${#samples[@]}"
num_betas="${#betas[@]}"
expected_tasks=$(( num_samples * num_betas ))
if (( task_id < 0 || task_id >= expected_tasks )); then
    echo "SLURM_ARRAY_TASK_ID=${task_id} is outside 0-$(( expected_tasks - 1 ))" >&2
    exit 1
fi

sample_index=$(( task_id % num_samples ))
beta_index=$(( task_id / num_samples ))

beta="${betas[beta_index]}"
samples_per_core="${samples[sample_index]}"

run_args=(
    --spinmodel=ISO
    --beta="${beta}"
    --JL=-2.4495
    --numImagTimeSteps=50
    --numRealTimeSteps=200
    --Tmax=15
    --gaussianFactorization=dense
    --numSamplesPerCore="${samples_per_core}"
    --project="${project}"
    --fileext="samples_per_core=${samples_per_core}"
    --spinInsertionStrategy=prefix
    --samplingStrategy=independent
    --cstype=C
    --Bname=z
    --Babs=0.5
)

echo "Starting task ${task_id}: beta=${beta}, samples/core=${samples_per_core}, h_z=0.5, JL=-2.4495"
mpirun -n 16 ./executable_DOUBLE.out "${run_args[@]}"
echo "Finished task ${task_id}: beta=${beta}, samples/core=${samples_per_core}, h_z=0.5, JL=-2.4495"
