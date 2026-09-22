#!/usr/bin/env bash
# 6 beta values x 5 unique grids x 3 physical cases = 90 jobs.
# The real-time branch varies N_R at fixed N_tau=80.  The imaginary-time
# branch varies N_tau at fixed N_R=200; its N_tau=80 baseline is shared with
# the real-time branch and is therefore submitted only once.
# Submit from either the repository root or spinDMFT_Keldysh/ with:
#   sbatch spinDMFT_Keldysh/slurm_beta_discretization_separate_scan.sh
#   sbatch slurm_beta_discretization_separate_scan.sh
# To submit only a subset, for example tasks 0 through 29:
#   sbatch --array=0-29 slurm_beta_discretization_separate_scan.sh
#SBATCH --job-name=keldysh_beta_discretization_separate
#SBATCH --array=0-89
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

betas=(0.2 0.5 1 1.5 2 2.5)
# N_R is varied while holding N_tau=80.
real_time_steps=(100 200 400)
# N_tau is varied while holding N_R=200.  N_tau=80 is omitted because the
# (N_tau, N_R)=(80, 200) baseline is already in real_time_steps above.
imag_time_steps=(40 160)
samples_per_core=100000

# Case order: h_z=0, JL=0; h_z=0.5, JL=0; h_z=0.5, JL=-2.4495.
projects=(
    Keldysh_stype_A_discretization_separate
    Keldysh_stype_C_hz0.5_discretization_separate
    Keldysh_stype_C_hz0.5_JLm2p4495_discretization_separate
)

task_id="${SLURM_ARRAY_TASK_ID}"
num_cases="${#projects[@]}"
num_real_grids="${#real_time_steps[@]}"
num_imag_grids="${#imag_time_steps[@]}"
num_grids=$(( num_real_grids + num_imag_grids ))
num_betas="${#betas[@]}"
expected_tasks=$(( num_cases * num_grids * num_betas ))
if (( task_id < 0 || task_id >= expected_tasks )); then
    echo "SLURM_ARRAY_TASK_ID=${task_id} is outside 0-$(( expected_tasks - 1 ))" >&2
    exit 1
fi

case_index=$(( task_id % num_cases ))
grid_index=$(( (task_id / num_cases) % num_grids ))
beta_index=$(( task_id / (num_cases * num_grids) ))

beta="${betas[beta_index]}"
project="${projects[case_index]}"
if (( grid_index < num_real_grids )); then
    scan=real_time
    num_imag=80
    num_real="${real_time_steps[grid_index]}"
else
    scan=imaginary_time
    imag_index=$(( grid_index - num_real_grids ))
    num_imag="${imag_time_steps[imag_index]}"
    num_real=200
fi

run_args=(
    --spinmodel=ISO
    --beta="${beta}"
    --numImagTimeSteps="${num_imag}"
    --numRealTimeSteps="${num_real}"
    --Tmax=15
    --gaussianFactorization=dense
    --numSamplesPerCore="${samples_per_core}"
    --project="${project}"
    --fileext="scan=${scan}__numImagTimeSteps=${num_imag}__numRealTimeSteps=${num_real}"
    --spinInsertionStrategy=prefix
    --samplingStrategy=independent
)

case "${case_index}" in
    0)
        run_args+=(--JL=0 --cstype=A)
        ;;
    1)
        run_args+=(--JL=0 --cstype=C --Bname=z --Babs=0.5)
        ;;
    2)
        run_args+=(--JL=-2.4495 --cstype=C --Bname=z --Babs=0.5)
        ;;
esac

echo "Starting task ${task_id}: beta=${beta}, scan=${scan}, numImagTimeSteps=${num_imag}, numRealTimeSteps=${num_real}, project=${project}"
mpirun -n 16 ./executable_DOUBLE.out "${run_args[@]}"
echo "Finished task ${task_id}: beta=${beta}, scan=${scan}, numImagTimeSteps=${num_imag}, numRealTimeSteps=${num_real}, project=${project}"
