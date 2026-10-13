#!/usr/bin/env bash
#SBATCH --job-name=staged-sweep
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=med
#SBATCH --array=0-17
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")"

betas=(0.2 0.5 1.0 1.5 2.0 2.5)
jls=(0 -2.4495)

# 18 array tasks: 12 stype C with Bz=0.5 (both JL values), then 6 stype A
# with no field and JL=0. Each group spans all six beta values.
task_id=${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}
if (( task_id < 0 || task_id > 17 )); then
    echo "SLURM_ARRAY_TASK_ID must be in 0..17 (got ${task_id})" >&2
    exit 2
fi

if (( task_id < 12 )); then
    cstype=C
    field_args=(--Bname=z --Babs=0.5 --dontrescaleB)
    local_id=$task_id
else
    cstype=A
    field_args=(--Bname=none)
    local_id=$((task_id - 12))
fi

beta=${betas[$((local_id % 6))]}
if [[ ${cstype} == C ]]; then
    JL=${jls[$((local_id / 6))]}
else
    JL=0
fi

echo "Task ${task_id}: cstype=${cstype}, beta=${beta}, JL=${JL}"
echo "MPI ranks: ${SLURM_NTASKS:-16}; samples per rank: 10000000"

srun --ntasks=16 --cpus-per-task=1 ./executable_DOUBLE.out \
    --spinmodel=ISO \
    --cstype="${cstype}" \
    "${field_args[@]}" \
    --JL="${JL}" \
    --beta="${beta}" \
    --numImagTimeSteps=50 \
    --numRealTimeSteps=200 \
    --Tmax=15 \
    --numSamplesPerCore=10000000 \
    --preparationSamplesPerCore=100000 \
    --iterlimit=30 \
    --realTimeEvolution=staggered-strang \
    --gaussianFactorization=dense
