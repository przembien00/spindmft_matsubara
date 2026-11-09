#!/usr/bin/env bash
# 4 configurations x 3 field cases x 11 beta values = 132 independent jobs.
# Recommended submission (creates logs before Slurm opens its output files):
#   bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit
# Inspect a command without running it:
#   bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --dry-run 0
#SBATCH --job-name=cluster_keldysh_cubic
#SBATCH --array=0-131
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=med
#SBATCH --mem-per-cpu=5G
#SBATCH --time=6:00:00
#SBATCH --output=logs/slurm_%x_%A_%a.out
#SBATCH --error=logs/slurm_%x_%A_%a.err

set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Slurm executes a copy from its spool directory, not the original script path.
if [[ ! -f "${script_dir}/benchmark_cubic_sweep.py" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    if [[ -f "${SLURM_SUBMIT_DIR}/CspinDMFT_Keldysh/benchmark_cubic_sweep.py" ]]; then
        script_dir="${SLURM_SUBMIT_DIR}/CspinDMFT_Keldysh"
    elif [[ -f "${SLURM_SUBMIT_DIR}/benchmark_cubic_sweep.py" ]]; then
        script_dir="${SLURM_SUBMIT_DIR}"
    else
        echo "Cannot find CspinDMFT_Keldysh under ${SLURM_SUBMIT_DIR}." >&2
        exit 1
    fi
fi
if [[ "${1:-}" == --submit ]]; then
    shift
    mkdir -p "${script_dir}/logs"
    cd "${script_dir}"
    plan="${SAMPLE_PLAN:-${script_dir}/Timing/cubic_beta_sweep/sample_plan.tsv}"
    [[ -r "$plan" ]] || { echo "Missing sample plan: $plan" >&2; exit 1; }
    # An explicit array selection uses the script's conservative memory default.
    # Otherwise reserve 2 GB/rank for N=2 and 5 GB/rank for N=4.
    for option in "$@"; do
        if [[ "$option" == --array || "$option" == --array=* || "$option" == -a || "$option" == -a* ]]; then
            exec sbatch "$@" "${script_dir}/slurm_cubic_beta_sweep.sh"
        fi
    done
    total_cases=$(( $(wc -l < "$plan") - 1 ))
    for size in 2 4; do
        case_indices=()
        index=0
        while IFS=$'\t' read -r configuration rest; do
            [[ "$configuration" == config ]] && continue
            if [[ "$configuration" == *"_N=${size}_"* ]]; then
                case_indices+=("$index")
            fi
            index=$(( index + 1 ))
        done < "$plan"
        (( ${#case_indices[@]} )) || continue
        indices=()
        for beta_index in {0..10}; do
            for case_index in "${case_indices[@]}"; do
                indices+=("$(( beta_index * total_cases + case_index ))")
            done
        done
        array_spec="$(IFS=,; echo "${indices[*]}")"
        memory=5G
        [[ "$size" != 2 ]] || memory=2G
        sbatch --array="$array_spec" --mem-per-cpu="$memory" \
            --job-name="cluster_keldysh_cubic_N${size}" "$@" \
            "${script_dir}/slurm_cubic_beta_sweep.sh"
    done
    exit 0
fi

dry_run=false
if [[ "${1:-}" == --dry-run ]]; then
    dry_run=true
    task_id="${2:?Usage: $0 --dry-run TASK_ID}"
else
    task_id="${SLURM_ARRAY_TASK_ID:?Use --submit or --dry-run TASK_ID}"
fi
plan="${SAMPLE_PLAN:-${script_dir}/Timing/cubic_beta_sweep/sample_plan.tsv}"
[[ -r "$plan" ]] || { echo "Missing sample plan: $plan. Run benchmark_cubic_sweep.py first." >&2; exit 1; }
total_cases=$(( $(wc -l < "$plan") - 1 ))
total_tasks=$(( total_cases * 11 ))
[[ "$task_id" =~ ^[0-9]+$ ]] && (( task_id < total_tasks )) || {
    echo "Task ID must be an integer from 0 through $(( total_tasks - 1 ))." >&2; exit 1;
}

betas=(0.2 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5)
beta="${betas[task_id / total_cases]}"
row_index=$(( task_id % total_cases + 2 ))
IFS=$'\t' read -r config physical_case symmetry field amplitude samples prep_samples prep_limit estimate rest \
    < <(sed -n "${row_index}p" "$plan")
rest="${rest%$'\r'}"
if [[ "$rest" != *$'\t'ready ]] && ! "$dry_run"; then
    echo "${config}/${physical_case}: sample plan exceeds its time or memory allowance. Recalibrate before submitting this case." >&2
    exit 1
fi
[[ "$samples" =~ ^[0-9]+$ && "$prep_samples" =~ ^[0-9]+$ && "$prep_limit" =~ ^[0-9]+$ ]] || {
    echo 'Invalid sample counts in plan.' >&2; exit 1;
}
[[ -f "${script_dir}/Configuration_Data/CubicLattice/${config}.hdf5" ]] || {
    echo "Missing configuration: $config" >&2; exit 1;
}
sign=AFM
[[ "$config" != *J=-* ]] || sign=FM
project="Cubic_beta_sweep/${sign}/${physical_case}"
ranks="${SLURM_NTASKS:-16}"
executable="${KELDYSH_EXECUTABLE:-${script_dir}/executable_DOUBLE.out}"
[[ -x "$executable" ]] || { echo "Missing executable: $executable. Build the solver on the cluster first." >&2; exit 1; }
[[ "$ranks" -eq 16 ]] || { echo 'This plan is intended for 16 MPI ranks.' >&2; exit 1; }

# A sample count is PER RANK. Do not divide it by the number of ranks:
# all 16 ranks sample in parallel and jointly provide 16 times the statistics.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export BLAZE_NUM_THREADS=1
args=(
    --config="$config" --project=CubicLattice --dstproject="$project"
    --spinspinmodel=ISO --beta="$beta" --cstype="$symmetry" "--${field}=${amplitude}"
    --numImagTimeSteps=50 --numRealTimeSteps=150 --Tmax=10
    --realTimeSubsteps=1 --timeTranslationInvariant=true
    --samplingStrategy=independent --numSamplesPerSet=1 --numBlocks=32
    --numSamplesPerCore="$samples" --preparationSamplesPerCore="$prep_samples"
    --iterlimit=10 --preparationIterlimit="$prep_limit" --reliterror=5
    --fileext="cubic_sweep__ni=50__nr=150__case=${physical_case}__samples=${samples}__prep=${prep_samples}"
)
echo "Task ${task_id}: ${config}, beta=${beta}, ${physical_case}; samples/rank=${samples}, preparation/rank=${prep_samples}"
echo "Conservative timing estimate for up to ${prep_limit} preparation + 10 real-time iterations: ${estimate} seconds."
if "$dry_run"; then
    printf '%q ' mpirun -n "$ranks" "$executable" "${args[@]}"
    printf '\n'
    exit 0
fi
cd "$script_dir"
set +e
mpirun -n "$ranks" "$executable" "${args[@]}"
status=$?
set -e
if (( status == 3 )); then
    echo 'Saved unconverged output after the 10-iteration real-time limit.'
    status=0
elif (( status == 2 )); then
    echo 'Preparation failed to converge; real-time stage was not started.' >&2
fi
exit "$status"
