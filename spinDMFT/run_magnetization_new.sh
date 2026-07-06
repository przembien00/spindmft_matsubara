#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "=== STAGE 1: h_z = 0.1 ==="
for beta in $(seq 2.5 0.1 3.5); do
    echo "RUNNING JOB WITH BETA = $beta, Babs=0.1"
    mpirun -n 4 ./executable_DOUBLE.out \
        --spinmodel=ISO --cstype=C --JL=-2 --beta=$beta \
        --numTimeSteps=99 --numSamplesPerCore=100000 \
        --mhStepSize=0.9 --mhBurnIn=100 --iterlimit=30 \
        --numBlocks=20 --critneg=100 --Bname=z --Babs=0.1 \
        --project=Magnetization_new
done

echo "=== STAGE 2: h_z = 0.001, using h_z=0.1 result as initial condition ==="
for beta in $(seq 2.5 0.1 3.5); do
    echo "RUNNING JOB WITH BETA = $beta, Babs=0.001 (initialized from Babs=0.1)"
    initfile="Magnetization_new/spinmodel=ISO__JL=-2__beta=${beta}__h=z_h_abs=0.1"
    mpirun -n 4 ./executable_DOUBLE.out \
        --spinmodel=ISO --cstype=C --JL=-2 --beta=$beta \
        --numTimeSteps=99 --numSamplesPerCore=100000 \
        --mhStepSize=0.9 --mhBurnIn=100 --iterlimit=30 \
        --numBlocks=20 --critneg=100 --Bname=z --Babs=0.001 \
        --loadinit --initcorrfile="$initfile" \
        --project=Magnetization_new
done
echo "DONE"
