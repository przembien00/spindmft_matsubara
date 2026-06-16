#!/bin/bash

EXE="./executable_MH_DOUBLE.out"
START=$(date +%s%N)

for beta in 1 3 5; do
    echo "======= Running beta=$beta ======="
    iter_start=$(date +%s%N)
    
    mpirun -n 1 $EXE \
        --numSamplesPerCore=2000 \
        --numTimeSteps=40 \
        --beta=$beta \
        --mhBurnIn=100 \
        --mhStepSize=0.15 \
        --noselfcons \
        --reliterror=10000 \
        2>&1 | tail -15
    
    iter_end=$(date +%s%N)
    iter_time=$(( (iter_end - iter_start) / 1000000 ))
    echo "beta=$beta took ${iter_time}ms"
    echo ""
done

END=$(date +%s%N)
TOTAL=$(( (END - START) / 1000000 ))
echo "Total runtime: ${TOTAL}ms ($(( TOTAL / 1000 ))s)"
