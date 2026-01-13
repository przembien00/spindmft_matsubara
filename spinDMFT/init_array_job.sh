# Import data from file:
file="J_array.txt"
if [ -e "$file" ]; then
    mapfile -t data < "$file"
else
    echo "File not found: $file"
    exit 0
fi

echo "STARTING ARRAY JOB : SPINDMFT CHANGING T"
len=${#data[@]}
o=19
for((i=0;i<$len;i++))
do
    d="${data[$i]}"
    echo "RUNNING JOB WITH BETA = "$d
    mpirun -n 4 executable_DOUBLE.out --beta=1 --JQ=$d --numSamplesPerCore=1000000 --numSamplesPerSet=100 --cstype=A --spinmodel=ISO --numTimeSteps=300 --loadinit --project="Iterative_Init" --initcorrfile="Iterative_Init/spinmodel=ISO__JQ="$o"__beta=1" --critneg=0.1 --iterlimit=30
    o=$d
done
echo "DONE"