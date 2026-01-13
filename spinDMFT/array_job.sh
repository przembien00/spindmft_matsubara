# Import data from file:
file="temperature_array.txt"
if [ -e "$file" ]; then
	mapfile -t data < "$file"
else
	echo "File not found: $file"
	exit 0
fi

echo "STARTING ARRAY JOB : SPINDMFT CHANGING T"
len=${#data[@]}
for((i=0;i<$len;i++))
do
	d="${data[$i]}"
	echo "RUNNING JOB WITH BETA = "$d
	mpirun -n 4 executable_DOUBLE.out --beta=$d --numSamplesPerCore=100000 --numSamplesPerSet=100 --cstype=A --spinmodel=ISO --numTimeSteps=50 --critneg=0.001 --JQ=1 --project="Matsubara"
wait
done
echo "DONE"
