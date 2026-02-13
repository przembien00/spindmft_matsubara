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
	mpirun -n 8 executable_DOUBLE.out --beta=$d --Bname=z --Babs=0.5 --JL=2 --numSamplesPerCore=500000 --numSamplesPerSet=100 --cstype=C --spinmodel=ISO --numTimeSteps=199 --critneg=0.1 --project="Paper"
wait
done
echo "DONE"


