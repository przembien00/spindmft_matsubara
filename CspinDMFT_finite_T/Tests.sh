echo -e "=== Creating Configs ==="
cd Configuration
make Tests
cd ../
echo -e "\n=== Performing Tests ==="
#fext="--fileext=Reference_Data"
fext="--fileext=Test_Data"
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --stype=B --numTimeSteps=200 --dt=0.1 --spinspinmodel=DRF --project=Tests --config=Test_1 $fext
echo ""
mpirun -n 4 executable_DOUBLE.out --adaptive --staterrtolerance=0.0025 --numSamplesPerCore=100 --numSamplesPerSet=100 --stype=B --numTimeSteps=200 --dt=0.1 --spinspinmodel=DRF --project=Tests --config=Test_2 $fext
echo ""
echo -e "\n=== Plotting ==="
PYTHONPATH=../python_libs python Tests.py
cd Data/Tests
rm -f *Test_Data*
