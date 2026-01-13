echo -e "=== Creating Configs ==="
cd Configuration
make Tests
cd ../
echo -e "\n=== Performing Tests ==="
#fext="--fileext=Reference_Data"
fext="--fileext=Test_Data"
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --stype=B --numTimeSteps=200 --dt=0.1 --spinspinmodel=DRF --spinmfmodel=DRF --initmode=import --impcorrfile="Tests/spinmodel=DRF_Reference_2" --project=Tests --config=Test_1 $fext
echo ""
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --stype=B --numTimeSteps=200 --dt=0.1 --spinspinmodel=DRF --spinmfmodel=Ising --initmode=import --impcorrfile="Tests/spinmodel=DRF_Reference_2" --project=Tests --config=Test_2 $fext
echo ""
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=200 --numSamplesPerSet=50 --stype=B --numTimeSteps=50 --dt=0.2 --extrapolate --spinspinmodel=DRF --spinmfmodel=DRF --initmode=import --impcorrfile="Tests/spinmodel=DRF_Reference_2" --project=Tests --config=Test_3 $fext
echo ""
echo -e "\n=== Plotting ==="
PYTHONPATH=../python_libs python Tests.py
cd Data/Tests
rm -f *Test_Data*