fext="Test"
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --cstype=A --spinmodel=ISO --numTimeSteps=200 --dt=0.02 --project=Tests --fileext=$fext"_1"
echo ""
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --cstype=B --spinmodel=DRF --numTimeSteps=200 --dt=0.1 --project=Tests --fileext=$fext"_2"
echo ""
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --cstype=C --spinmodel=DRF --numTimeSteps=200 --dt=0.1 --Bname=z --project=Tests --fileext=$fext"_3"
echo ""
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numSamplesPerSet=100 --cstype=D --spinmodel=DRF --numTimeSteps=200 --dt=0.1 --project=Tests --fileext=$fext"_4"
echo -e "\n\n=== Plotting ==="
PYTHONPATH=../python_libs python Tests.py
cd Data/Tests
rm *Test_*