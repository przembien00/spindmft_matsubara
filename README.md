# Collection of algorithms concerning spin dynamic mean-field theory

## Description:
This repository contains implementations of the numerical methods spinDMFT, CspinDMFT and nl-spinDMFT. The codes are written entirely in C++.

## 1. Installation of required programms and libraries
The installation is yet only explained for Linux operating systems.

| programm/library      | description                                                                                                       |
| --------------------- | ----------------------------------------------------------------------------------------------------------------- |
| GCC 13.2.0            | C++ compiler collection, C++17 has to be available                                                                |
| MPI 4.1.6             | for parallelized computations                                                                                     |
| Blaze 3.8             | efficient header-only C++ library for matrices and vectors, see https://bitbucket.org/blaze-lib/blaze/src/master/ |
| HDF5 1.10.10          | C++ library for handling hdf5 files in C++, hdf5: efficient file format allowing to store meta data               |
| Boost 1.83            | C++ library, required here for parameter handling at runtime via with boost programm options                      |
| BLAS & LAPACK 3.12.0  | C++ libraries, required here for blaze (e.g., for diagonalization)                                                |

If BLAS & LAPACK are not available, diagonalizations can alternatively be carried out with eigen (another C++ library, see https://eigen.tuxfamily. org/index.php?title=Main_Page). This may affect the performance, see https://bitbucket.org/blaze-lib/blaze/src/master/ for more information. How to employ eigen will be addressed below under "Altering CMakeLists.txt".

## 2. Compiling the codes regularly
Once all requirements have been installed, the codes can be compiled using cmake. The procedure for this is explained on the example of spinDMFT in the following, but it works exactly the same for CspinDMFT and nl-spinDMFT. Open a terminal and navigate to the directory "spinDMFT/Algorithm". Execute the following commands:
```sh
cp CMakeLists.txt_ CMakeLists.txt
mkdir build
cd build
cmake ..
make
```
The executable file "spinDMFT/executable_DOUBLE.out" should have been generated, if the compilation was successful.

FFAST MATH?!

## 3. Altering CMakeLists.txt
The compilation can be adjusted by altering the file "CMakeLists.txt". The options under the bullet point "# things to be adapted" may be changed if desired: One can switch from the data type DOUBLE to FLOAT under "set(USE_TYPE DOUBLE)" which reduces the computation time and numerical accuracy. If BLAS & LAPACK are not available, Blaze cannot be installed properly. However, Blaze itself is a header-only library so that most of its functionality can still be used if the repository is just downloaded. Depending on the location of the Blaze directory, it may be automatically found in the compilation process. If not, one can provide the path to the blaze directory manually under "include_directories(/path/to/blaze)".Specifically, diagonalizations cannot be performed by blaze, if LAPACK is not available. In this case, Eigen may be installed as an alternative. To force the use of Eigen for the necessary diagonalizations, one has to enable "add_definitions(-DEIGEN)".

## 4. Testing the codes
The main directory of each code contains a simple testing script named "Tests.sh", which may be carried out via 
```sh
./Tests.sh
```
once the executable file "executable_DOUBLE.out" exists. The script runs several quick simulations and plots their output using python and matplotlib (which obviously need to be installed for this to work). If the code is changed or extended, the testing script may serve as a quick check. The plots show a comparison between fixed reference data and newly generated test data. The deviation between reference and test has to be small (but not exactly zero because the Monte-Carlo simulation implies a statistical deviation). Note that the executed tests do not cover all options and parameters of the code and, therefore, do not guarantee the absense of any bugs. Note also that the Test data are automatically deleted at the end of the script.

## 5. Usage of the code spinDMFT
The executables "executable_XXX.out" with XXX=FLOAT or DOUBLE can be carried out using mpirun. At first one may test the help command via
```sh
mpirun executable_XXX.out --help
```
A help menu should pop up in the terminal. It shows the available options and parameters to run the code. The syntax for running an actual simulation is 
```sh
mpirun -n Y executable_XXX.out --parameter1=value1 --parameter2=value2 --option1 --option2
```
Here, "Y" should be replaced by the desired number of cores to be used, and "parameter", "value" and "option" accordingly. Several parameters and options can be added in this way. They may be added in arbitrary order to the mpirun command. As a specific example, one may run
```sh
mpirun -n 4 executable_DOUBLE.out --numSamplesPerCore=2500 --numTimeSteps=100 --dt=0.02 --cstype=A --spinmodel=ISO
```
This command starts a simulation of spinDMFT on 4 cores with 2500 samples per core, 100 time steps, a time step width of 0.02 considering an isotropic Heisenberg model. A small description of the simulation should appear followed by some progress bars. Once the simulation is finished, one should find an hdf5-file containing the corresponding simulation data in the Data directory. This file contains the bare simulation results as well as a lot of meta data including numerical details, time measurements of the simulation and the initially set parameters. 

### 5.1 Reading the data
The data of the hdf5-file can be read out using 
```sh
h5dump filename.hdf5
```
if h5dump is installed. For a better overview, one can alternatively use python. To this end, navigate to the main directory of the repository and execute
```sh
python quickreader.py <filename.hdf5>
```
replacing <filename.hdf5> by the actual name of the hdf5-file including the correct relative path. Note that the python library h5py is required for the script to work. The script prints all attributes and datasets to the terminal. For the datasets, only the shape is shown. One can execute
```sh
python quickreader.py <filename.hdf5> showdsets
```
if one desires to see the content of the datasets.

## 6. Usage of the codes CspinDMFT and nl-spinDMFT
In contrast to spinDMFT, the simulations of CspinDMFT and nl-spinDMFT are more complex and cannot be started immediately. Only the help command is easily accessible via 
```sh
mpirun executable_XXX.out --help
```
replacing "XXX" by FLOAT or DOUBLE. For an actual simulation, one requires specific configuration files that carry information about the spin cluster and the spin-mean-field couplings to be considered. These files are stored in the subfolder "Configuration_Data" of the corresponding directory of each code. More information on these files and how to generate them can be found under [6.1 Creating configuration files](#61-creating-configuration-files). Once a configuration file exists, the simulations for CspinDMFT or nl-spinDMFT can be started the same way as for spinDMFT via the syntax
```sh
mpirun -n Y executable_XXX.out --config=<config_file> --parameter1=value1 --parameter2=value2 --option1 --option2
```
replacing <config_file> by the corresponding file without ".hdf5"-ending and without the path (if the configuration file is located in a subdirectory of "Configuration_Data", this is specified via the "--project" option). The simulation will then be carried out utilizing the data in the configuration file. Note that providing the flag "--config=<config_file>" is mandatory.

To see how the data from CspinDMFT and nl-spinDMFT are stored, one may consider the Reference Data required for the test scripts "Tests.sh" stored in "Data/Tests". To this end, one can use h5dump or python as described in [5.1 Reading the data](#51-reading-the-data).

### 6.1 Creating configuration files



## 8. Functionality of the code CspinDMFT
in CspinDMFT 4D tensor JCRA, ij, kl
why? -> k and l from 1 to numSpins have a direct physical meaning; in nl-spinDMFT environment and cluster are not necessarily related


## Roadmap
Adding a full documentation of the code CspinDMFT and nl-spinDMFT is planned.

add local extra int to physics.h

spinmodels always global (not yet adapted in C and nl-spinDMFT)

python_libs?

potential efficiency boost covmatrix diagonalization not parallelized
still many copies between CspinDMFT and nl-spinDMFT
better iterator usage in computation of correlations S_X S_Y and S_Z
extend tests?

MPI overhead time not measured

spin flt to str should also accept 3o2 e.g.

adaptive discretization!

adaptive num samples per core not implemented in spinDMFT code, why?

TEST CODE ON CLUSTER?!


## Project status
The project is updated frequently. 

## Contributing

## Support
Feel free to contact the author by sending an e-mail to "timo.graesser@tu-dortmund.de" for any questions or suggestions. Note that the code has been written by a physicist not software engineer. Any suggestions for more clarity or higher efficiency are appreciated. 

## License



## Extending the Codes

### Generally
how to include new parameters:
1) add it in the header file Parameter_Space.h
2) add it in Parameter_Space.cpp, so that it can be adjusted with boost program options
3a) perhaps add it to the essential parameters print function
3b) add it to the filename string if its a physical parameter
4) add it to void HDF5_Storage::store_main in Storage_Concept.cpp
5) include the parameter in the code

----- Options -----

efficiency test with: perf record, perf report



Guideline for Naming:
classes, structs:           'UpperCamelCase'

functions, methods:         'lower_case'

variable instances:         varying
-> int, uint:               'lowerCamelCase' (often 'num_UpperCamelCase')
-> double:                  'lower_case'
-> string, char:            'lower_case'
-> class/struct inst:       'lower_case' or 'my_lower_case'

global variables:     'UPPER_CASE'

error functions:            'UPPER_CASE'

name spaces:                'Upper_case'






nl-spinDMFT
1) Importing mean-fields 
Tha algorithm imports spin correlations, which can be transformed into mean-field correlations...

explanation:
in case of loaded correlations: correlation tensor cluster should in the first step be independent
of numSpins (because it relates only to the environment)
1) import desired set of correlations (e.g. "11", "12", ...)
2) give each correlation an index k (k could be just the position of the correlation)
3) compute mean-field ViVj from transformation that maps from k to ij (new routine for this, not the standard self-consistency)
