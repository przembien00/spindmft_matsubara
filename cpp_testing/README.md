# Shared checks

Algorithm source and build files may use their own directory and common
libraries (`cpp_libs`, `python_libs`, external numerical libraries and shared
test support). They must not include or compile files from another algorithm,
including archived algorithms. Put reusable implementation in the common
libraries and reusable test fixtures here.

Run the repository-wide dependency check with:

```sh
python3 cpp_testing/check_algorithm_dependencies.py
```

It checks C++ includes, CMake/Makefile paths and source symlinks, including
common libraries that might indirectly introduce solver dependencies. It
excludes generated builds, data and analysis output. Optional initialization
from another simulation's data and comparison plots do not make the solvers
depend on each other's implementation.

The Keldysh and finite-temperature cluster CTest suites also run this check.
To run it independently, copy `CMakeLists.txt_` to `CMakeLists.txt`, configure
this directory with CMake and run CTest.
`common/test_correlation_weights.cpp` is compiled separately against each
finite-temperature cluster solver's own mean-field implementation.

Cross-algorithm configuration-generator and migration checks also live here:

```sh
python3 -B -m unittest discover -s cpp_testing/common -p test_configuration_weights.py
```

These integration checks require NumPy and h5py and load each generator from
its owning algorithm. No algorithm uses another algorithm's generator.
