# CMake generated Testfile for 
# Source directory: /Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm
# Build directory: /Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(complex_gaussian "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_complex_gaussian")
set_tests_properties(complex_gaussian PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;130;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
add_test(contour_kernel "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_contour_kernel")
set_tests_properties(contour_kernel PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;134;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
add_test(ratio_statistics "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_ratio_statistics")
set_tests_properties(ratio_statistics PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;138;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
add_test(keldysh_propagation "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_keldysh_propagation")
set_tests_properties(keldysh_propagation PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;146;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
add_test(full_contour_measurement "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_full_contour_measurement")
set_tests_properties(full_contour_measurement PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;154;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
add_test(full_contour_measurement_antithetic_option "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_full_contour_measurement" "--samplingStrategy=independent" "--antitheticPairs" "--numSamplesPerCore=2")
set_tests_properties(full_contour_measurement_antithetic_option PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;155;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
add_test(full_contour_measurement_antithetic_pcn_option "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/build_keldysh/test_full_contour_measurement" "--samplingStrategy=pcn" "--antitheticPairs" "--numSamplesPerCore=1")
set_tests_properties(full_contour_measurement_antithetic_pcn_option PROPERTIES  _BACKTRACE_TRIPLES "/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;159;add_test;/Users/przembien/Projects/spindmft_matsubara/spinDMFT_Keldysh/Algorithm/CMakeLists.txt;0;")
