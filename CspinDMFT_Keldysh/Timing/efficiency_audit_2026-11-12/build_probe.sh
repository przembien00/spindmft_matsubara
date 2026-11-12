#!/bin/sh
# Run from repository root, after building CspinDMFT_Keldysh/Algorithm/build.
# These include/library paths describe the local macOS installation used in this audit.
set -eu
python3 CspinDMFT_Keldysh/Timing/efficiency_audit_2026-11-12/make_probe.py
/usr/bin/c++ -O3 -DNDEBUG -std=gnu++17 -DUSE_DOUBLE \
 -I CspinDMFT_Keldysh/Algorithm -I cpp_libs \
 -isystem /opt/homebrew/include -isystem /opt/homebrew/Cellar/open-mpi/5.0.9/include \
 CspinDMFT_Keldysh/Timing/efficiency_audit_2026-11-12/probe_main.cpp \
 CspinDMFT_Keldysh/Timing/efficiency_audit_2026-11-12/probe_functions.cpp \
 CspinDMFT_Keldysh/Algorithm/build/libcluster_keldysh.a \
 -Wl,-rpath,/opt/homebrew/lib /opt/homebrew/Cellar/open-mpi/5.0.9/lib/libmpi.dylib \
 /opt/homebrew/lib/libhdf5.320.1.1.dylib -framework Accelerate -lm -ldl \
 /opt/homebrew/lib/libboost_program_options.dylib /opt/homebrew/lib/libboost_container.dylib \
 -o CspinDMFT_Keldysh/Timing/efficiency_audit_2026-11-12/probe.out
