#!/bin/bash
# Environment for the W2 MPI tier (workup_mpi.jl, mpio_smoke.jl).
# Source before ANY julia --project=<this dir> invocation:
#   source workup/mpi_env/setup_mpi_env.sh
#
# This mirrors the stack acasey's Python MPI workup ran with (his setup.sh):
#   cephtweaks (LD_PRELOAD libcephtweaks.so) + CEPHTWEAKS_LAZYIO=1
#   openmpi/5.0.6 + hdf5/mpi-1.14.5
# LocalPreferences.toml pins MPI.jl and HDF5.jl to exactly these module
# libraries (system binaries, never the JLL artifacts).
#
# HDF5_USE_FILE_LOCKING: ceph supports POSIX locks, but HDF5 file locking
# across MPI ranks on ceph is the known hazard; the mpio_smoke.jl matrix
# tests locking on/off. Default here: FALSE (safe on ceph, and what the
# smoke results support).

module load cephtweaks
module load openmpi/5.0.6
module load hdf5/mpi-1.14.5

export CEPHTWEAKS_LAZYIO=1
export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}

echo "mpi_env ready: $(which mpiexec)"
echo "  HDF5_USE_FILE_LOCKING=$HDF5_USE_FILE_LOCKING CEPHTWEAKS_LAZYIO=$CEPHTWEAKS_LAZYIO"
