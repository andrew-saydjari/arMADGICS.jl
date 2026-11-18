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

#
# Safe to source from a `set -euo pipefail` shell: the site lmod scripts and the
# `module` function read unset variables, so they run with nounset/errexit off and
# the caller's options are restored afterwards (scripts/lib_lmod.sh).

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../scripts" && pwd)/lib_lmod.sh"

lmod_relax_opts
if ! type module &>/dev/null; then
    for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
        # shellcheck disable=SC1090
        [ -r "$f" ] && source "$f" && break
    done
fi
lmod_restore_opts

lmod_module load cephtweaks
lmod_module load openmpi/5.0.6
lmod_module load hdf5/mpi-1.14.5

export CEPHTWEAKS_LAZYIO=1
export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}

echo "mpi_env ready: $(which mpiexec)"
echo "  HDF5_USE_FILE_LOCKING=$HDF5_USE_FILE_LOCKING CEPHTWEAKS_LAZYIO=$CEPHTWEAKS_LAZYIO"
