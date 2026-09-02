# workup/mpi_env — MPI tier environment (W2)

Separate Julia project binding MPI.jl + HDF5.jl to the CLUSTER MODULE stack
(system OpenMPI + parallel HDF5), for `workup_mpi.jl` and `mpio_smoke.jl`.
It is deliberately isolated from the main arMADGICS.jl project and from
`workup/Project.toml`: those keep using the serial HDF5_jll artifact, and
mixing module libraries into their manifests would poison them.

## Stack (pins recorded 2026-09-02, ccalin051 / rocky9, modules/2.4-20250724)

| component | version | source |
|---|---|---|
| Julia | 1.11.6 (juliaup) | |
| OpenMPI | 5.0.6 | `module load openmpi/5.0.6` → `/mnt/sw/nix/store/k1aj11vmsmwf411qyp4l57lk9di2x929-openmpi-5.0.6` |
| HDF5 (parallel) | 1.14.5 | `module load hdf5/mpi-1.14.5` → `/mnt/sw/nix/store/mh5b4q8vqxcdcmirjnr18v3r6g4ix8by-hdf5-1.14.5` |
| cephtweaks | LD_PRELOAD `libcephtweaks.so`, `CEPHTWEAKS_LAZYIO=1` | `module load cephtweaks` |
| MPI.jl | 0.20.27 | Manifest.toml |
| HDF5.jl | 0.17.3 | Manifest.toml |
| MPIPreferences.jl | 0.1.12 | Manifest.toml |

This is the same module stack acasey's Python MPI workup used (his
`setup.sh`: cephtweaks + `CEPHTWEAKS_LAZYIO=1`, openmpi/5.0.6,
hdf5/mpi-1.14.5, `HDF5_MPI=ON CC=mpicc` for the h5py source build — the
"sneaky module loading" needed for ceph).

## How the binding works (the exact incantation)

1. `MPIPreferences.use_system_binary()` with `openmpi/5.0.6` loaded →
   `[MPIPreferences]` section of `LocalPreferences.toml`
   (binary="system", abi="OpenMPI", libmpi="libmpi", mpiexec="mpiexec").
2. HDF5.jl must use the module's PARALLEL libhdf5 instead of the serial
   HDF5_jll artifact. Two subtleties (both hit in practice):
   * `HDF5.API.set_libraries!` can't run while `using HDF5` fails, and it
     fails because MPIPreferences abi=OpenMPI makes HDF5_jll resolve its
     `mpi+openmpi` artifact variant, whose `OpenMPI_jll.__init__` dlopens the
     ARTIFACT `libmpi.so` — which crashes against the module's
     `LD_LIBRARY_PATH` (`undefined symbol: opal_single_threaded`). So the
     JLLs themselves must be overridden to the system libraries via the
     JLLWrappers `<product>_path` preferences — see the `[OpenMPI_jll]` and
     `[HDF5_jll]` sections of `LocalPreferences.toml`.
   * Preferences for TRANSITIVE deps only take effect if the package appears
     in this project's `[extras]` (see Project.toml).
3. Verify with:
   `julia --project=. -e 'using MPI, HDF5; @show MPI.identify_implementation() HDF5.has_parallel()'`
   → `("OpenMPI", v"5.0.6")`, `has_parallel = true`.

If the module system moves the nix store paths (module version bump), rerun
step 1 and update the paths in `LocalPreferences.toml`.

## Usage

```bash
source workup/mpi_env/setup_mpi_env.sh          # modules + locking env
mpiexec -n 4 julia --project=workup/mpi_env workup/mpi_env/mpio_smoke.jl \
    --out /mnt/ceph/.../smoke.h5
mpiexec -n 4 julia --project=workup/mpi_env workup/workup_mpi.jl \
    --rawdir RAW --outdir OUT --fibers 1:1
```

`HDF5_USE_FILE_LOCKING`: HDF5 file locking on ceph is the known hazard for
multi-process HDF5; `setup_mpi_env.sh` defaults it to `FALSE`. The
`mpio_smoke.jl` matrix (locking TRUE/FALSE/unset) documents what actually
works on this filesystem — see the W2 PR text for measured results.
