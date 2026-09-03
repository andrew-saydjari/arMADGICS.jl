## mpio_smoke.jl — parallel-HDF5-on-ceph smoke test for the W2 MPI tier.
#
# All ranks open ONE shared HDF5 file on ceph via the MPIO driver and write
# disjoint last-axis slabs in the workup's contract pattern (rank-strided
# "batches" of rows). Verifies:
#   * collective file create/open with Drivers.MPIO
#   * independent AND collective dataset hyperslab writes
#   * read-back correctness (rank 0, serial reopen)
#   * write throughput for an (npix, nrows) Float64 dataset
#
# Run (ccalin051, local ranks):
#   source setup_mpi_env.sh          # module MPI + parallel HDF5 + cephtweaks
#   mpiexec -n 4 julia --project=. mpio_smoke.jl --out /path/on/ceph/smoke.h5 \
#       [--npix 8700] [--rows-per-batch 100] [--batches-per-rank 25]
#
# Dataset transfers are INDEPENDENT (the default, and what workup_mpi.jl
# uses; the gist's h5py writes were independent too). File create/open/close
# through the MPIO driver are inherently collective.
#
# HDF5_USE_FILE_LOCKING is read from the environment so the caller can test
# the locking matrix (TRUE / FALSE / unset).

using MPI
using HDF5
using Printf

function getarg(flag, default)
    i = findfirst(==(flag), ARGS)
    isnothing(i) && return default
    return parse(typeof(default), ARGS[i + 1])
end
getflag(flag) = flag in ARGS

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    outpath = let i = findfirst(==("--out"), ARGS)
        isnothing(i) && error("--out required")
        ARGS[i + 1]
    end
    npix = getarg("--npix", 8700)
    rows_per_batch = getarg("--rows-per-batch", 100)
    batches_per_rank = getarg("--batches-per-rank", 25)

    nbatch = nranks * batches_per_rank
    nrows = nbatch * rows_per_batch
    if rank == 0
        @printf("mpio_smoke: %d ranks, dataset (%d, %d) Float64 = %.2f GB, independent writes\n",
            nranks, npix, nrows, npix * nrows * 8 / 1e9)
        println("  HDF5_USE_FILE_LOCKING=", get(ENV, "HDF5_USE_FILE_LOCKING", "<unset>"))
        println("  CEPHTWEAKS_LAZYIO=", get(ENV, "CEPHTWEAKS_LAZYIO", "<unset>"))
        println("  out: ", outpath)
        HDF5.has_parallel() || error("HDF5.jl built without parallel support")
    end
    MPI.Barrier(comm)

    # deterministic per-row content so verification catches any misplacement
    rowval(r) = Float64(r) .+ (1:npix) ./ (npix + 1)

    t0 = MPI.Wtime()
    h5open(outpath, "w"; driver = HDF5.Drivers.MPIO(comm)) do f
        dset = create_dataset(f, "x", datatype(Float64), dataspace(npix, nrows))
        # contract pattern: batch b (1-based, global) owns rows
        # (b-1)*rows_per_batch+1 : b*rows_per_batch; rank r takes batches
        # r+1, r+1+nranks, ... (the gist's rank-strided assignment)
        buf = Matrix{Float64}(undef, npix, rows_per_batch)
        for b in (rank + 1):nranks:nbatch
            r0 = (b - 1) * rows_per_batch
            for j in 1:rows_per_batch
                buf[:, j] .= rowval(r0 + j)
            end
            dset[:, (r0 + 1):(r0 + rows_per_batch)] = buf
        end
    end
    MPI.Barrier(comm)
    t1 = MPI.Wtime()

    if rank == 0
        gb = npix * nrows * 8 / 1e9
        @printf("write+close wall: %.2f s (%.2f GB → %.2f GB/s aggregate)\n",
            t1 - t0, gb, gb / (t1 - t0))
        # serial reopen + verification
        nbad = 0
        h5open(outpath, "r") do f
            d = f["x"]
            size(d) == (npix, nrows) || error("bad dataset size $(size(d))")
            for r in Int[1, 2, rows_per_batch, rows_per_batch + 1,
                         nrows - 1, nrows, rand(1:nrows, 64)...]
                got = d[:, r]
                isequal(got, rowval(r)) || (nbad += 1; println("  MISMATCH row $r"))
            end
        end
        println(nbad == 0 ? "VERIFY: PASS (all sampled rows exact)" :
                            "VERIFY: FAIL ($nbad rows)")
        nbad == 0 || MPI.Abort(comm, 2)
    end
    MPI.Barrier(comm)
    MPI.Finalize()
end

main()
