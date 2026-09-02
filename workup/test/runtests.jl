## Unit tests for RowContract.jl (W1) — synthetic fixtures only, no real data.
using Test
using HDF5

include(joinpath(@__DIR__, "..", "RowContract.jl"))
using .RowContract

# ---------------------------------------------------------------- fixtures --
# Tiny synthetic corpus: batchsize 3; fibers 2 (7 rows), 5 (3 rows), 7 (5 rows)
# => batches: (2,1),(2,4),(2,7), (5,1), (7,1),(7,4); nsamp = 15.
const BSZ = 3
const FIBSPEC = [(2, 7), (5, 3), (7, 5)]

function synth_identity()
    tele = String[]; mjd = String[]; expnum = Int[]; fib = Int[]; sid = Int[]
    for (f, n) in FIBSPEC, i in 1:n
        push!(tele, isodd(i) ? "apo" : "lco")
        push!(mjd, string(60000 + i))
        push!(expnum, 10 * f + i)
        push!(fib, f)
        push!(sid, 1000 * f + i)
    end
    return (tele = tele, mjd = mjd, expnum = expnum, adjfiberindx = fib, sdss_id = sid)
end

function write_full_list_info(path, ident)
    h5open(path, "w") do io
        io["tele"] = ident.tele
        io["mjd"] = ident.mjd
        io["expnum"] = ident.expnum
        io["adjfiberindx"] = ident.adjfiberindx
        io["sdss_id"] = ident.sdss_id
    end
end

function write_batch_info(path, ident)
    open(path, "w") do io
        println(io, "# Batch information for arMADGICS pipeline")
        println(io, "# Format: linear_index, tele, mjd, expnum, adjfiberindx ")
        println(io, "# Total batches: 6")
        println(io, "")
        li = 0; prevfib = -1
        for r in 1:length(ident.tele)
            f = ident.adjfiberindx[r]
            li = (f == prevfib) ? li + 1 : 1
            prevfib = f
            println(io, "$li, $(ident.tele[r]), $(ident.mjd[r]), $(ident.expnum[r]), $f")
        end
    end
end

function write_batch_file(dir, id::BatchId, nrow; npix = 11)
    sub = joinpath(dir, lpad(id.adjfiberindx, 3, "0"))
    mkpath(sub)
    p = joinpath(sub, batch_filename(id))
    h5open(p, "w") do io
        io["flux"] = rand(npix, nrow)             # stored {nrow, npix} in file
        io["snr"] = rand(nrow)
        io["RV_flag"] = rand(Int, nrow)
        io["adjfiberindx"] = fill(id.adjfiberindx, nrow)
        io["hdr"] = "This is only a header"
        attrs = HDF5.attributes(io["hdr"])
        attrs["git_commit"] = "deadbeef"
        attrs["pipeline"] = "arMADGICS.jl"
    end
    return p
end

function synth_corpus(dir; drop = BatchId[], extra = BatchId[])
    ident = synth_identity()
    write_full_list_info(joinpath(dir, "full_list_info.h5"), ident)
    write_batch_info(joinpath(dir, "batch_info.txt"), ident)
    fli = load_full_list_info(joinpath(dir, "full_list_info.h5"))
    fidx = verify_fiber_blocks(fli; batchsize = BSZ)
    for id in expected_batches(fidx)
        id in drop && continue
        write_batch_file(dir, id, length(batch_within_range(fidx, id)))
    end
    for id in extra
        n = id.adjfiberindx <= 600 && fidx.counts[id.adjfiberindx] > 0 ?
            BSZ : BSZ
        write_batch_file(dir, id, n)
    end
    return fli, fidx
end

# ------------------------------------------------------------------- tests --

@testset "filename parsing" begin
    id = parse_batch_filename("arMADGICS_fiber_007_batch_0000301.h5")
    @test id == BatchId(7, 301)
    @test parse_batch_filename(batch_filename(BatchId(600, 1234567))) ==
          BatchId(600, 1234567)
    @test isnothing(parse_batch_filename("full_list_info.h5"))
    @test isnothing(parse_batch_filename("arMADGICS_fiber_7_batch_301.h5"))
    @test batch_relpath(BatchId(7, 301)) ==
          joinpath("007", "arMADGICS_fiber_007_batch_0000301.h5")
end

@testset "identity load + fiber-block verification" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        @test fli.nsamp == 15
        @test tele_name(fli.tele[1]) == "apo"
        @test fli.mjd[1] == 60001
        @test fidx.counts[2] == 7 && fidx.counts[5] == 3 && fidx.counts[7] == 5
        @test fidx.offsets[2] == 0 && fidx.offsets[5] == 7 && fidx.offsets[7] == 10
        @test fiber_row_range(fidx, 7) == 11:15
        # global row ranges from identity
        @test batch_row_range(fidx, BatchId(2, 4)) == 4:6
        @test batch_row_range(fidx, BatchId(2, 7)) == 7:7      # short last batch
        @test batch_row_range(fidx, BatchId(7, 4)) == 14:15    # short last batch
        @test batch_of_linear_index(fidx, 7, 5) == BatchId(7, 4)
        @test_throws ErrorException batch_row_range(fidx, BatchId(2, 5)) # off-grid
        @test_throws ErrorException batch_row_range(fidx, BatchId(3, 1)) # empty fiber
        @test_throws ErrorException batch_row_range(fidx, BatchId(2, 10)) # past end
    end
end

@testset "non-contiguous fiber blocks rejected" begin
    ident = synth_identity()
    bad = collect(ident.adjfiberindx)
    bad[2] = 7   # fiber 7 row inside the fiber-2 block
    mktempdir() do dir
        p = joinpath(dir, "full_list_info.h5")
        write_full_list_info(p, (tele = ident.tele, mjd = ident.mjd,
            expnum = ident.expnum, adjfiberindx = bad, sdss_id = ident.sdss_id))
        fli = load_full_list_info(p)
        @test_throws ErrorException verify_fiber_blocks(fli; batchsize = BSZ)
    end
end

@testset "batch_info crosscheck" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        cc = crosscheck_batch_info(joinpath(dir, "batch_info.txt"), fli, fidx)
        @test cc.nrows == 15
        @test cc.n_mismatch == 0
        @test cc.header_total_batches == 6
        # corrupt one row: wrong expnum + wrong linear_index
        lines = readlines(joinpath(dir, "batch_info.txt"))
        lines[10] = "99, apo, 60001, 999, 5"
        open(joinpath(dir, "batch_info.txt"), "w") do io
            foreach(l -> println(io, l), lines)
        end
        cc2 = crosscheck_batch_info(joinpath(dir, "batch_info.txt"), fli, fidx)
        @test cc2.n_mismatch == 1
        @test occursin("row 6", cc2.first_mismatches[1])
    end
end

@testset "reconcile: complete corpus" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        disc = discover_batches(dir)
        @test length(disc) == 6
        rec = reconcile_batches(disc, expected_batches(fidx))
        @test length(rec.present) == 6
        @test isempty(rec.missing) && isempty(rec.extra)
    end
end

@testset "reconcile: missing batch" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir; drop = [BatchId(2, 4)])
        disc = discover_batches(dir)
        exp = expected_batches(fidx)
        # hard error by default, message lists the missing file
        err = try
            reconcile_batches(disc, exp)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("arMADGICS_fiber_002_batch_0000004.h5", err.msg)
        @test occursin("MISSING", err.msg)
        # allow_missing: mask + sentinel spec instead of silent shift
        rec = reconcile_batches(disc, exp; allow_missing = true)
        @test rec.missing == [BatchId(2, 4)]
        msk = missing_rows_mask(fidx, rec.missing)
        @test findall(msk) == [4, 5, 6]
        ki = discover_keys(rec.paths[BatchId(2, 1)])
        fs = fill_spec(ki)
        @test isnan(fs["flux"]) && isnan(fs["snr"])
        @test fs["RV_flag"] == typemin(Int)
    end
end

@testset "reconcile: extra batch always fatal" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir; extra = [BatchId(9, 1)])
        disc = discover_batches(dir)
        err = try
            reconcile_batches(disc, expected_batches(fidx); allow_missing = true)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("arMADGICS_fiber_009_batch_0000001.h5", err.msg)
        @test occursin("extra", err.msg)
    end
end

@testset "out-of-order / shuffled discovery is irrelevant" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        disc = discover_batches(dir)
        # rebuild the discovery dict in scrambled insertion order
        ids = collect(keys(disc))
        scrambled = Dict{BatchId, String}()
        for id in reverse(sort(ids))
            scrambled[id] = disc[id]
        end
        rec = reconcile_batches(scrambled, expected_batches(fidx))
        @test rec.present == sort(ids)
        # row ranges depend only on identity, never on discovery order
        @test [batch_row_range(fidx, id) for id in rec.present] ==
              [1:3, 4:6, 7:7, 8:10, 11:13, 14:15]
    end
end

@testset "key discovery + hdr dataset attrs + integrity" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        disc = discover_batches(dir)
        p = disc[BatchId(2, 1)]
        ki = discover_keys(p)
        @test Set(keys(ki)) == Set(["flux", "snr", "RV_flag", "adjfiberindx"])
        @test ki["flux"].shape == (11, 3) && ki["flux"].eltype == Float64
        @test ki["snr"].shape == (3,)
        hdr = read_hdr_attrs(p)
        @test hdr["git_commit"] == "deadbeef"
        @test hdr["pipeline"] == "arMADGICS.jl"
        # clean batch
        @test isempty(check_batch_integrity(p, BatchId(2, 1), fidx; ref_keyinfo = ki))
        # short last batch checks against its own expected length
        p7 = disc[BatchId(7, 4)]
        @test isempty(check_batch_integrity(p7, BatchId(7, 4), fidx; ref_keyinfo = ki))
        # wrong-length batch detected
        bad = write_batch_file(dir, BatchId(5, 1), 2)  # fiber 5 batch should have 3 rows
        probs = check_batch_integrity(bad, BatchId(5, 1), fidx; ref_keyinfo = ki)
        @test any(occursin("last-axis length 2 != expected batch length 3", s) for s in probs)
    end
end

@testset "fiber-filtered discovery + nrow-based integrity method" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        disc7 = discover_batches(dir; fibers = 7:7)
        @test sort(collect(keys(disc7))) == [BatchId(7, 1), BatchId(7, 4)]
        disc25 = discover_batches(dir; fibers = 2:5)
        @test length(disc25) == 4
        # nrow-based method agrees with the FiberIndex-based one
        p = disc7[BatchId(7, 4)]
        ki = discover_keys(p)
        @test check_batch_integrity(p, 2; ref_keyinfo = ki) ==
              check_batch_integrity(p, BatchId(7, 4), fidx; ref_keyinfo = ki)
        @test !isempty(check_batch_integrity(p, 3; ref_keyinfo = ki))
    end
end

println("All RowContract tests passed.")

include(joinpath(@__DIR__, "test_workup_serial.jl"))
