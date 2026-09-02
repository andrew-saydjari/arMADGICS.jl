## Unit tests for WorkupSerial.jl (W2 serial tier) — synthetic mini-corpora
## only. Included from runtests.jl AFTER the fixture helpers are defined.

include(joinpath(@__DIR__, "..", "WorkupSerial.jl"))
using .WorkupSerial
using Distributed

const QUIET = s -> nothing

"Read one batch fixture's full contents for spot checks."
function read_fixture(dir, id::BatchId)
    p = joinpath(dir, lpad(id.adjfiberindx, 3, "0"), batch_filename(id))
    out = Dict{String, Any}()
    h5open(p, "r") do f
        for k in keys(f)
            k == "hdr" && continue
            out[k] = read(f[k])
        end
    end
    return out
end

@testset "workup_serial: full synthetic pipeline (inline)" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        out = joinpath(dir, "wu")
        plan = run_workup(dir, out; batchsize = BSZ, log = QUIET)
        @test plan.nsamp_out == 15
        @test plan.row_offset == 0
        @test isempty(plan.missing)
        # every key present, right shape, and rows land at identity ranges
        for (k, v) in plan.keyinfo
            p = out_file_path(out, k)
            @test isfile(p)
            h5open(p, "r") do f
                @test size(f[k]) == (v.shape[1:(end - 1)]..., 15)
                @test eltype(f[k]) == v.eltype
                # hdr provenance copied onto the hdr dataset
                @test read_attribute(f["hdr"], "git_commit") == "deadbeef"
                @test read_attribute(f["hdr"], "pipeline") == "arMADGICS.jl"
                @test read_attribute(f, "workup_row_offset") == 0
                # spot-check rows: batch (2,4) → rows 4:6; (7,4) → rows 14:15
                for (id, rng) in ((BatchId(2, 4), 4:6), (BatchId(7, 4), 14:15),
                                  (BatchId(5, 1), 8:10), (BatchId(2, 7), 7:7))
                    fx = read_fixture(dir, id)[k]
                    got = ndims(f[k]) == 1 ? f[k][rng] : f[k][:, rng]
                    @test isequal(got, fx)
                end
                @test !haskey(f, "missing_row_mask")   # complete corpus
            end
        end
        # checkpoint removed on success
        @test !isfile(joinpath(out, "workup_serial.ckpt"))
    end
end

@testset "workup_serial: missing batch — hard error vs allow_missing mask" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir; drop = [BatchId(2, 4)])
        out = joinpath(dir, "wu")
        # default: hard error, nothing silently shifted
        @test_throws ErrorException run_workup(dir, out; batchsize = BSZ, log = QUIET)
        # allow_missing: sentinel fill + mask
        out2 = joinpath(dir, "wu2")
        plan = run_workup(dir, out2; batchsize = BSZ, allow_missing = true, log = QUIET)
        @test plan.missing == [BatchId(2, 4)]
        h5open(out_file_path(out2, "flux"), "r") do f
            msk = read(f["missing_row_mask"])
            @test findall(msk .== 1) == [4, 5, 6]
            fx = read(f["flux"])
            @test all(isnan, fx[:, 4:6])
            # neighbors untouched
            @test isequal(fx[:, 1:3], read_fixture(dir, BatchId(2, 1))["flux"])
            @test isequal(fx[:, 7:7], read_fixture(dir, BatchId(2, 7))["flux"])
        end
        h5open(out_file_path(out2, "RV_flag"), "r") do f
            @test all(read(f["RV_flag"])[4:6] .== typemin(Int))
        end
        h5open(out_file_path(out2, "snr"), "r") do f
            @test all(isnan, read(f["snr"])[4:6])
        end
    end
end

@testset "workup_serial: truncated / corrupt batch detection aborts" begin
    # (a) byte-level truncation → unreadable
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        p = joinpath(dir, "005", batch_filename(BatchId(5, 1)))
        bytes = read(p)
        open(io -> write(io, bytes[1:200]), p, "w")
        err = try
            run_workup(dir, joinpath(dir, "wu"); batchsize = BSZ, log = QUIET)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("integrity", err.msg)
        @test occursin(batch_filename(BatchId(5, 1)), err.msg)
        @test occursin("unreadable", err.msg)
    end
    # (b) mid-write kill signature: valid HDF5 file but missing keys (the real
    #     087/arMADGICS_fiber_087_batch_0048201.h5 failure mode)
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        p = joinpath(dir, "007", batch_filename(BatchId(7, 4)))
        h5open(p, "w") do io
            io["flux"] = rand(11, 1)             # wrong row count too (expect 2)
            io["hdr"] = "This is only a header"
        end
        err = try
            run_workup(dir, joinpath(dir, "wu"); batchsize = BSZ, log = QUIET)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin(batch_filename(BatchId(7, 4)), err.msg)
        @test occursin("missing key snr", err.msg)
        @test occursin("last-axis length 1 != expected batch length 2", err.msg)
    end
end

@testset "workup_serial: contiguous fiber window" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        out = joinpath(dir, "wu7")
        plan = run_workup(dir, out; batchsize = BSZ, fibers = 7:7, log = QUIET)
        @test plan.nsamp_out == 5
        @test plan.row_offset == 10
        h5open(out_file_path(out, "flux"), "r") do f
            @test size(f["flux"]) == (11, 5)
            @test isequal(f["flux"][:, 1:3], read_fixture(dir, BatchId(7, 1))["flux"])
            @test isequal(f["flux"][:, 4:5], read_fixture(dir, BatchId(7, 4))["flux"])
            @test read_attribute(f, "workup_row_offset") == 10
            @test read_attribute(f, "workup_fibers") == [7, 7]
        end
        # window containing fibers 2+5: rows 1:10
        out2 = joinpath(dir, "wu25")
        plan2 = run_workup(dir, out2; batchsize = BSZ, fibers = 2:5, log = QUIET)
        @test plan2.nsamp_out == 10 && plan2.row_offset == 0
        h5open(out_file_path(out2, "snr"), "r") do f
            @test isequal(f["snr"][8:10], read_fixture(dir, BatchId(5, 1))["snr"])
        end
    end
end

@testset "workup_serial: resume skips checkpointed batches, rewrites the rest" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        out = joinpath(dir, "wu")
        run_workup(dir, out; batchsize = BSZ, log = QUIET)
        truth = Dict(k => h5read(out_file_path(out, k), k)
                     for k in ("flux", "snr", "RV_flag", "adjfiberindx"))
        # simulate a crash: pretend all batches except (2,4) were checkpointed,
        # and clobber rows of a checkpointed batch (5,1)→rows 8:10 (must stay
        # clobbered = skipped) and of the pending batch (2,4)→rows 4:6 (must be
        # rewritten with real data)
        allids = [BatchId(2, 1), BatchId(2, 4), BatchId(2, 7), BatchId(5, 1),
                  BatchId(7, 1), BatchId(7, 4)]
        open(joinpath(out, "workup_serial.ckpt"), "w") do io
            for id in allids
                id == BatchId(2, 4) && continue
                println(io, batch_filename(id))
            end
        end
        h5open(out_file_path(out, "snr"), "r+") do f
            f["snr"][4:6] = fill(-999.0, 3)
            f["snr"][8:10] = fill(-777.0, 3)
        end
        run_workup(dir, out; batchsize = BSZ, resume = true, log = QUIET)
        snr = h5read(out_file_path(out, "snr"), "snr")
        @test isequal(snr[4:6], truth["snr"][4:6])       # rewritten
        @test all(snr[8:10] .== -777.0)                  # skipped (checkpointed)
        @test !isfile(joinpath(out, "workup_serial.ckpt"))
    end
end

@testset "workup_serial: distributed readers produce identical output" begin
    mktempdir() do dir
        fli, fidx = synth_corpus(dir)
        out_inline = joinpath(dir, "wu_inline")
        run_workup(dir, out_inline; batchsize = BSZ, log = QUIET)

        procs_added = addprocs(2; exeflags = "--project=$(Base.active_project())")
        try
            rcp = abspath(joinpath(@__DIR__, "..", "RowContract.jl"))
            wsp = abspath(joinpath(@__DIR__, "..", "WorkupSerial.jl"))
            for p in procs_added
                remotecall_fetch(Base.include, p, Main, rcp)
                remotecall_fetch(Base.include, p, Main, wsp)
            end
            out_dist = joinpath(dir, "wu_dist")
            plan = run_workup(dir, out_dist; batchsize = BSZ, log = QUIET)
            for k in keys(plan.keyinfo)
                a = h5read(out_file_path(out_inline, k), k)
                b = h5read(out_file_path(out_dist, k), k)
                @test isequal(a, b)
            end
        finally
            rmprocs(procs_added)
        end
    end
end

println("All WorkupSerial tests passed.")
