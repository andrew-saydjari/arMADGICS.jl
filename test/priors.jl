# Unit tests for src/priors.jl (pass-1 runtime integration wiring points).
# Path-construction tests are unconditional; tests that open real prior files are
# guarded on ceph availability (CI runs without the prior store).

@testset "priors: build_prior_dict paths + env overrides" begin
    prior_dir_tst = "/tmp/arm_prior_root/"

    # defaults (no overrides): pass-1b starCont layout, split by telescope
    for k in ("ARM_STARCONT_PRIOR_DIR",)
        delete!(ENV, k)
    end
    pd = build_prior_dict(prior_dir_tst)
    @test pd["starCont_apo"] ==
          joinpath(prior_dir_tst, "2026_09_05/prior_outputs/starCont_pass1c", "built_apo", "APOGEE_starcont_svd_60_f")
    @test pd["starCont_lco"] ==
          joinpath(prior_dir_tst, "2026_09_05/prior_outputs/starCont_pass1c", "built_lco", "APOGEE_starcont_svd_60_f")
    @test haskey(pd, "chebmsk")
    @test haskey(pd, "starLines_refLSF")

    # env override redirects the starCont root (E6-style prior-swap runs)
    ENV["ARM_STARCONT_PRIOR_DIR"] = "/tmp/arm_alt_starcont"
    pd2 = build_prior_dict(prior_dir_tst)
    @test pd2["starCont_apo"] == "/tmp/arm_alt_starcont/built_apo/APOGEE_starcont_svd_60_f"
    delete!(ENV, "ARM_STARCONT_PRIOR_DIR")

    # per-fiber path completion + loud failure on a missing file
    @test_throws ErrorException per_fiber_prior_file("/tmp/arm_definitely_missing/pfx_f", 85)
    tdir = mktempdir()
    pfx = joinpath(tdir, "APOGEE_starcont_svd_60_f")
    touch(pfx * "085.h5")
    @test per_fiber_prior_file(pfx, 85) == pfx * "085.h5"
end

@testset "priors: ddstaronly is refused loudly (audit item 5)" begin
    # must throw at prior-load time regardless of file availability (the guard runs
    # before any file is opened) — on pre-integration main this configuration read
    # msk_starCor from a CLOSED file handle
    pd = build_prior_dict("/tmp/arm_prior_root/")
    @test_throws ErrorException load_fiber_priors(pd, 85; ddstaronly=true)
    err = try
        load_fiber_priors(pd, 85; ddstaronly=true)
    catch e
        sprint(showerror, e)
    end
    @test occursin("E7", err) # error message must point at the E7 deliverable
end

# contract tests against the real prior store (skipped off-cluster)
prior_dir_real = "/mnt/ceph/users/sdssv/work/asaydjari/"
pd_real = build_prior_dict(prior_dir_real)
if isfile(pd_real["chebmsk"]) && isfile(pd_real["starCont_apo"] * "085.h5") &&
   isfile(pd_real["skycont"] * "085.h5") && isfile(pd_real["skyLines_faint"] * "085.h5")
    @testset "priors: load_fiber_priors contract (real files)" begin
        pv = load_fiber_priors(pd_real, 85)
        @test length(pv) == 10
        chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont,
        V_starlines_refLSF, V_starlines, msk_starCor, V_skycont, V_skyline_faint = pv
        npix = 8700
        @test length(chebmsk_exp) == npix
        @test eltype(chebmsk_exp) == Bool
        @test size(V_starcont) == (npix, 60)
        @test all(isfinite, V_starcont)
        @test size(V_skycont) == (npix, 30)
        @test all(isfinite, V_skycont)
        @test size(V_skyline_faint) == (npix, 120)
        @test all(isfinite, V_skyline_faint)
        @test size(V_starlines_refLSF, 1) == npix
        @test V_starlines === V_starlines_refLSF # E7 not yet wired (TODO tracked)
        @test length(msk_starCor) == npix && all(msk_starCor)
        # masks are Bool and skymsk can only shrink chebmsk
        for m in (skymsk_bright, skymsk_faint, skymsk)
            @test eltype(m) == Bool && length(m) == npix
            @test all(m .<= chebmsk_exp)
        end
        @test skymsk == skymsk_faint # bright lines fully masked (no bright submask)
        # per-telescope mask switch (audit item 3): runtime chebmsk is the 2026_04_25
        # per-telescope mask (apo 7742 / lco 7833), not the 2025 global (7783)
        apo_msk, lco_msk = h5open(pd_real["chebmsk"]) do f
            convert.(Bool, read(f["apo"])), convert.(Bool, read(f["lco"]))
        end
        @test chebmsk_exp == apo_msk
        @test count(apo_msk) == 7742 && count(lco_msk) == 7833
        @test apo_msk != lco_msk
        # adjfiberindx range guard + telescope routing (lco file may not exist here;
        # routing is exercised via the error message path name)
        @test_throws ErrorException load_fiber_priors(pd_real, 0)
        @test_throws ErrorException load_fiber_priors(pd_real, 601)
    end
else
    @info "priors: skipping real-file contract tests (prior store not available)"
end
