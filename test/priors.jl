# Unit tests for src/priors.jl (pass-1 runtime integration wiring points).
# Path-construction tests are unconditional; tests that open real prior files are
# guarded on ceph availability (CI runs without the prior store).

@testset "priors: build_prior_dict paths + env overrides" begin
    prior_dir_tst = "/tmp/arm_prior_root/"

    # defaults (no overrides): pass-1c starCont layout, split by telescope.
    # Save/restore any caller-set override rather than clobbering it (the real-file
    # tests below honor the caller's environment).
    saved = Dict(k => get(ENV, k, nothing) for k in ("ARM_STARCONT_PRIOR_DIR",))
    for k in keys(saved)
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
    for (k, v) in saved
        isnothing(v) || (ENV[k] = v)
    end

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

@testset "priors: prior_support_mask unit behavior" begin
    npix = 100
    msk = falses(npix); msk[11:90] .= true
    V = zeros(npix, 3)
    V[msk, :] .= 1.0                       # unit row norm over the mask
    # uniform support: nothing IN the mask is dropped. Pixels OUTSIDE it are zero by
    # construction and are correctly reported unsupported, so the guard can only ever
    # shrink chebmsk_exp, never grow it.
    @test all(prior_support_mask(V, msk)[msk])
    @test !any(prior_support_mask(V, msk)[.!msk])
    # a pixel at 1e-3 of the median row norm is rejected at the default 1e-2 reltol,
    # and kept if the threshold is loosened below it
    V[20, :] .= 1e-3
    supp = prior_support_mask(V, msk)
    @test !supp[20]
    @test count(msk .& .!supp) == 1
    @test prior_support_mask(V, msk; reltol=1e-4)[20]
    # 3-D priors (npix, ncomp, nsub) are accepted (row norm over all trailing dims)
    V3 = ones(npix, 3, 2); V3[20, :, :] .= 1e-3
    @test !prior_support_mask(V3, trues(npix))[20]
    # loud failures rather than silent nonsense
    @test_throws ErrorException prior_support_mask(V, falses(npix)) # no power on mask
    @test_throws ErrorException prior_support_mask(V, trues(npix + 1)) # length mismatch
    # env plumbing
    @test prior_support_reltol() == 1e-2
    withenv("ARM_PRIOR_SUPPORT_RELTOL" => "0") do
        @test prior_support_reltol() == 0.0
    end
end

# contract tests against the real prior store (skipped off-cluster)
prior_dir_real = "/mnt/ceph/users/sdssv/work/asaydjari/"
pd_real = build_prior_dict(prior_dir_real)
if isfile(pd_real["chebmsk"]) && isfile(pd_real["starCont_apo"] * "085.h5") &&
   isfile(pd_real["skycont"] * "085.h5") && isfile(pd_real["skyLines_faint"] * "085.h5") &&
   isfile(pd_real["starLines_LSF"] * "085.h5")
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
        # E7: fit basis is the PER-FIBER prior (same shape family as refLSF),
        # report basis stays refLSF; distinct arrays in default mode
        @test size(V_starlines) == (npix, 50, 10)
        @test V_starlines !== V_starlines_refLSF
        @test V_starlines != V_starlines_refLSF
        @test all(isfinite, V_starlines)
        @test length(msk_starCor) == npix && all(msk_starCor)
        # hack-mode fallback must restore the pre-E7 behavior exactly
        ENV["ARM_STARLINES_REFLSF_HACK"] = "1"
        pv_hack = load_fiber_priors(pd_real, 85)
        delete!(ENV, "ARM_STARLINES_REFLSF_HACK")
        @test pv_hack[7] === pv_hack[6] # V_starlines aliases V_starlines_refLSF
        @test pv_hack[6] == V_starlines_refLSF
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
        # prior-support guard is a NO-OP at APO: bit-identical to the raw mask
        @test chebmsk_exp == apo_msk
        @test count(chebmsk_exp) == 7742
        @test count(apo_msk) == 7742 && count(lco_msk) == 7833
        @test apo_msk != lco_msk

        # THE INVARIANT THAT WOULD HAVE CAUGHT THE LCO EDGE BUG: it is not enough for
        # Vmat to be zero OFF its mask (the old check) — it must be non-negligible ON
        # it, or the solve fits pixels whose prior pins the component to ~0.
        relpow(V, m) = (rn = vec(sqrt.(sum(abs2, reshape(V, size(V, 1), :), dims=2)));
                        rn ./ median(view(rn, m)))
        for (nm, V) in (("starCont", V_starcont), ("skyCont", V_skycont))
            r = relpow(V, chebmsk_exp)
            @test minimum(r[chebmsk_exp]) >= 1e-2      # the guard's contract
            @test minimum(r[apo_msk]) > 0.1            # APO is clean on the RAW mask
            @test all(iszero, r[.!apo_msk])            # still zero off-mask
        end
        # adjfiberindx range guard + telescope routing
        @test_throws ErrorException load_fiber_priors(pd_real, 0)
        @test_throws ErrorException load_fiber_priors(pd_real, 601)
    end

    @testset "priors: prior-support guard drops the LCO dead-edge pixels" begin
        # LCO regression (2026_09_07): 67 in-mask px carry ~1e-3 of the median row
        # norm because the priors are built over the STATIC per-telescope mask. They
        # are identical on all 300 LCO fibers; 63 of them survive submsk_faint, so
        # only this guard removes them.
        adj_lco = 301
        if isfile(pd_real["starCont_lco"] * lpad(adj_lco, 3, "0") * ".h5") &&
           isfile(pd_real["skycont"] * lpad(adj_lco, 3, "0") * ".h5") &&
           isfile(pd_real["starLines_LSF"] * lpad(adj_lco, 3, "0") * ".h5")
            DEAD = vcat(276:320, 3638:3655, 8511:8514)
            @test length(DEAD) == 67
            lco_raw = h5open(pd_real["chebmsk"]) do f
                convert.(Bool, read(f["lco"]))
            end

            pv_lco = load_fiber_priors(pd_real, adj_lco)
            cheb_g, _, skymsk_faint_g, skymsk_g = pv_lco[1], pv_lco[2], pv_lco[3], pv_lco[4]
            # exactly the 67 known pixels are removed, and nothing else
            @test count(lco_raw) - count(cheb_g) == 67
            @test findall(lco_raw .& .!cheb_g) == DEAD
            @test !any(cheb_g[DEAD])
            @test !any(skymsk_g[DEAD]) && !any(skymsk_faint_g[DEAD])

            # ... and with the guard OFF, submsk_faint alone leaves 63 of them in the
            # solve mask: skymsk/simplemsk do NOT save it (measured on real exposures)
            withenv("ARM_PRIOR_SUPPORT_RELTOL" => "0") do
                pv_off = load_fiber_priors(pd_real, adj_lco)
                @test pv_off[1] == lco_raw            # guard disabled -> raw mask
                @test count(pv_off[4][DEAD]) == 63    # skymsk keeps 63 dead pixels
            end

            # threshold robustness: the population is bimodal with a wide empty band,
            # so anything in ~3e-3..0.19 gives the identical mask
            for rt in ("3e-3", "5e-2", "1.5e-1")
                withenv("ARM_PRIOR_SUPPORT_RELTOL" => rt) do
                    @test findall(lco_raw .& .!load_fiber_priors(pd_real, adj_lco)[1]) == DEAD
                end
            end
        else
            @info "priors: skipping LCO prior-support guard test (lco priors not available)"
        end
    end
else
    @info "priors: skipping real-file contract tests (prior store not available)"
end
