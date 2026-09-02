# Run the arM ingest + core per-spectrum pipeline path on the M1/M2/M3 fixture
# (built by build_M123_fixture.jl) and record per-fiber outcomes, so the same
# driver documents behavior BEFORE and AFTER the M1/M2/M3 fixes.
#
# Usage:
#   julia --project=<env> scripts/validation/run_M123_fixture.jl <src_dir> <fixture_dir> <out_h5>
# where <src_dir> is the arMADGICS checkout whose src/ should be exercised
# (allows running the pre-fix and post-fix code with the identical driver),
# <fixture_dir> is the fixture reduxBase, and <out_h5> collects the results.
#
# Priors (read-only): the production 2025_07_31 set from pipeline.jl.
# Almanac (read-only): the 2026_05_01 run's allobs file.

using LinearAlgebra, Statistics, StatsBase, HDF5, JLD2, FITSIO
using EllipsisNotation, ShiftedArrays, SparseArrays, Interpolations, LowRankOps
BLAS.set_num_threads(4)

length(ARGS) >= 3 || error("usage: run_M123_fixture.jl <src_dir> <fixture_dir> <out_h5>")
src_dir = abspath(ARGS[1]) * "/"
fixture_dir = abspath(ARGS[2])
out_h5 = abspath(ARGS[3])

include(joinpath(src_dir, "src/utils.jl"))
include(joinpath(src_dir, "src/gridSearch.jl"))
include(joinpath(src_dir, "src/componentAndPosteriors.jl"))
include(joinpath(src_dir, "src/fileNameHandling.jl"))
include(joinpath(src_dir, "src/ingest.jl"))
include(joinpath(src_dir, "src/lowRankPrescription.jl"))
include(joinpath(src_dir, "src/marginalizeEW.jl"))
include(joinpath(src_dir, "src/spectraInterpolation.jl"))
include(joinpath(src_dir, "src/chi2Wrappers.jl"))
include(joinpath(src_dir, "src/pipelineCore.jl"))

## globals expected by pipeline_single_spectra (mirrors pipeline.jl)
prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
reduxBase = fixture_dir
almanacFile = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/almanac/allobs_57600_61160.h5"
cache_dir = ""
inject_cache_dir = ""
refine_iters = 5
minRVres = 1 // 10
RV_err_step = 4
lvl1 = -70:1//2:70
lvl2 = -8:2//10:8
lvl3 = -3:1//10:3
slvl_tuple = (lvl1, lvl2, lvl3)
logUniWaveAPOGEE = 10 .^ range((start = 4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)

prior_dict = Dict{String,String}()
prior_dict["chebmsk"] = prior_dir * "2025_07_31/prior_dump/chebmsk_exp.h5"
prior_dict["starCont"] = prior_dir * "2025_07_31/prior_dump/APOGEE_starcont_svd_60_rough.h5"
prior_dict["starLines_refLSF"] = prior_dir * "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5"

## load priors exactly as multi_spectra_batch does
f = h5open(prior_dict["starLines_refLSF"]); V_starlines_refLSF = read(f["Vmat"]); close(f)
f = h5open(prior_dict["chebmsk"]); chebmsk_exp = read(f["chebmsk_exp"]); close(f)
skymsk_bright = chebmsk_exp
skymsk_faint = chebmsk_exp
skymsk = chebmsk_exp
f = h5open(prior_dict["starCont"]); V_starcont = read(f["Vmat"]); close(f)
V_starlines = V_starlines_refLSF # hack (matches production)
msk_starCor = ones(Bool, length(chebmsk_exp))
prior_vec = (chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont, V_starlines_refLSF, V_starlines, msk_starCor)

## fibers to run (see build_M123_fixture.jl)
run_fibers = [85, 86, 87, 88, 89, 90, 92, 93, 100, 101, 103]
fiber_notes = Dict(90 => "dead fiber_id 211", 100 => "injected all-masked", 101 => "injected NaN-flux/good-mask (A4)", 103 => "injected partial-NaN flux (good msk)")

tele, mjd, expnum = "apo", "58588", 11

## independent, code-version-agnostic diagnostics straight from the fixture file
raw = h5open(get_1Duni_name(reduxBase, tele, mjd, expnum)) do fh
    (flux=read(fh["flux_1d"]), ivar=read(fh["ivar_1d"]), msk=read(fh["mask_1d"]))
end

nfib = length(run_fibers)
status = fill("", nfib)
snr_pipe = fill(NaN, nfib)
snr_oldform = fill(NaN, nfib)
snr_newform = fill(NaN, nfib)
starscale0_v = fill(NaN, nfib)
ingestBit_v = fill(-999, nfib)
skyBit_v = fill(-999, nfib)
nSkyFibers_v = fill(-999, nfib)
rv_pixoff = fill(NaN, nfib)
rv_flag = fill(-999, nfib)
chi2res_v = fill(NaN, nfib)
final_pix_cnt = fill(-1, nfib)
apVisit_absmax = fill(NaN, nfib)
apVisit_nfloored = fill(-1, nfib)
errmsg = fill("", nfib)

med_or_nan(x) = isempty(x) ? NaN : median(x)

outs = Vector{Any}(nothing, nfib) # keep full outputs for the batch-extraction contract check

for (i, fib) in enumerate(run_fibers)
    fs = raw.flux[:, fib]; iv = raw.ivar[:, fib]; mk = Bool.(raw.msk[:, fib])
    # old (inverted, unmasked) and new (masked flux*sqrt(ivar)) snr formulas, computed
    # directly from the file so both numbers exist regardless of code version
    old_vals = filter(x -> !isnan(x) && !iszero(x), fs ./ sqrt.(iv))
    snr_oldform[i] = med_or_nan(old_vals) # NOTE: +/-Inf deliberately NOT filtered (old behavior)
    new_vals = filter(x -> isfinite(x) && !iszero(x), (fs .* sqrt.(iv))[mk])
    snr_newform[i] = med_or_nan(new_vals)

    argtup = (tele=tele, mjd=mjd, expnum=expnum, adjfiberindx=fib)
    t0 = time()
    local out
    try
        out = pipeline_single_spectra(argtup, prior_vec)
        outs[i] = out
        status[i] = "ok"
    catch e
        status[i] = "CRASH"
        errmsg[i] = sprint(showerror, e)
        println("fiber $fib CRASHED after $(round(time()-t0, digits=1))s: ", errmsg[i][1:min(end, 300)])
        continue
    end
    meta = out[1]
    snr_pipe[i] = meta[10]
    starscale0_v[i] = meta[2]
    ingestBit_v[i] = length(meta) >= 11 ? meta[11] : -999
    skyBit_v[i] = length(meta) >= 12 ? meta[12] : -999
    nSkyFibers_v[i] = Int(meta[9])
    rv = out[2][1]
    rv_pixoff[i] = Float64(rv[1])
    rv_flag[i] = Int(rv[6])
    chi2res_v[i] = out[3][1]
    final_pix_cnt[i] = Int(out[3][3])
    apV = out[4][9]; fmsk_out = out[4][10]
    finite_apV = filter(isfinite, apV)
    apVisit_absmax[i] = isempty(finite_apV) ? NaN : maximum(abs.(finite_apV))
    apVisit_nfloored[i] = count(Bool.(fmsk_out) .& .!isfinite.(apV))
    println("fiber $fib done in $(round(time()-t0, digits=1))s: snr=$(round(snr_pipe[i], sigdigits=4)) ingestBit=$(ingestBit_v[i]) rv=$(round(rv_pixoff[i], sigdigits=4)) flag=$(rv_flag[i])")
    flush(stdout)
end

rm(out_h5; force=true)
mkpath(dirname(out_h5))
h5open(out_h5, "w") do fh
    fh["fiberindx"] = run_fibers
    fh["status"] = status
    fh["snr_pipe"] = snr_pipe
    fh["snr_oldformula_unmasked"] = snr_oldform
    fh["snr_newformula_masked"] = snr_newform
    fh["starscale0"] = starscale0_v
    fh["ingestBit"] = ingestBit_v
    fh["skyBit"] = skyBit_v
    fh["nSkyFibers"] = nSkyFibers_v
    fh["RV_pixoff_final"] = rv_pixoff
    fh["RV_flag"] = rv_flag
    fh["chi2res"] = chi2res_v
    fh["final_pix_cnt"] = final_pix_cnt
    fh["apVisit_absmax"] = apVisit_absmax
    fh["apVisit_nfloored_in_finalmsk"] = apVisit_nfloored
    fh["errmsg"] = errmsg
    attrs(fh)["src_dir"] = src_dir
    attrs(fh)["fixture_dir"] = fixture_dir
end

println("\n=== summary (src_dir = $src_dir) ===")
println(rpad("fib", 5), rpad("note", 34), rpad("status", 7), rpad("snr_pipe", 13), rpad("snr_old_frm", 13), rpad("snr_new_frm", 13), rpad("starscale0", 13), rpad("ingestBit", 10), rpad("skyBit", 8), rpad("nSky", 6), rpad("RVpix", 10), rpad("apV|max|", 12), "err")
for (i, fib) in enumerate(run_fibers)
    println(rpad(fib, 5), rpad(get(fiber_notes, fib, "sci"), 34), rpad(status[i], 7),
        rpad(round(snr_pipe[i], sigdigits=4), 13), rpad(round(snr_oldform[i], sigdigits=4), 13),
        rpad(round(snr_newform[i], sigdigits=4), 13), rpad(round(starscale0_v[i], sigdigits=4), 13),
        rpad(ingestBit_v[i], 10), rpad(skyBit_v[i], 8), rpad(nSkyFibers_v[i], 6), rpad(round(rv_pixoff[i], sigdigits=4), 10),
        rpad(round(apVisit_absmax[i], sigdigits=4), 12), errmsg[i][1:min(end, 80)])
end
println("results written to $out_h5")

## Batch-extraction contract check (post-fix code only): run the multi_spectra_batch
## extraction lambdas (kept in sync with pipeline.jl RVextract) over the mixed
## success/failure outs to catch shape/type mismatches that would crash a real batch save.
if all(x -> x !== nothing, outs) && length(outs[1][1]) >= 12
    metai = 1
    RVind, RVchi, RVcom, strpo = 2, 3, 4, 5
    adjfiberindx = 0 # dummy for the contract check
    RVextract = [
        (x -> x[metai][1], "data_pix_cnt"),
        (x -> x[metai][2], "starscale"),
        (x -> x[metai][3], "skyscale0"),
        (x -> x[metai][4], "flux"),
        (x -> x[metai][5], "fluxivar"),
        (x -> x[metai][6], "flux_nans"),
        (x -> x[metai][7], "fluxerr2_nans"),
        (x -> convert(Vector{Int}, x[metai][8]), "simplemsk"),
        (x -> x[metai][9], "nSkyFibers"),
        (x -> x[metai][10], "snr"),
        (x -> x[metai][11], "ingestBit"),
        (x -> x[metai][12], "skyBit"),
        (x -> adjfiberindx, "adjfiberindx"),
        (x -> Float64.(x[RVind][1][1]), "RV_pixoff_final"),
        (x -> Float64.(x[RVind][1][3]), "RV_pixoff_disc_final"),
        (x -> x[RVind][1][2], "RV_minchi2_final"),
        (x -> x[RVind][1][6], "RV_flag"),
        (x -> x[RVind][1][7], "RV_pix_var"),
        (x -> x[RVchi][1], "RVchi2_residuals"),
        (x -> x[RVchi][2], "avg_flux_conservation"),
        (x -> x[RVchi][3], "final_pix_cnt"),
        (x -> x[RVchi][4], "starscale1"),
        (x -> x[RVchi][5], "skyscale1"),
        (x -> x[RVind][2][1][3], "RV_p5delchi2_lvl1"),
        (x -> x[RVind][2][2][3], "RV_p5delchi2_lvl2"),
        (x -> x[RVind][2][3][3], "RV_p5delchi2_lvl3"),
        (x -> x[RVcom][1], "x_residuals_z_v0"),
        (x -> x[RVcom][2], "x_residuals_v0"),
        (x -> x[RVcom][3], "x_skyLines_faint_v0"),
        (x -> x[RVcom][4], "x_skyContinuum_v0"),
        (x -> x[RVcom][5], "x_starContinuum_v0"),
        (x -> x[RVcom][6], "x_starLines_v0"),
        (x -> x[RVcom][7], "x_starLineCof_v0"),
        (x -> x[RVcom][8], "tot_p5chi2_v0"),
        (x -> x[RVcom][9], "apVisit_v0"),
        (x -> Int.(x[RVcom][10]), "finalmsk"),
        (x -> x[RVcom][11], "x_starLines_restFrame_v0"),
        (x -> x[strpo], "x_starLines_err_v0"),
    ]
    function extractor(x, elemap, elename, savename)
        len = length(x)
        exobj = elemap(x[1])
        outmat = zeros(eltype(exobj), size(exobj)..., len)
        for i = 1:len
            outmat[.., i] .= elemap(x[i])
        end
        h5write(savename, elename, outmat)
    end
    batchfile = replace(out_h5, ".h5" => "_batchstyle.h5")
    rm(batchfile; force=true)
    for elelst in RVextract
        extractor(outs, elelst[1], elelst[2], batchfile)
    end
    println("batch-style extraction contract check PASSED -> $batchfile")
else
    println("batch-style extraction contract check SKIPPED (pre-fix code or crashed fibers)")
end
