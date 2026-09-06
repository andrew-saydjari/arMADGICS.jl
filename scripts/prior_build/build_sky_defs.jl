## Core definitions for build_skyCont.jl and build_skyLines.jl (E5): turn the sky
## samples written by sample_sky.jl into the skyCont / skyLines SVD priors.
#
# E5 producer/consumer contract: sample_sky.jl (via sample_sky_defs.jl) writes
# per-fiber HDF5 files into a sample directory:
#   skyflux_NNN.h5  : "skyflux"  (nwave x nsamp, Float64)
#   skyivar_NNN.h5  : "skyivar"  (nwave x nsamp, Float64, INVERSE variance; 0 = masked)
#   skymsk_NNN.h5   : "skymsk"   (nwave x nsamp, Bool, true = good)
#   skycont_NNN.h5  : "skycont"  (nwave x nsamp, Float64, smooth continuum component)
#   skyline_NNN.h5  : "skyline"  (nwave x nsamp, Float64, line component, 0 off-mask)
# The builders below read exactly these files/datasets from a `sample_dir` argument
# (the .jdat-era builders deserialized Serialization files from hardcoded date paths,
# and build_skyLines expected a skyvar_* VARIANCE file that the sampler never wrote;
# see gspice_ivar_from_skyivar for the variance-vs-ivar determination).
#
# Requires (from the includer): HDF5, LinearAlgebra, plus scripts/prior_build/prior_utils.jl
# (nanzeromedian, expand_msk) and, for build_skyLines, scripts/prior_build/gspice.jl
# (module gspice, resolved at call time).
# Author - Andrew Saydjari

wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat

"""
    get_config_arg(flag; env_key=nothing, default=nothing, args=ARGS)

Tiny CLI/env config reader for the prior-build drivers: returns the value of
`--flag=value` from `args`, else `ENV[env_key]`, else `default`; errors when the
argument is required (no default) and absent. Used so the sample directory is an
explicit input instead of a hardcoded date path (E5).
"""
function get_config_arg(flag::String; env_key=nothing, default=nothing, args=ARGS)
    pref = flag * "="
    for a in args
        if startswith(a, pref)
            return String(a[nextind(a, length(pref)):end])
        end
    end
    if !isnothing(env_key) && haskey(ENV, env_key)
        return ENV[env_key]
    end
    if isnothing(default)
        error("Required argument $flag=<value> (or ENV[\"$env_key\"]) not provided")
    end
    return default
end

sky_sample_path(sample_dir, prefix, adjfibindx) =
    joinpath(sample_dir, prefix * "_" * lpad(adjfibindx, 3, "0") * ".h5")

"""
    read_sky_sample(sample_dir, prefix, adjfibindx)

Read one per-fiber sample matrix written by sample_sky (dataset `prefix` in
`sample_dir/prefix_NNN.h5`) and check the wavelength dimension.
"""
function read_sky_sample(sample_dir, prefix, adjfibindx)
    path = sky_sample_path(sample_dir, prefix, adjfibindx)
    isfile(path) || error("read_sky_sample: $path not found; expected the $(prefix)_NNN.h5 files written by sample_sky.jl in sample_dir=$sample_dir")
    out = h5read(path, prefix)
    (size(out, 1) == length(wavetarg)) || error("read_sky_sample: $path dataset $prefix has leading dimension $(size(out,1)), expected length(wavetarg)=$(length(wavetarg))")
    return out
end

function read_chipgap_msk(chipgap_msk_path, adjfibindx)
    return h5open(chipgap_msk_path, "r") do f
        if adjfibindx > 300
            read(f["lco"])
        else
            read(f["apo"])
        end
    end
end

"""
    gspice_ivar_from_skyivar(skyivar_sub)

E5 SEMANTIC CHECK (variance vs inverse variance at the gspice boundary):

`gspice.gspice_covar_iter_mask(flux, ivar, mask; ...)` consumes INVERSE variance —
inside `gspice_standard_scale` it computes `meanivar = sum(ivar.*wt,...)` and scales
each spectrum by `refscale = sqrt.(meanivar)`, i.e. into S/N units, which is only
correct for ivar. The .jdat-era `skyvar_*` files held VARIANCE, hence the historical
`ivar = 1 ./ fluxvar` in build_skyLines.jl. The sampler (sample_sky_defs.jl) writes
`skyivar_*` = the redux `ivar_1d` columns (already inverse variance; it is used as
`Ainv = Diagonal(fivar)` there), so it passes through UNINVERTED here. Masked pixels
carry skyivar = 0, which is valid zero-weight gspice input (the old path produced the
same 0 from 1/Inf); inverting skyivar instead would inject Infs at those pixels.

Input is the (npix_sub x nsamp_sub) slice of the sampler's skyivar matrix; output is
the (nsamp_sub x npix_sub) Float64 ivar array in gspice's (nspec, npix) convention.
"""
function gspice_ivar_from_skyivar(skyivar_sub::AbstractMatrix)
    ivar = convert.(Float64, collect(skyivar_sub'))
    all(isfinite, ivar) || error("gspice_ivar_from_skyivar: non-finite skyivar entries; the sampler should write finite inverse variances (0 for masked pixels)")
    all(>=(0), ivar) || error("gspice_ivar_from_skyivar: negative skyivar entries")
    return ivar
end

"""
    build_skyCont(adjfibindx; sample_dir, chipgap_msk_path, out_dir="sky_priors/", nsub=30)

Build the sky-continuum SVD prior for one (adjusted) fiber index from the
`skycont_NNN.h5` samples in `sample_dir`. Writes
`out_dir/APOGEE_skycont_svd_<nsub>_fNNN.h5` (checkpointed: skipped if it exists)
and returns that filename. Math unchanged from the .jdat-era builder; only the
sample I/O contract changed (E5). The .jdat-era code also loaded skymsk_* here but
never used it (the skycont samples are dense smooth fits, so no per-sample pixel
mask enters the continuum covariance) — that dead read is dropped.
"""
function build_skyCont(adjfibindx; sample_dir, chipgap_msk_path,
        out_dir="sky_priors/", nsub=30)
    fname = joinpath(out_dir, "APOGEE_skycont_svd_"*string(nsub)*"_f"*lpad(adjfibindx,3,"0")*".h5")
    if !isfile(fname)
        skycont = read_sky_sample(sample_dir, "skycont", adjfibindx)
        chipgapmsk = read_chipgap_msk(chipgap_msk_path, adjfibindx)

        specsum = dropdims(sum(skycont,dims=1),dims=1)
        Vred = skycont[chipgapmsk,specsum.>0];
        # weights = ones(size(Vred,2));
        # Vred .*= reshape(weights,1,:);
        nsamp = size(Vred,2)
        # norm_weights = weights'*weights
        Csky = Vred*Vred'
        # Csky./=norm_weights
        Csky./=nsamp

        SF = svd(Csky);
        EVEC = zeros(length(wavetarg),size(SF.U,2))
        EVEC[chipgapmsk,:].=SF.U;

        dirName = splitdir(fname)[1]
        if !ispath(dirName)
            mkpath(dirName)
        end
        h5write(fname,"Vmat",EVEC[:,1:nsub]*Diagonal(sqrt.(SF.S[1:nsub])))
        h5write(fname,"λv",SF.S[1:nsub])
    end
    return fname
end

# bright/faint sky-line split threshold: APO fibers (adjfibindx<=300) vs LCO (301-600)
# INHERITED FROM apMADGICS (src/prior_build/build_skyLines.jl:56-62), where the samples
# were DRP apCframe fluxes. arM samples ar1Dunical, whose flux scale is ~58x (APO) /
# ~81x (LCO) smaller, so these constants flag ZERO pixels here — a silent no-op.
# See prior_outputs/sky_pass1/THRESHOLD_FINDING.md (finding #35). Do not "fix" the
# numbers without a policy decision; use bright_policy= (below), which is guarded.
default_thresh_bright_faint(adjfibindx) = (adjfibindx <= 300) ? 2000 : 645

"""
    resolve_bright_threshold(policy, adjfibindx, median_sky, legacy_fun) -> (thresh, desc)

Resolve the bright/faint split threshold under a threshold POLICY, so the choice is a
parameter rather than a hardcoded constant. `policy` is `nothing` (legacy: use
`legacy_fun(adjfibindx)`, i.e. the inherited 2000/645) or a Dict with `:mode`:

  - `:off`       — no split at all; returns `nothing`. Faint prior spans every pixel
                   (this is exactly what the inherited constants do today, but declared).
  - `:absolute`  — `policy[:apo]` / `policy[:lco]`, in the sample's own flux units.
  - `:quantile`  — `policy[:bright_frac]`: solve (bisection on the POST-`expand_msk`
                   fraction, so the target is the realized bright fraction) for the
                   threshold flagging that fraction of pixels. Unit-free: immune to any
                   future change of flux normalization.

Returns the threshold (or `nothing` for `:off`) and a human-readable description that
is recorded in the built prior for provenance.
"""
function resolve_bright_threshold(policy, adjfibindx, median_sky, legacy_fun=default_thresh_bright_faint)
    if isnothing(policy)
        t = legacy_fun(adjfibindx)
        return float(t), "legacy inherited absolute (apo 2000 / lco 645; see THRESHOLD_FINDING.md)"
    end
    mode = get(policy, :mode, :absolute)
    if mode == :off
        return nothing, "split off (single unified sky-line prior over all pixels)"
    elseif mode == :absolute
        v = (adjfibindx <= 300) ? policy[:apo] : policy[:lco]
        return float(v), "absolute apo=$(policy[:apo]) lco=$(policy[:lco])"
    elseif mode == :quantile
        target = float(policy[:bright_frac])
        (0 < target < 1) || error("resolve_bright_threshold: :bright_frac must be in (0,1), got $target")
        vals = filter(isfinite, median_sky)
        isempty(vals) && error("resolve_bright_threshold: median_sky has no finite entries for adjfibindx=$adjfibindx")
        npix = length(median_sky)
        lo, hi = max(minimum(vals), 1e-6), max(maximum(vals), 1e-3)
        thr = hi
        for _ in 1:60
            mid = sqrt(lo * hi)
            fr = count(.!expand_msk(median_sky .< mid, rad=4)) / npix
            fr > target ? (lo = mid) : (hi = mid)
            thr = sqrt(lo * hi)
        end
        return thr, "quantile target_bright_frac=$(target) -> thresh=$(round(thr, digits=4))"
    end
    error("resolve_bright_threshold: unknown policy mode $(repr(mode)) (expected :off, :absolute, :quantile)")
end

"""
    resolve_bright_mask(policy, adjfibindx, median_sky, submsk, legacy_fun) -> (mskflux, desc, thresh)

Bright-pixel mask in the builder's COMPRESSED (submsk) space. Handles the scalar-threshold
policies via `resolve_bright_threshold`, plus the calibrated line DETECTOR
(`:linedetect`, AKS 2026-09-06), which cannot be expressed as any single flux value.

The detector runs on the FULL wavetarg grid (median_sky scattered back through `submsk`,
NaN elsewhere) because a moving window over the compressed space would splice chip gaps
together. Steps, and what the calibration against DR17's actual bright mask said about
each (see prior_outputs/sky_pass1/BRIGHT_MASK_CALIBRATION.md):
  1. local continuum: `cont_window` (DEFAULT 0 = OFF). Measured to be unnecessary — the
     samples are already continuum-subtracted upstream (sample_sky_defs.jl:69), so the
     baseline is ~0 and subtracting a local median changes pixel IoU by <0.01.
  2. locally-adaptive SCALE: running MAD over `scale_window` px, then flag z = resid/scale
     > `k`. This is where the throughput adaptation actually happens (AKS's concern that we
     have not throughput-divided): a GLOBAL scale under-flags the red end, where an
     undivided throughput suppresses line amplitudes — red-end line recall 0.73 (global)
     vs 0.92 (running, APO).
  3. DILATION by `dilation` px to cover line wings. The calibration independently selects
     4, which is exactly DR17's own `expand_msk(..., rad=4)`.

Requires scripts/prior_build/e5_bright_line_detect.jl to be included by the caller.
"""
function resolve_bright_mask(policy, adjfibindx, median_sky, submsk,
        legacy_fun=default_thresh_bright_faint)
    if !isnothing(policy) && get(policy, :mode, :absolute) == :linedetect
        sw = get(policy, :scale_window, 2001)
        k = float(get(policy, :k, 90.0))
        dil = get(policy, :dilation, 4)
        cw = get(policy, :cont_window, 0)
        x = fill(NaN, length(submsk))
        x[submsk] .= median_sky
        cont = cw == 0 ? zeros(length(x)) : running_median_nan(x, cw)
        resid = x .- cont
        scale = running_spread_fast(resid, sw; kind=:mad, stride=50)
        base = falses(length(x))
        @inbounds for i in eachindex(x)
            if isfinite(resid[i]) && isfinite(scale[i]) && scale[i] > 0 && resid[i] > k * scale[i]
                base[i] = true
            end
        end
        full = dilate_msk(base, dil)
        @inbounds for i in eachindex(full)
            isfinite(x[i]) || (full[i] = false)
        end
        desc = "linedetect scale_window=$sw k=$k dilation=$dil cont_window=$cw (calibrated vs DR17 bright mask)"
        return full[submsk], desc, NaN
    end
    thr, desc = resolve_bright_threshold(policy, adjfibindx, median_sky, legacy_fun)
    mskflux = if isnothing(thr)
        falses(length(median_sky))
    else
        .!expand_msk(median_sky .< thr, rad=4)
    end
    return mskflux, desc, isnothing(thr) ? nothing : thr
end

"""
    check_bright_fraction(frac, adjfibindx, thresh, desc; bounds, guard, split_off)

Sanity guard on the realized bright fraction of the bright/faint split (finding #35:
the inherited thresholds silently flagged 0% of pixels in ar1Dunical units, while the
deployed DR17-era priors split ~8.3%). Fires when `frac` falls outside `bounds`
(default 1%-20%); `guard` is `:error`, `:warn` (default), or `:off`. Skipped when
`split_off` is true, where a 0% fraction is the declared intent (and is asserted).
"""
function check_bright_fraction(frac, adjfibindx, thresh, desc;
        bounds=(0.01, 0.20), guard::Symbol=:warn, split_off::Bool=false)
    if split_off
        iszero(frac) || error("check_bright_fraction: split declared OFF for adjfibindx=$adjfibindx but $(100frac)% of pixels were flagged bright")
        return frac
    end
    lo, hi = bounds
    if !(lo <= frac <= hi)
        msg = "build_skyLines bright/faint split flagged $(round(100frac, digits=4))% of pixels " *
              "for adjfibindx=$adjfibindx (threshold=$thresh, policy: $desc), outside the expected " *
              "[$(100lo)%, $(100hi)%] band. 0% means the threshold is a SILENT NO-OP in these flux " *
              "units — the failure mode of finding #35 (THRESHOLD_FINDING.md): the inherited " *
              "apo 2000 / lco 645 flag nothing on ar1Dunical samples, where the deployed DR17-era " *
              "priors split ~8.3%. Set bright_policy=Dict(:mode=>:off) to declare a no-split build."
        if guard == :error
            error(msg)
        elseif guard == :warn
            @warn msg
            println("WARNING: " * msg); flush(stdout)
        end
    end
    return frac
end

"""
    build_skyLines(adjfibindx; sample_dir, chipgap_msk_path, out_dir="sky_priors/",
                   nsub_faint=120, nsub_bright=120, nsigma_schedule=[20, 8, 6],
                   usamp_factor=7, maxbadpix=650, reg_eps=1e-3, min_obscnt=10,
                   thresh_bright_faint=default_thresh_bright_faint, silent=true)

Build the faint sky-line SVD priors (plain + GSPICE-masked) for one (adjusted)
fiber index from the `skyline_NNN.h5` / `skymsk_NNN.h5` / `skyivar_NNN.h5` samples
in `sample_dir`. Writes `out_dir/APOGEE_skyline_faint_svd_<nsub_faint>_fNNN.h5` and
`out_dir/APOGEE_skyline_faint_GSPICE_svd_<nsub_bright>_fNNN.h5` (checkpointed) and
returns the two filenames. Bright priors are still not produced (unrecoverable, per
the original in-file comment); their names only gate the outer checkpoint check.
Math unchanged from the .jdat-era builder except the skyivar wiring — see
gspice_ivar_from_skyivar (the old code read a VARIANCE skyvar_* file and inverted
it; the sampler writes inverse variance, used directly). Requires gspice.jl to be
included by the caller.
"""
function build_skyLines(adjfibindx; sample_dir, chipgap_msk_path,
        out_dir="sky_priors/", nsub_faint=120, nsub_bright=120,
        nsigma_schedule=[20, 8, 6], usamp_factor=7, maxbadpix=650, reg_eps=1e-3,
        min_obscnt=10, thresh_bright_faint=default_thresh_bright_faint, silent=true,
        bright_policy=nothing, bright_frac_bounds=(0.01, 0.20), bright_guard::Symbol=:warn)
    fnameFaint = joinpath(out_dir, "APOGEE_skyline_faint_svd_"*string(nsub_faint)*"_f"*lpad(adjfibindx,3,"0")*".h5")

    fnameFaintGSPICE = joinpath(out_dir, "APOGEE_skyline_faint_GSPICE_svd_"*string(nsub_bright)*"_f"*lpad(adjfibindx,3,"0")*".h5")

    # M4: fnameBright/fnameBrightGSPICE were referenced but never defined (UndefVarError);
    # bright priors are not currently produced (see comment below), so these names only
    # gate the outer checkpoint check
    fnameBright = joinpath(out_dir, "APOGEE_skyline_bright_svd_"*string(nsub_bright)*"_f"*lpad(adjfibindx,3,"0")*".h5")
    fnameBrightGSPICE = joinpath(out_dir, "APOGEE_skyline_bright_GSPICE_svd_"*string(nsub_bright)*"_f"*lpad(adjfibindx,3,"0")*".h5")

    if !(isfile(fnameBright) & isfile(fnameBrightGSPICE) & isfile(fnameFaint) & isfile(fnameFaintGSPICE))
        skyline = read_sky_sample(sample_dir, "skyline", adjfibindx)
        # E5: skymsk round-trips as Bool from the sampler's HDF5; keep the Float64
        # weighting semantics of the .jdat-era builder for the masked matmuls below
        skymsk = convert.(Float64, read_sky_sample(sample_dir, "skymsk", adjfibindx))
        skyivar = read_sky_sample(sample_dir, "skyivar", adjfibindx)

        chipgapmsk = read_chipgap_msk(chipgap_msk_path, adjfibindx)

        mkpath(out_dir) # h5write does not create the output directory

        # Sep Bright/Faint
        specsum = dropdims(sum(skyline,dims=1),dims=1)
        obscnt = dropdims(sum(skymsk,dims=2),dims=2);
        submsk = (obscnt.>=min_obscnt) .& chipgapmsk;
        Vred = skyline[submsk,specsum.>0];
        skymsked = skymsk[submsk,specsum.>0]
        Vred .*= skymsked;

        median_sky = dropdims(nanzeromedian(Vred,2),dims=2);

        submsk_bright = copy(submsk)
        submsk_faint = copy(submsk)
        # E5 finding #35: the split threshold is now a guarded POLICY, not a bare
        # constant. Default (bright_policy=nothing) reproduces the inherited
        # thresh_bright_faint behaviour exactly; the guard reports when that behaviour
        # is a no-op instead of failing silently.
        mskflux, bright_desc, bright_thresh = resolve_bright_mask(
            bright_policy, adjfibindx, median_sky, submsk, thresh_bright_faint)
        bright_frac = count(mskflux) / length(mskflux)
        split_off = !isnothing(bright_policy) && get(bright_policy, :mode, :absolute) == :off
        check_bright_fraction(bright_frac, adjfibindx, bright_thresh, bright_desc;
            bounds=bright_frac_bounds, guard=bright_guard, split_off=split_off)
        mskflux_big = zeros(Bool,length(submsk))
        mskflux_big[submsk].=mskflux

        submsk_bright[mskflux_big].&= true
        submsk_bright[.!mskflux_big].&= false

        submsk_faint[mskflux_big].&= false
        submsk_faint[.!mskflux_big].&= true;

        # not making bright priors for now, unrecoverable
        # Faint
        if !(isfile(fnameFaint) & isfile(fnameFaintGSPICE))
            specsum = dropdims(sum(skyline,dims=1),dims=1)
            obscnt = dropdims(sum(skymsk,dims=2),dims=2);
            submsk = (obscnt.>=min_obscnt) .& chipgapmsk .& submsk_faint;
            Vred = skyline[submsk,specsum.>0];
            skymsked = skymsk[submsk,specsum.>0];
            Vred .*= skymsked
            if !isfile(fnameFaint)
                norm_weights = skymsked*skymsked';
                Csky = Vred*Vred'
                Csky./=(norm_weights .+ (norm_weights.==0));

                SF = svd(Csky);
                EVEC = zeros(length(wavetarg),size(SF.U,2))
                EVEC[submsk,:].=SF.U;

                h5write(fnameFaint,"Vmat",EVEC[:,1:nsub_faint]*Diagonal(sqrt.(SF.S[1:nsub_faint])))
                h5write(fnameFaint,"λv",SF.S[1:nsub_faint])
                h5write(fnameFaint,"submsk",convert.(Int,submsk)) # different for bright/faint skylines
                # E5 finding #35 provenance: what the split actually did on this fiber
                h5write(fnameFaint,"bright_thresh",isnothing(bright_thresh) ? NaN : float(bright_thresh))
                h5write(fnameFaint,"bright_frac",bright_frac)
                h5write(fnameFaint,"bright_policy",bright_desc)
            end

            # Faint GSPICE
            if !isfile(fnameFaintGSPICE)
                flux = collect(Vred');
                # E5: sampler writes INVERSE variance; gspice consumes ivar, so no
                # inversion (the .jdat-era skyvar_* file held variance and was inverted
                # here). See gspice_ivar_from_skyivar docstring for the determination.
                ivar = gspice_ivar_from_skyivar(skyivar[submsk,specsum.>0]);
                mask = collect((skymsked.==0)'); #invert to 0/1 encoding of GSPICE

                out = gspice.gspice_covar_iter_mask(flux, ivar, mask; nsigma=nsigma_schedule, maxbadpix=maxbadpix, usamp_factor=usamp_factor, reg_eps=reg_eps, silent=silent);

                Vred_1 = copy(Vred)
                skymsked_1 = convert.(Float64,.!out[2]')
                Vred_1 .*= skymsked_1;
                norm_weights_1 = skymsked_1*skymsked_1';
                Csky = Vred_1*Vred_1'
                Csky./=(norm_weights_1 .+ (norm_weights_1.==0));

                SF = svd(Csky);
                EVEC = zeros(length(wavetarg),size(SF.U,2))
                EVEC[submsk,:].=SF.U;

                h5write(fnameFaintGSPICE,"Vmat",EVEC[:,1:nsub_bright]*Diagonal(sqrt.(SF.S[1:nsub_bright])))
                h5write(fnameFaintGSPICE,"λv",SF.S[1:nsub_bright])
                h5write(fnameFaintGSPICE,"submsk",convert.(Int,submsk)) # different for bright/faint skylines
                h5write(fnameFaintGSPICE,"bright_thresh",isnothing(bright_thresh) ? NaN : float(bright_thresh))
                h5write(fnameFaintGSPICE,"bright_frac",bright_frac)
                h5write(fnameFaintGSPICE,"bright_policy",bright_desc)
            end
        end
    end
    return fnameFaint, fnameFaintGSPICE
end
