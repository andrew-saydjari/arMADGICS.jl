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
default_thresh_bright_faint(adjfibindx) = (adjfibindx <= 300) ? 2000 : 645

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
        min_obscnt=10, thresh_bright_faint=default_thresh_bright_faint, silent=true)
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
        mskflux = .!expand_msk(median_sky .< thresh_bright_faint(adjfibindx),rad=4)
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
            end
        end
    end
    return fnameFaint, fnameFaintGSPICE
end
