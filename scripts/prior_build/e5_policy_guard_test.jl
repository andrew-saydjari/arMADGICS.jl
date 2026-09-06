# Verify the finding-#35 guard + policy on REAL cached median_sky (fibers 10 apo, 600 lco).
using HDF5, Statistics, Printf
ARM="/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
isnanorzeroorinf(x) = isnan(x) || iszero(x) || isinf(x)
nanzeromedian(x) = all(isnanorzeroorinf, x) ? NaN : median(filter(!isnanorzeroorinf, x))
nanzeromedian(x, y) = mapslices(nanzeromedian, x, dims=y)
function expand_msk(msk; rad=1)
    lmsk=length(msk); msknew=ones(Bool,lmsk)
    for i in 1:lmsk
        if any(.!view(msk, maximum((1,i-rad)):minimum((i+rad,lmsk)))); msknew[i]=false; end
    end
    msknew
end
src = read(ARM*"scripts/prior_build/build_sky_defs.jl", String)
i1 = first(findfirst("# bright/faint sky-line split threshold", src))
i2 = first(findfirst("\"\"\"\n    build_skyLines(", src))
include_string(Main, src[i1:i2-1], "policy_defs.jl")

SC="/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/screens"
for (f,tele) in [(10,"apo"),(600,"lco")]
    ms = h5read(joinpath(SC,"median_sky_$(lpad(f,3,"0")).h5"), "median_sky")
    npix = length(ms)
    for spec in ["legacy","off","abs","quantile"]
        pol = spec=="legacy" ? nothing :
              spec=="off" ? Dict(:mode=>:off) :
              spec=="abs" ? (tele=="apo" ? Dict(:mode=>:absolute,:apo=>35.0,:lco=>8.0) : Dict(:mode=>:absolute,:apo=>35.0,:lco=>8.0)) :
              Dict(:mode=>:quantile,:bright_frac=>0.083)
        thr, desc = resolve_bright_threshold(pol, f, ms)
        mskflux = isnothing(thr) ? falses(npix) : .!expand_msk(ms .< thr, rad=4)
        frac = count(mskflux)/npix
        # guard behaviour under :warn (never throws) and :error (throws when out of band)
        threw = false
        try
            check_bright_fraction(frac, f, thr, desc; guard=:error, split_off=isnothing(thr))
        catch; threw = true; end
        @printf("f%03d %-8s thr=%-8s bright=%6.3f%%  guard(:error) %s\n", f, spec,
            isnothing(thr) ? "none" : string(round(thr,digits=2)), 100frac,
            threw ? "FIRED (build would stop)" : "passed")
    end
end
