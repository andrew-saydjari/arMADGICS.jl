# E5: end-to-end test of the COMBINED bright-mask policy (AKS 2026-09-07) on REAL cached
# median_sky, through the SAME code path the builder uses (resolve_bright_mask +
# check_bright_fraction from build_sky_defs.jl, and e5_parse_thresh_policy from
# e5_sky_pass1_defs.jl).
#
# The property that matters, and that this asserts: the delivered mask is
# FIBER-INDEPENDENT. For every fiber, the bright pixels must be exactly
# combined_mask ∩ that fiber's submsk — no fiber may contribute a pixel of its own.
# Usage: julia --project=<arM-E5b> e5_combined_policy_test.jl
using HDF5, Statistics, Printf

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")   # dilate_msk, running_spread_fast

isnanorzeroorinf(x) = isnan(x) || iszero(x) || isinf(x)
function expand_msk(msk; rad=1)
    lmsk = length(msk); msknew = ones(Bool, lmsk)
    for i in 1:lmsk
        if any(.!view(msk, maximum((1, i - rad)):minimum((i + rad, lmsk)))); msknew[i] = false; end
    end
    msknew
end
# pull in just the policy block of build_sky_defs.jl (as e5_policy_guard_test.jl does)
src = read(ARM * "scripts/prior_build/build_sky_defs.jl", String)
i1 = first(findfirst("# bright/faint sky-line split threshold", src))
i2 = first(findfirst("\"\"\"\n    build_skyLines(", src))
include_string(Main, src[i1:i2-1], "policy_defs.jl")
# and the parser
src2 = read(ARM * "scripts/prior_build/e5_sky_pass1_defs.jl", String)
j1 = first(findfirst("const E5_BRIGHT_COMBINED_DEFAULT", src2))
j2 = first(findfirst("\n\"\"\"\n    e5_link_skycont", src2))
include_string(Main, src2[j1:j2], "policy_parse.jl")

const SC = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/screens"
const TESTFIB = [1, 10, 76, 150, 295, 300, 301, 388, 450, 519, 570, 600]

nfail = 0
for spec in ["combined", "combined:union", "combined:majority", "combined:apo_drop12"]
    pol, tag = e5_parse_thresh_policy(spec)
    ref = Bool.(h5read(pol[:path], "mask_" * pol[:variant]))
    @printf("\n=== %-24s tag=%-22s (%d px in the combined mask)\n", spec, tag, count(ref))
    for f in TESTFIB
        n = lpad(f, 3, "0")
        ms = h5read(joinpath(SC, "median_sky_$n.h5"), "median_sky")
        sub = Bool.(h5read(joinpath(SC, "median_sky_$n.h5"), "submsk"))
        mskflux, desc, thr = resolve_bright_mask(pol, f, ms, sub)
        frac = count(mskflux) / length(mskflux)
        # THE fiber-independence assertion
        expect = ref[sub]
        ok = mskflux == expect
        ok || (global nfail += 1)
        threw = false
        try
            check_bright_fraction(frac, f, thr, desc; bounds=(0.04, 0.15), guard=:error)
        catch; threw = true; end
        @printf("  f%s (%s): bright %6.3f%%  fiber-independent=%-5s  guard(:error) %s\n",
            n, f <= 300 ? "apo" : "lco", 100frac, ok,
            threw ? "FIRED (build would STOP)" : "passed")
    end
end

# a bad variant and a bad spec must fail LOUDLY at parse time, not mid-build
for bad in ["combined:nope", "combined:"]
    caught = false
    try; e5_parse_thresh_policy(bad); catch e; caught = true; end
    @printf("\nparse %-18s -> %s\n", repr(bad), caught ? "REJECTED at parse time (good)" : "ACCEPTED (BAD!)")
    caught || (global nfail += 1)
end

println(nfail == 0 ? "\nALL PASS: the combined mask is fiber-independent and the guard is live." :
        "\n$nfail FAILURE(S)")
exit(nfail == 0 ? 0 : 1)
