## E5: CHECKED ASSERTION that skyCont products are invariant to the bright/faint
## threshold policy, i.e. reusable across every option in THRESHOLD_FINDING.md
## (finding #35) and never worth rebuilding.
#
# Two independent checks:
#  (1) STATIC: build_skyCont's source references no threshold/policy symbol, and its
#      method signature exposes no threshold knob. The split is computed inside
#      build_skyLines only.
#  (2) EMPIRICAL: rebuild skyCont from the same samples into a scratch dir and assert
#      the product is BITWISE identical (Vmat and λv) to the shipped one in built/.
#      Determinism + threshold-independence => the existing files are reusable as-is.
#
# Usage: julia --project=<arM-E5b> e5_assert_skycont_invariant.jl [adjfib ...]
#        (default: one APO and one LCO fiber that are already built)
using HDF5, LinearAlgebra, Printf

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
const OUT = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
const CHIPGAP = get(ENV, "E5_CHIPGAP", "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5")

include(ARM * "scripts/prior_build/build_sky_defs.jl")

failures = String[]

## (1) static check
src = read(ARM * "scripts/prior_build/build_sky_defs.jl", String)
body_start = findfirst("function build_skyCont(", src)
body_stop = findfirst("# bright/faint sky-line split threshold", src)
body = src[first(body_start):first(body_stop)]
for sym in ["thresh_bright_faint", "bright_policy", "bright_thresh", "median_sky",
    "submsk_faint", "submsk_bright", "resolve_bright_threshold"]
    if occursin(sym, body)
        push!(failures, "STATIC: build_skyCont body references threshold symbol `$sym`")
    end
end
kws = Base.kwarg_decl(first(methods(build_skyCont)))
for k in kws
    if occursin("thresh", String(k)) || occursin("bright", String(k)) || occursin("policy", String(k))
        push!(failures, "STATIC: build_skyCont exposes threshold-like kwarg `$k`")
    end
end
@printf("(1) STATIC: build_skyCont kwargs = %s -> %s\n", kws,
    isempty(failures) ? "no threshold dependence" : "FAILED")

## (2) empirical bitwise rebuild check
fibers = isempty(ARGS) ? Int[] : parse.(Int, ARGS)
if isempty(fibers)
    for f in vcat(1:300, 301:600)
        isfile(joinpath(OUT, "built", "APOGEE_skycont_svd_30_f" * lpad(f, 3, "0") * ".h5")) || continue
        isfile(joinpath(OUT, "samples", "skycont_" * lpad(f, 3, "0") * ".h5")) || continue
        push!(fibers, f)
        (any(<=(300), fibers) && any(>(300), fibers)) && break
        length(fibers) >= 2 && all(<=(300), fibers) && f > 300 && break
    end
    fibers = unique(vcat(first(filter(<=(300), fibers), 1), first(filter(>(300), fibers), 1)))
end

scratch = mktempdir()
for adjfib in fibers
    n = lpad(adjfib, 3, "0")
    ship = joinpath(OUT, "built", "APOGEE_skycont_svd_30_f$n.h5")
    isfile(ship) || (push!(failures, "EMPIRICAL: shipped skyCont missing for f$n"); continue)
    fresh = build_skyCont(adjfib; sample_dir=joinpath(OUT, "samples"),
        chipgap_msk_path=CHIPGAP, out_dir=scratch, nsub=30)
    V1 = h5read(ship, "Vmat"); V2 = h5read(fresh, "Vmat")
    l1 = h5read(ship, "λv"); l2 = h5read(fresh, "λv")
    okV = size(V1) == size(V2) && V1 == V2
    okl = l1 == l2
    @printf("(2) EMPIRICAL f%s: Vmat bitwise identical=%s  λv bitwise identical=%s  (maxabs dV=%.3e)\n",
        n, okV, okl, size(V1) == size(V2) ? maximum(abs.(V1 .- V2)) : NaN)
    (okV && okl) || push!(failures, "EMPIRICAL: skyCont rebuild differs for f$n")
end
rm(scratch, recursive=true, force=true)

println()
if isempty(failures)
    println("ASSERTION PASSED: skyCont products are threshold-policy invariant and reusable")
    println("  => a bright/faint threshold rebuild must redo skyLines ONLY; the existing")
    println("     APOGEE_skycont_svd_30_fNNN.h5 files are valid under options A/B/C/D alike.")
else
    println("ASSERTION FAILED:")
    foreach(f -> println("  - ", f), failures)
    exit(1)
end
