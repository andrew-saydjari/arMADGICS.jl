#= 01_almanac_scan.jl
   CHAIN POSITION: offtarget step 1 of 3. Needed only for the VISIT-level flag (step 2); the
     exposure-level `offtarget` column is computed inside the Julia allExposure builder (stage
     S_OFFTARGET, same predicate), which replaced 0.4.0's s2_derive_flag.jl.
   INPUTS:  ARGS[1] the almanac file.   OUTPUTS: ARGS[2] almanac_fps_flags.h5
   SOURCE:  work/2026_09_16/offtarget_sidecar/s1_almanac_scan.jl (0.4.0). Logic unchanged; paths became ARGS.
   s1: scan the pass-1 almanac for the FPS off-target predicate.
   For every FPS-era (obs, mjd, config_id) fiber table, record per almanac fiber index:
     assigned, on_target, valid, disabled, and the category string code.
   Output: almanac_fps_flags.h5 with parallel arrays + a (nfib x nconfig) matrix per field.
   READ-ONLY on the almanac. =#
using HDF5, Printf

const A = ARGS[1]
const OUT = ARGS[2]
say(s) = (println(s); flush(stdout))

function main()
    t0 = time()
    h = h5open(A, "r")
    teles = Int8[]; mjds = Int32[]; cfgs = Int32[]
    ASG = Vector{Vector{Bool}}(); ONT = Vector{Vector{Bool}}()
    VAL = Vector{Vector{Bool}}(); DIS = Vector{Vector{Bool}}()
    CAT = Vector{Vector{UInt8}}()
    catcodes = Dict{String,UInt8}("" => UInt8(0))
    catlist = String[""]
    nbadlen = 0
    for (ti, obs) in enumerate(("apo", "lco"))
        nights = sort(parse.(Int, collect(keys(h["raw/$obs"]))))
        nfps = 0
        for m in nights
            g = h["raw/$obs/$m"]
            haskey(g, "fibers") || continue
            gf = g["fibers"]
            ks = collect(keys(gf))
            isempty(ks) && continue
            haskey(gf[ks[1]], "on_target") || continue   # plate era: no such column
            nfps += 1
            for cfg in ks
                gg = gf[cfg]
                haskey(gg, "on_target") || continue
                asg = read(gg["assigned"]); ont = read(gg["on_target"])
                val = read(gg["valid"]);    dis = read(gg["disabled"])
                cat = read(gg["category"])
                n = length(asg)
                if !(length(ont) == n && length(val) == n && length(dis) == n && length(cat) == n)
                    nbadlen += 1; continue
                end
                cc = Vector{UInt8}(undef, n)
                for i in 1:n
                    s = String(strip(String(cat[i])))
                    c = get(catcodes, s, UInt8(255))
                    if c == 255
                        push!(catlist, s); c = UInt8(length(catlist) - 1); catcodes[s] = c
                    end
                    cc[i] = c
                end
                push!(teles, Int8(ti)); push!(mjds, Int32(m)); push!(cfgs, Int32(parse(Int, cfg)))
                push!(ASG, Vector{Bool}(asg)); push!(ONT, Vector{Bool}(ont))
                push!(VAL, Vector{Bool}(val)); push!(DIS, Vector{Bool}(dis)); push!(CAT, cc)
            end
        end
        say("$obs: $(length(nights)) nights, $nfps FPS-era nights, cumulative configs = $(length(cfgs))  [$(round(time()-t0,digits=1)) s]")
    end
    close(h)
    say("bad-length fiber tables skipped: $nbadlen")
    lens = unique(length.(ASG))
    say("fiber-table lengths seen: $lens")
    NF = maximum(length.(ASG)); NC = length(cfgs)
    say("configs = $NC, matrix $(NF) x $(NC)")
    function pack(V, T)
        M = zeros(T, NF, NC)
        for j in 1:NC
            v = V[j]
            @inbounds for i in 1:length(v); M[i, j] = T(v[i]); end
        end
        M
    end
    Masg = pack(ASG, UInt8); Mont = pack(ONT, UInt8)
    Mval = pack(VAL, UInt8); Mdis = pack(DIS, UInt8); Mcat = pack(CAT, UInt8)
    nlen = Int16[length(v) for v in ASG]
    h5open(OUT, "w") do f
        f["tele"] = teles; f["mjd"] = mjds; f["config_id"] = cfgs; f["nfiber"] = nlen
        f["assigned"] = Masg; f["on_target"] = Mont; f["valid"] = Mval
        f["disabled"] = Mdis; f["category_code"] = Mcat
        attrs(f)["categories"] = join(catlist, "|")
        attrs(f)["note"] = "FPS-era almanac fiber flags, indexed [almanac_fiber_index, config]"
    end
    say("categories: $(catlist)")
    # quick prevalence, corpus-wide over the almanac itself
    off = (Masg .== 1) .& (Mont .== 0)
    for (ci, cn) in enumerate(catlist)
        sel = Mcat .== UInt8(ci - 1)
        na = count(sel .& (Masg .== 1)); no = count(sel .& off)
        na == 0 && continue
        @printf("  cat %-16s assigned=%9d  off-target=%8d  (%.3f %%)\n", cn, na, no, 100no/na)
    end
    @printf("ALL assigned=%d off=%d (%.3f %%)\n", count(Masg .== 1), count(off), 100count(off)/count(Masg .== 1))
    say("wrote $OUT  [$(round(time()-t0,digits=1)) s]")
end
main()
