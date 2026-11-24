# 05_apply_allvisitrv.jl
#
# CHAIN POSITION: fiber offsets step 5 of 5 (see ../README.md).
# INPUTS:  ARGS[1] WORKDIR (fiber_offsets.csv, build.h5), ARGS[2] the 21-column allVisitRV, ARGS[3] output path, ARGS[4] optional text naming the input in the `provenance` root attribute (default ARGS[2]; the 0.4.0 file names its scifest path)
# OUTPUTS: ARGS[3] allVisitRV with the three v_rad_corr_* columns appended (24 datasets), WORKDIR/corr.h5
# SOURCE:  work/2026_09_16/fiber_offsets/s4_apply.jl (0.4.0). Logic unchanged; hard-coded paths became
#          ARGS (WORKDIR = the fiber-offset working directory shared by steps 1-5).
#          The exposure-level half of s4 (staged_allExposure_v_rad_corr.h5) was superseded in 0.4.0
#          by s10_final_build.jl and is now the S_FIBEROFFSET stage of the Julia allExposure builder;
#          its file write was removed here (the per-row values are still computed: the visit
#          reduction needs them).
#
# s4_apply.jl -- apply the offset table and STAGE the corrected products.
#
# Writes two files, both in this directory, both staged for review (0.4.0 wording).  Nothing is released
# and neither input is touched.
#
#  1. staged_allExposure_v_rad_corr.h5   -- the exposure-level delta.  One velocity column
#     (v_rad_corr) plus the per-row offset actually subtracted and a flag, ROW-ALIGNED with
#     allExposure_pass1_dr21.h5: row i here is row i there, all 26,715,540 of them.
#
#  2. staged_allVisitRV_pass1_dr21_corr.h5 -- the full 21-column visit table carried across
#     verbatim, with three columns appended:
#       v_rad_corr_mean          plain unweighted mean of the members' v_rad_corr
#       v_rad_corr_scatter       sample std (n-1) of the members' v_rad_corr
#       v_rad_corr_scatter_jack  jackknifed (Eq. 8 on Eq. 10), NaN at n_exp == 2
#     Same membership, same plain-mean convention, same NaN-at-n_exp==2 rule as shipped.
#
# NaN POLICY.  A non-calibratable block has offset = NaN, so its exposures get
# v_rad_corr = NaN and any visit containing one gets NaN for all three corrected columns.
# That is the DR17 paper's own choice for its three uncalibratable LCO fibers, and it is why
# the count of NaN rows is reported here rather than hidden.
using HDF5, CSV, DataFrames, Statistics, Printf, Base.Threads
using Jackknife

const DIR  = ARGS[1]
const AVR  = ARGS[2]
const AVR_LABEL = length(ARGS) >= 4 ? ARGS[4] : AVR
say(s) = (println(s); flush(stdout))

function main()
    t0 = time()
    off = fill(NaN, 1200); cal = falses(1200)
    d = CSV.read(joinpath(DIR, "fiber_offsets.csv"), DataFrame; comment = "#")
    for r in eachrow(d)
        b = r.fiber + (r.tele == "lco" ? 300 : 0) + 600*(r.era - 1)
        off[b] = r.offset_kms; cal[b] = r.calibratable
    end
    say("offset table: $(count(cal)) calibratable blocks")

    h = h5open(joinpath(DIR, "build.h5"), "r")
    blk = Int.(read(h["all/block"])); msk = read(h["all/mask"]) .== 0x01
    v   = read(h["all/v_rad"])
    memb_off = read(h["visit/memb_off"]); memb_row = read(h["visit/memb_row"])
    close(h)
    N = length(v)

    # ------------------------------------------------------------ exposure level
    vc = fill(NaN, N); ov = fill(NaN, N); fl = zeros(UInt8, N)
    @threads for i in 1:N
        if !msk[i]
            fl[i] = 0x02                      # row outside the release default mask
        elseif !cal[blk[i]]
            fl[i] = 0x01                      # masked-in, but its block is not calibratable
        else
            ov[i] = off[blk[i]]; vc[i] = v[i] - ov[i]
        end
    end
    say(@sprintf("exposure rows: %d total ; corrected %d ; NaN because not calibratable %d ; NaN because outside the mask %d",
        N, count(isfinite, vc), count(==(0x01), fl), count(==(0x02), fl)))
    say(@sprintf("  offsets applied: median %+.1f m/s, range %+.1f .. %+.1f m/s",
        1000median(filter(isfinite, ov)), 1000minimum(filter(isfinite, ov)), 1000maximum(filter(isfinite, ov))))

    # ------------------------------------------------------------ visit level
    ND = length(memb_off) - 1
    vm = fill(NaN, ND); vs = fill(NaN, ND); vj = fill(NaN, ND)
    nbad = Threads.Atomic{Int}(0)
    @threads for i in 1:ND
        lo = memb_off[i]; hi = memb_off[i+1] - 1
        n = hi - lo + 1
        x = Vector{Float64}(undef, n)
        ok = true
        @inbounds for k in 1:n
            xv = vc[memb_row[lo+k-1]]
            isfinite(xv) || (ok = false; break)
            x[k] = xv
        end
        ok || (Threads.atomic_add!(nbad, 1); continue)
        vm[i] = mean(x)
        vs[i] = std(x)                                    # (n-1) normalised, Eq. 10
        n > 2 && (vj[i] = Jackknife.estimate(std, x))     # Eq. 8 on Eq. 10
    end
    say(@sprintf("visit rows: %d ; corrected %d ; NaN (a member on a non-calibratable block) %d",
        ND, count(isfinite, vm), nbad[]))
    say(@sprintf("  v_rad_corr_scatter      median %.4f km/s (delivered v_rad_scatter median 0.0696)", median(filter(isfinite, vs))))
    say(@sprintf("  v_rad_corr_scatter_jack median %.4f km/s on %d rows (delivered 0.137 on 1,967,077)",
        median(filter(isfinite, vj)), count(isfinite, vj)))

    say("copying the visit table and appending ...")
    out2 = ARGS[3]
    isfile(out2) && rm(out2)
    src = h5open(AVR, "r"); dst = h5open(out2, "w")
    # copy_object preserves dtype (incl. the |S3 `tele`), chunking, filters and attributes
    for k in sort(collect(keys(src))); copy_object(src, k, dst, k); end
    for (k, val, desc) in (
        ("v_rad_corr_mean", vm, "unweighted mean of the member v_rad_corr, km/s, barycentric. Same membership and same plain-mean convention as v_rad_mean. NaN if any member sits on a non-calibratable (tele, fiber, era) block."),
        ("v_rad_corr_scatter", vs, "sample standard deviation (n-1) of the member v_rad_corr, km/s. Defined wherever v_rad_corr_mean is."),
        ("v_rad_corr_scatter_jack", vj, "jackknifed (debiased) v_rad_corr_scatter, Saydjari et al. 2025 Eq. 8 applied to Eq. 10, km/s. NaN when n_exp == 2, exactly as v_rad_scatter_jack."))
        dst[k, chunk = (65536,), shuffle = (), deflate = 9] = val
        attrs(dst[k])["description"] = desc
    end
    attrs(dst)["provenance"] = "STAGED, NOT RELEASED. Columns 1-21 are byte-identical copies of " * AVR_LABEL *
        "; the three v_rad_corr_* columns were added by 2026_09_16/fiber_offsets/s4_apply.jl using fiber_offsets.csv."
    close(src); close(dst)
    chmod(out2, 0o644)

    # verify the carried columns are unchanged
    a = h5open(AVR, "r"); b = h5open(out2, "r")
    bad = String[]
    for k in sort(collect(keys(a)))
        isequal(read(a[k]), read(b[k])) || push!(bad, k)
    end
    close(a); close(b)
    say(isempty(bad) ? "  verified: all 21 carried columns identical to the released file" :
                       "  MISMATCH in: " * join(bad, ", "))
    @assert isempty(bad)

    h5open(joinpath(DIR, "corr.h5"), "w") do g
        g["visit/v_corr_mean"] = vm; g["visit/v_corr_scatter"] = vs; g["visit/v_corr_jack"] = vj
    end
    chmod(joinpath(DIR, "corr.h5"), 0o644)
    say("wrote $out2  ($(round(time()-t0,digits=1)) s)")
end

main()
