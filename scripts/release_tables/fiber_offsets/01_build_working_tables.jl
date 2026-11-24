# 01_build_working_tables.jl
#
# CHAIN POSITION: fiber offsets step 1 of 5 (see ../README.md).
# INPUTS:  ARGS[1] WORKDIR, ARGS[2] allExposure (without v_rad_corr), ARGS[3] the 21-column allVisitRV, ARGS[4] seqmap.jls
# OUTPUTS: WORKDIR/build.h5, WORKDIR/build_counts.txt
# SOURCE:  work/2026_09_16/fiber_offsets/s1_build.jl (0.4.0). Logic unchanged; hard-coded paths became
#          ARGS (WORKDIR = the fiber-offset working directory shared by steps 1-5).
#
# s1_build.jl -- working tables for the per-fiber RV zero-point ("ubercal") derivation.
#
# Builds, from the READ-ONLY release products:
#   * an exposure working table: release default mask applied, plus the block index
#     b = adjfiberindx + 600*(era-1)  in 1..1200, which is the parameter index of the
#     (telescope, fiber, era) zero point.  adjfiberindx is 1..300 at APO, 301..600 at LCO,
#     so ONE parameter vector absorbs the fiber term, the telescope term and the era step,
#     which is what was asked for.
#   * the visit membership map (CSR) linking every delivered allVisitRV row to its member
#     exposure rows, rebuilt with the SAME code path as allvisitrv/s2_build.jl and then
#     required to reproduce the delivered v_rad_mean bitwise before anything is written.
#
# Membership uses the release default mask ONLY (no S/N requirement), which is what the
# delivered table was built with; the S/N requirement enters only in the offset solve.
#
# READ-ONLY outside 2026_09_16/fiber_offsets and 2026_09_16/plots/fiber_offsets.
using HDF5, Statistics, Printf, Serialization, Base.Threads

const DIR  = ARGS[1]
const PROD = ARGS[2]
const AVR  = ARGS[3]
const SEQ  = ARGS[4]

say(s) = (println(s); flush(stdout))

function main()
    t0 = time()
    say("reading allExposure ...")
    f = h5open(PROD, "r")
    tele   = read(f["tele"]);              mjd    = read(f["mjd"])
    expnum = read(f["expnum"]);            sid    = read(f["sdss_id"])
    mjdmid = read(f["mjd_mid_exposure"]);  v_rad  = read(f["v_rad"])
    snr    = read(f["snr"]);               era    = read(f["era"])
    cat    = read(f["category"]);          fib    = read(f["adjfiberindx"])
    dith   = read(f["dithered_pixels"]);   evr    = read(f["e_v_rad"])
    qa     = read(f["qa_bitmask"]);        ing    = read(f["ingestBit"])
    rvf    = read(f["RV_flag"]);           bf     = read(f["bary_flag"])
    close(f)
    N = length(mjd)
    say("  N = $N   ($(round(time()-t0,digits=1)) s)")

    mask = (qa .== 0) .& (ing .== 0) .& (rvf .== 0) .& (bf .== 0) .& isfinite.(v_rad)
    say(@sprintf("  release default mask: %d rows (%.3f %%)", count(mask), 100count(mask)/N))
    usable = mask .& isfinite.(snr) .& (snr .> 0)
    say("  ... and finite S/N > 0: $(count(usable))")

    # block index 1..1200 = adjfiberindx + 600*(era-1)
    blk = Vector{Int16}(undef, N)
    @threads for i in 1:N
        blk[i] = Int16(fib[i] + 600 * (Int(era[i]) - 1))
    end

    # ---------------------------------------------------------------- visit membership
    say("mapping onto almanac object sequences ...")
    seqmap = deserialize(SEQ)
    lut = Dict{Tuple{Int8,Int32}, Tuple{Vector{Int32},Vector{Int32},Vector{Int32}}}()
    gid = Int32(0)
    for k in sort(collect(keys(seqmap)))
        v = seqmap[k]
        st = Int32[a for (a,b) in v]; en = Int32[b for (a,b) in v]
        o = sortperm(st); st = st[o]; en = en[o]
        ids = Int32[(gid += Int32(1)) for _ in st]
        lut[k] = (st, en, ids)
    end
    seqid = zeros(Int32, N)
    @threads for i in 1:N
        e = get(lut, (tele[i], mjd[i]), nothing)
        e === nothing && continue
        st, en, ids = e
        j = searchsortedlast(st, expnum[i])
        (j >= 1 && expnum[i] <= en[j]) && (seqid[i] = ids[j])
    end
    say("  $gid sequences; rows inside one = $(count(>(0), seqid))")

    dlab = Vector{Float32}(undef, N)
    @threads for i in 1:N
        dlab[i] = (!isfinite(dith[i]) || dith[i] == 0f0) ? Float32(NaN) : dith[i]
    end

    valid = (seqid .> 0) .& (sid .> 0)
    idx = findall(valid)
    key = [(seqid[i], sid[i]) for i in idx]
    ord = sortperm(key); idx = idx[ord]; key = key[ord]
    NR = length(idx)
    bnd = Int[1]; for k in 2:NR; key[k] == key[k-1] || push!(bnd, k); end
    push!(bnd, NR + 1)
    NV = length(bnd) - 1
    say("  candidate visits = $NV")

    # reconstruct the balanced visits, keeping the member row indices
    R_sid = Int64[]; R_tele = Int8[]; R_mjd = Int32[]; R_mid = Float64[]; R_v = Float64[]
    memb_off = Int32[1]; memb = Int32[]
    buf = Int[]
    for b in 1:NV
        rg = bnd[b]:(bnd[b+1]-1)
        empty!(buf)
        for k in rg; r = idx[k]; mask[r] && push!(buf, r); end
        isempty(buf) && continue
        any(r -> isnan(dlab[r]), buf) && continue
        vals = unique(dlab[r] for r in buf)
        length(vals) == 2 || continue
        c1 = count(r -> dlab[r] == vals[1], buf)
        c1 == length(buf) - c1 || continue
        sort!(buf, by = r -> mjdmid[r])
        r0 = buf[1]
        push!(R_sid, sid[r0]); push!(R_tele, tele[r0]); push!(R_mjd, mjd[r0])
        push!(R_mid, mean(mjdmid[r] for r in buf))
        push!(R_v,   mean(v_rad[r]  for r in buf))
        append!(memb, Int32.(buf)); push!(memb_off, Int32(length(memb) + 1))
    end
    NRB = length(R_sid)
    say("  balanced visits reconstructed = $NRB   ($(round(time()-t0,digits=1)) s)")

    # ---------------------------------------------------------------- join onto delivered
    g = h5open(AVR, "r")
    D_sid = read(g["sdss_id"]); D_mjd = read(g["mjd"]); D_mid = read(g["mjd_mid"])
    D_v = read(g["v_rad_mean"]); D_tele_raw = read(g["tele"]); D_n = read(g["n_exp"])
    close(g)
    ND = length(D_sid)
    D_tele = Vector{Int8}(undef, ND)
    for i in 1:ND
        s = D_tele_raw[i]
        D_tele[i] = (s isa String ? s : String(copy(reinterpret(UInt8, [s])))) == "apo" ? Int8(1) : Int8(2)
    end
    kd = Dict{Tuple{Int8,Int32,Int64,Int64}, Int}()
    for i in 1:ND; kd[(D_tele[i], D_mjd[i], D_sid[i], round(Int64, D_mid[i]*1e8))] = i; end
    # delivered row -> reconstructed visit
    d2r = zeros(Int32, ND)
    for i in 1:NRB
        j = get(kd, (R_tele[i], R_mjd[i], R_sid[i], round(Int64, R_mid[i]*1e8)), 0)
        j > 0 && (d2r[j] = Int32(i))
    end
    nhit = count(>(0), d2r)
    say(@sprintf("  delivered rows matched: %d / %d (%.5f %%)", nhit, ND, 100nhit/ND))
    @assert nhit == ND "visit membership rebuild is not complete -- $(ND-nhit) delivered rows unmatched"
    dv = maximum(abs(R_v[d2r[i]] - D_v[i]) for i in 1:ND)
    dn = maximum(abs(Int(memb_off[d2r[i]+1]-memb_off[d2r[i]]) - Int(D_n[i])) for i in 1:ND)
    say(@sprintf("  max |v_reconstructed - v_delivered| = %.3e km/s ;  max |n_exp diff| = %d", dv, dn))
    @assert dv == 0.0 && dn == 0 "rebuild does not reproduce the delivered table"

    # re-emit membership in DELIVERED row order
    off = Int32[1]; mm = Int32[]
    for i in 1:ND
        r = d2r[i]
        for k in memb_off[r]:(memb_off[r+1]-1); push!(mm, memb[k]); end
        push!(off, Int32(length(mm)+1))
    end
    say("  membership entries = $(length(mm))")

    # ---------------------------------------------------------------- write
    E = findall(usable)
    out = joinpath(DIR, "build.h5")
    isfile(out) && rm(out)
    h5open(out, "w") do h
        h["exp/row"]      = Int32.(E)          # 1-based row in allExposure_pass1_dr21.h5
        h["exp/sdss_id"]  = sid[E]
        h["exp/tele"]     = tele[E]
        h["exp/era"]      = Vector{Int8}(era[E])
        h["exp/block"]    = blk[E]
        h["exp/fiber"]    = fib[E]
        h["exp/mjd"]      = mjd[E]
        h["exp/mjd_mid"]  = mjdmid[E]
        h["exp/v_rad"]    = v_rad[E]
        h["exp/e_v_rad"]  = evr[E]
        h["exp/snr"]      = snr[E]
        h["exp/category"] = cat[E]
        # full-length block index so the corrected column can be written for every masked row
        h["all/block"]    = blk
        h["all/mask"]     = Vector{UInt8}(mask)
        h["all/v_rad"]    = v_rad
        # visit membership, delivered row order
        h["visit/memb_off"] = off
        h["visit/memb_row"] = mm
    end
    chmod(out, 0o644)
    say("wrote $out  ($(round(time()-t0,digits=1)) s)")

    open(joinpath(DIR, "build_counts.txt"), "w") do io
        println(io, "exposure_rows=$N")
        println(io, "exposure_rows_default_mask=$(count(mask))")
        println(io, "exposure_rows_usable_snr=$(count(usable))")
        println(io, "delivered_visit_rows=$ND")
        println(io, "visit_membership_entries=$(length(mm))")
        println(io, "max_abs_v_recon_minus_delivered_kms=$dv")
        println(io, "max_abs_nexp_diff=$dn")
    end
    chmod(joinpath(DIR, "build_counts.txt"), 0o644)
end

main()
