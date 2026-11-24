#= 02_sidecar_and_visits.jl
   CHAIN POSITION: offtarget step 2 of 3 (after the allExposure build and 01_almanac_scan.jl,
     before 03_stage_allvisitrv.jl).
   INPUTS:  ARGS[1] final allExposure (with `offtarget`), ARGS[2] allVisitRV (24 columns),
            ARGS[3] almanac_fps_flags.h5 from step 1
   OUTPUTS: ARGS[4] offtarget_exposures.csv (released sidecar), ARGS[5] offtarget_visit_flag.h5
   SOURCE:  work/2026_09_16/offtarget_sidecar/s3_sidecar_and_visits.jl (0.4.0). Logic unchanged,
            except the exposure flag is read from the allExposure `offtarget` column instead of
            s2's offtarget_exposure_flag.h5 (same values: the builder computes the s2 predicate).
   s3: (a) write the offtarget_exposures.csv sidecar, (b) derive the visit-level flag and
   verify it is constant within every visit.
   Visit-level path is INDEPENDENT of the exposure-level path: it keys the almanac by
   (tele, mjd, config_id) and almanac index 301 - fiber, where the visit file's `fiber` is
   adjfiberindx - 300*(tele-1). Agreement between the two paths is the cross-check. =#
using HDF5, Printf, Statistics

const AEXP = ARGS[1]
const AVIS = ARGS[2]
const ALM  = ARGS[3]
const CSVOUT = ARGS[4]
const VISOUT = ARGS[5]
say(s) = (println(s); flush(stdout))

function main()
    t0 = time()
    flag = h5open(AEXP, "r") do g; read(g["offtarget"]); end
    f = h5open(AEXP, "r")
    tele = read(f["tele"]); mjd = read(f["mjd"]); expnum = read(f["expnum"])
    adj = read(f["adjfiberindx"]); sid = read(f["sdss_id"]); cfg = read(f["config_id"])
    close(f)
    N = length(mjd)
    telestr = ["apo", "lco"]

    # ---------- (a) sidecar CSV ----------
    idx = findall(flag .== 1)
    say("flagged exposure rows: $(length(idx))")
    open(CSVOUT, "w") do io
        println(io, "# APOGEE DR21 pass-1 (0.4.0): exposure x fiber rows whose FPS robot did not reach its assigned target.")
        println(io, "# predicate: almanac assigned == 1 and on_target == 0; the spectrum is of whatever sky the fiber reached, not of sdss_id.")
        println(io, "# FPS era only (MJD >= 59558); the plate era has no assigned/on_target columns and is not flagged here.")
        println(io, "tele,mjd,expnum,adjfiberindx,sdss_id")
        for i in idx
            println(io, telestr[tele[i]], ",", mjd[i], ",", expnum[i], ",", adj[i], ",", sid[i])
        end
    end
    say("wrote offtarget_exposures.csv")

    # ---------- (b) visit level ----------
    a = h5open(ALM, "r")
    at = read(a["tele"]); am = read(a["mjd"]); ac = read(a["config_id"])
    Masg = read(a["assigned"]); Mont = read(a["on_target"])
    close(a)
    lut = Dict{Tuple{Int8,Int32,Int32},Int32}()
    for j in 1:length(ac); lut[(at[j], am[j], ac[j])] = Int32(j); end

    v = h5open(AVIS, "r")
    vtele_raw = read(v["tele"]); vmjd = read(v["mjd"]); vcfg = read(v["config_id"])
    vfib = read(v["fiber"]); vsid = read(v["sdss_id"]); vera = read(v["era"])
    vnexp = read(v["n_exp"]); vsnr = read(v["snr_mean"])
    close(v)
    NV = length(vmjd)
    vtele = Int8[String(strip(String(x))) == "apo" ? Int8(1) : Int8(2) for x in vtele_raw]
    say("allVisitRV rows = $NV")

    vflag = zeros(Int8, NV); vmatched = falses(NV)
    @inbounds for i in 1:NV
        j = get(lut, (vtele[i], vmjd[i], vcfg[i]), Int32(0))
        j == 0 && continue
        k = 301 - Int(vfib[i])
        (1 <= k <= 300) || continue
        vmatched[i] = true
        vflag[i] = (Masg[k, j] == 1 && Mont[k, j] == 0) ? Int8(1) : Int8(0)
    end
    nv1 = count(vflag .== 1)
    @printf("visits with an FPS config match: %d (era==2 visits: %d)\n", count(vmatched), count(vera .== 2))
    @printf("FLAGGED visits: %d  (%.4f %% of all visits, %.4f %% of FPS visits)\n",
            nv1, 100nv1/NV, 100nv1/max(1, count(vera .== 2)))
    say("era==2 visits without a config match: $(count((vera .== 2) .& .!vmatched))")
    say("era!=2 visits with a config match: $(count((vera .!= 2) .& vmatched))")

    # ---------- constancy verification against the exposure-level flag ----------
    # group exposure rows by (tele, mjd, config_id, fiber, sdss_id); the visit definition
    # (contiguous object sequence sharing config_id, one star) is a refinement of this key.
    say("\n== verifying offtarget is constant within every visit ==")
    gmap = Dict{Tuple{Int8,Int32,Int32,Int16,Int64},Int8}()
    nsplit = 0
    @inbounds for i in 1:N
        fib = Int16(Int(adj[i]) - 300 * (Int(tele[i]) - 1))
        k = (tele[i], mjd[i], cfg[i], fib, sid[i])
        p = get(gmap, k, Int8(-1))
        if p == -1
            gmap[k] = flag[i]
        elseif p != flag[i]
            nsplit += 1
        end
    end
    say("exposure groups (tele,mjd,config_id,fiber,sdss_id): $(length(gmap))")
    say(nsplit == 0 ? "  offtarget is CONSTANT within every such group (0 split groups)" :
                      "  *** $nsplit exposure rows disagree with their group -- SPLIT VISITS ***")

    nmiss = 0; ndis = 0
    @inbounds for i in 1:NV
        k = (vtele[i], vmjd[i], vcfg[i], vfib[i], vsid[i])
        p = get(gmap, k, Int8(-1))
        if p == -1; nmiss += 1
        elseif p != vflag[i]; ndis += 1
        end
    end
    @printf("visits whose key has no exposure group: %d\n", nmiss)
    @printf("visits where the visit-level flag disagrees with the exposure-level flag: %d\n", ndis)
    (nsplit == 0 && ndis == 0 && nmiss == 0) || error("visit constancy verification FAILED")
    say("  VISIT CONSTANCY: PASS")

    # a few descriptive numbers
    for (lbl, s) in (("on-target", (vera .== 2) .& (vflag .== 0)), ("off-target", vflag .== 1))
        x = filter(isfinite, vsnr[s])
        @printf("  FPS visits %-11s n=%8d  median snr_mean=%8.2f  median n_exp=%.1f\n",
                lbl, count(s), isempty(x) ? NaN : median(x), median(vnexp[s]))
    end

    # task-94 star at visit level
    say("\n== task-94 star at visit level ==")
    s = findall(vsid .== 85604486)
    for i in sort(s, by = j -> vmjd[j])
        @printf("  MJD %d cfg %6d fiber %3d n_exp %d  offtarget=%d\n",
                vmjd[i], vcfg[i], vfib[i], vnexp[i], vflag[i])
    end

    h5open(VISOUT, "w") do g
        g["offtarget"] = vflag
    end
    say("\nwrote offtarget_visit_flag.h5  [$(round(time()-t0,digits=1)) s]")
end
main()
