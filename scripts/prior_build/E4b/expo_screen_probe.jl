using HDF5, Statistics, StatsBase, Printf
for tele in ("apo", "lco")
    aud = tele == "apo" ?
        "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5" :
        "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/tfunlist_audit_lco.h5"
    mf = h5open(aud, "r") do f; read(f["medflux"]); end
    N = size(mf, 2)
    em = [median(@view mf[:, i]) for i in 1:N]
    m = median(em); s = mad(em, normalize=true)
    for k in (5, 8, 10, 15, 20)
        thr = m + k * s
        @printf("[%s] med %.0f MAD %.0f  k=%2d thr %.0f  exposures caught: %d\n",
            tele, m, s, k, thr, count(>(thr), em))
    end
    thr10 = m + 10 * s
    caught = findall(>(thr10), em)
    top = sort(em, rev=true)[1:min(12, N)]
    println("[$tele] top-12 exposure medians: ", join(round.(Int, top), " "))
    if tele == "apo"
        nsub = sum(count(x -> 400 < x <= 10_000, @view mf[:, i]) for i in caught; init=0)
        @printf("[apo] k=10 catches %d exposures; their (400,10k] entries (NEW removals beyond ceiling): %d\n",
            length(caught), nsub)
        # per caught exposure: how many entries <=10k
        subs = sort([count(x -> 400 < x <= 10_000, @view mf[:, i]) for i in caught], rev=true)
        println("[apo] per-caught-exposure sub-10k entry counts (desc): ", join(subs[1:min(15, end)], " "))
    else
        println("[lco] k=10 caught rows: ", caught, " medians ", join(round.(Int, em[caught]), " "))
    end
end
