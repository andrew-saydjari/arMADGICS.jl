using CairoMakie, Printf, HDF5, Statistics, StatsBase, Serialization, Random
black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")
PLOTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/plots/c2_cut_analysis"
mf = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5","r") do f; read(f["medflux"]); end
chif = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5","r") do f; Float64.(read(f["chi_sq_fiber"])); end
N = size(mf,2); base = mf .> 400.0
fpisel = falses(300,N); fpisel[[76,226],:] .= true
rng = MersenneTwister(11)

begin
    fig = Figure(size=(1150, 500))
    ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="medflux (per fiber,exposure)",
        ylabel="chi2_fiber", title="APO: the >10k population separates by CLASS, not brightness")
    nfidx = findall(vec(base .& .!fpisel))
    keep = randsubseq(rng, nfidx, 30_000/length(nfidx))
    scatter!(ax, vec(mf)[keep], max.(vec(chif)[keep], 0.5), color=(:deepskyblue, 0.25), markersize=2, label="non-FPI fibers (subsampled)")
    # always include ALL non-FPI >10k points (the bad class)
    bad = findall(vec(base .& .!fpisel .& (mf .> 10_000)))
    scatter!(ax, vec(mf)[bad], max.(vec(chif)[bad], 0.5), color=(:deepskyblue, 0.5), markersize=2.5)
    fi = findall(vec(base .& fpisel))
    fkeep = randsubseq(rng, fi, min(1.0, 15_000/length(fi)))
    scatter!(ax, vec(mf)[fkeep], max.(vec(chif)[fkeep], 0.5), color=(:orangered, 0.5), markersize=2.5, label="FPI fibers 76/226")
    vlines!(ax, [10_000], color=:magenta, linestyle=:dash, label="old 10k ceiling")
    hlines!(ax, [58328.468], color=:yellow, linestyle=:dashdot, label="C3 gate (p99.9)")
    text!(ax, 11_000, 1.5, text="FPI bright mode: chi2 ~ 6 (pristine)", color=:orangered, fontsize=13)
    text!(ax, 11_000, 3000.0, text="39 bad exposures: chi2 ~ 5e4", color=:deepskyblue, fontsize=13)
    axislegend(ax, position=:lt, framevisible=false)
    save(joinpath(PLOTD, "fig_apo_class_money.png"), fig, px_per_unit=2)
    println("money saved")
end
begin
    fig = Figure(size=(1150, 420))
    for (i, fb) in enumerate((76, 226))
        ax = Axis(fig[1,i], xscale=log10, xlabel="medflux", ylabel="entries",
            title="APO fiberindx $fb (FPI): bimodal — bright mode is NORMAL for it")
        v = mf[fb, base[fb,:]]
        hist!(ax, v, bins=10 .^ range(log10(401), log10(1e5), 80), color=(:orangered, 0.8))
        vlines!(ax, [10_000], color=:magenta, linestyle=:dash)
        vlines!(ax, [2*median(v)], color=:cyan, linestyle=:dot)
    end
    Label(fig[0,:], "APO FPI-fiber bimodality (fbright 17%, q90/q50 ~5; next-highest fiber 2.9%) — magenta: old 10k ceiling; cyan: 2x fiber median", fontsize=14)
    save(joinpath(PLOTD, "fig_apo_bimodal.png"), fig, px_per_unit=2)
    println("bimodal saved")
end
begin
    pr = deserialize("apo_E2_probe.jdat")
    fig = Figure(size=(1000, 420))
    ax = Axis(fig[1,1], yscale=log10, xlabel="exposure rank (by non-FPI median chi2)",
        ylabel="non-FPI median chi2", title="APO exposure fit-quality tail: catastrophic family vs continuum")
    s = sort(pr.ec, rev=true)[1:200]
    scatter!(ax, 1:200, s, color=:deepskyblue, markersize=5)
    bset = Set(pr.bright); ord = sortperm(pr.ec, rev=true)[1:200]
    bmask = [o in bset for o in ord]
    scatter!(ax, (1:200)[bmask], s[bmask], color=:orangered, markersize=6, label="bright-median (>10k) exposures")
    hlines!(ax, [10_000], color=:springgreen, linestyle=:dash, label="proposed E2 threshold (~650x fleet med 15.5)")
    axislegend(ax, position=:rt, framevisible=false)
    save(joinpath(PLOTD, "fig_apo_expchi2_tail.png"), fig, px_per_unit=2)
    println("tail saved")
end
