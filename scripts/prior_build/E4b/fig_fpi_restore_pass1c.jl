# pass-1c FPI restoration figure (APO fibers 76/226).
# LEFT PANELS: what the policy actually restores — the fiber's kept-entry list,
# with the chi2 character of the restored entries (pristine).
# RIGHT PANELS: the built prior's eigen-spectrum, pass-1 ceiling build vs the
# final pass-1c build. NOTE: unlike the LCO pass-1b case, lambda2 DECREASES
# here (x0.60 / x0.47) while lambda1 is unchanged — the restored bright-mode
# draws are pristine and homogeneous, so the subleading direction carries less
# variance. Fleet-wide lambda2 moves BOTH ways (median ratio 1.0001, 18 fibers
# <0.9, 20 >1.1, uncorrelated with gate removals), so this is a subleading-mode
# reshuffle, NOT a variance-restoration headline.
using CairoMakie, ColorSchemes, Printf, HDF5, Serialization, Statistics
black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")
P = "/mnt/ceph/users/sdssv/work/asaydjari/"
PLOTD = P * "2026_09_05/plots/c2_cut_analysis"
NEWB = P * "2026_09_05/prior_outputs/starCont_20260905_final/built"
OLDB = P * "2026_09_04/prior_outputs/starCont_pass1/built_apo"
new = deserialize(P * "2026_09_05/tfunlists_final/20260905_apo_tfunlist.jdat")
old = deserialize(P * "2026_09_03/tfunlists_refit20260902/20260902_apo_tfunlist.jdat")
mf = h5open(P * "2026_09_05/tfunlists_final/tfunlist_audit_apo.h5", "r") do f
    read(f["medflux"])
end
chif = h5open(P * "2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5", "r") do f
    Float64.(read(f["chi_sq_fiber"]))
end
lp(p) = h5open(p, "r") do f; read(f["λv"]); end

fig = Figure(size=(1250, 720))
for (row, fb) in ((1, 76), (2, 226))
    res = collect(setdiff(Set(new[fb]), Set(old[fb])))
    kept_old = collect(Set(old[fb]))
    ax = Axis(fig[row, 1], xscale=log10, xlabel="entry medflux", ylabel="count",
        title=@sprintf("fiberindx %d: kept %d → %d (+%d)", fb, length(old[fb]), length(new[fb]),
            length(new[fb]) - length(old[fb])))
    hist!(ax, mf[fb, kept_old], bins=10 .^ range(2.6, 5.2, 60), color=(:deepskyblue, 0.65),
        label="pass-1 kept (10k ceiling)")
    hist!(ax, mf[fb, res], bins=10 .^ range(2.6, 5.2, 60), color=(:springgreen, 0.8),
        label=@sprintf("RESTORED %d (χ² med %.1f)", length(res), median(chif[fb, res])))
    vlines!(ax, [10_000.0], color=:orangered, linestyle=:dash, linewidth=2)
    axislegend(ax, position=:lt, framevisible=false, labelsize=10)
    lo = lp(joinpath(OLDB, "APOGEE_starcont_svd_60_f" * lpad(fb, 3, "0") * ".h5"))
    ln = lp(joinpath(NEWB, "APOGEE_starcont_svd_60_f" * lpad(fb, 3, "0") * ".h5"))
    ax2 = Axis(fig[row, 2], yscale=log10, xlabel="component k", ylabel=L"\lambda_k",
        title=@sprintf("built prior (λ1 ratio %.4f, λ2 ratio %.3f)", ln[1] / lo[1], ln[2] / lo[2]))
    lines!(ax2, 1:60, lo, color=:orangered, linewidth=2, label="pass-1 build (10k ceiling)")
    lines!(ax2, 1:60, ln, color=:springgreen, linewidth=2, label="pass-1c build (final policy)")
    axislegend(ax2, position=:rt, framevisible=false, labelsize=10)
end
Label(fig[0, :], "pass-1c: the struck 10k ceiling returns the APO FPI fibers' pristine bright mode (χ² med ~5 vs fleet-kept 37); λ1 unchanged, λ2 reshuffled",
    fontsize=14)
save(joinpath(PLOTD, "fig_fpi_restore_pass1c.png"), fig, px_per_unit=2)
println("saved")
