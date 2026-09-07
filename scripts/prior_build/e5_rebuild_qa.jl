## E5 rebuild QA: verify a policy-rebuilt skyLines prior set and compare it to the
## baseline built/ products.
#
# Answers three questions, in order:
#   A. INTEGRITY  — are all 600x3 products present, complete, and are the reused
#                   skyCont priors real, unbroken symlinks into the baseline?
#   B. MASK       — did the bright mask actually take effect, is it fiber-independent,
#                   and is it exactly `combined_mask ∩ that fiber's baseline submsk`?
#   C. SUBSPACE   — what did removing the bright pixels do to the priors?
#
# METHODOLOGY NOTE (read before adding a metric here).
# Prior SVD outputs are NOT bitwise reproducible: column signs flip run to run and
# rotations inside degenerate subspaces are arbitrary gauge. A bitwise or elementwise
# comparison reports enormous phantom differences. Principal angles are gauge-invariant
# but only where the spectrum is non-degenerate: at DEEP truncation the identity of the
# last retained mode swaps freely between builds (measured for these priors:
# lambda30/lambda1 = 1.7e-05, lambda30/lambda31 = 1.04), so a "89 deg at k=30" is
# truncation gauge and means nothing. See e5_qa_gauge_diag.jl and commit c216e10.
# Therefore this script LEADS with energy capture, then eigenvalues, and reports
# principal angles ONLY for low k where the modes carry real variance.
#
# Env: E5_OUT, E5_BUILT_NEW (basename of the rebuilt dir), E5_BRIGHT_COMBINED,
#      E5_MASK_VARIANT (dataset suffix in the combined file), E5_PLOTDIR, E5_QA_FIBERS
# Author - Andrew Saydjari (E5 pass 1)

using HDF5, Statistics, StatsBase, Printf, LinearAlgebra, Dates
using CairoMakie, ColorSchemes
black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")
BLAS.set_num_threads(2)

e5_out = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
new_tag = get(ENV, "E5_BUILT_NEW", "built_combined_telemaj_union")
plot_dir = get(ENV, "E5_PLOTDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e5_sky")
combined_path = get(ENV, "E5_BRIGHT_COMBINED", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined/e5_bright_combined.h5")
variant = get(ENV, "E5_MASK_VARIANT", "telemaj_union")
qa_fibers = [parse(Int, t) for t in split(get(ENV, "E5_QA_FIBERS", "10,76,295,351,460,519,600"), ",")]
mkpath(plot_dir)

base = joinpath(e5_out, "built")
new = joinpath(e5_out, new_tag)
report = String[]
say(s) = (push!(report, s); println(s); flush(stdout))

wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)

skycont_name(d, f) = joinpath(d, "APOGEE_skycont_svd_30_f" * lpad(f, 3, "0") * ".h5")
skyfaint_name(d, f) = joinpath(d, "APOGEE_skyline_faint_svd_120_f" * lpad(f, 3, "0") * ".h5")
skygspice_name(d, f) = joinpath(d, "APOGEE_skyline_faint_GSPICE_svd_120_f" * lpad(f, 3, "0") * ".h5")

say("E5 rebuild QA — $(Dates.now())")
say("baseline = $base")
say("rebuild  = $new")
say("mask     = $combined_path :: mask_$variant")

bright = h5read(combined_path, "mask_$variant")
say(@sprintf("combined bright mask: %d of %d grid pixels flagged", sum(bright), length(bright)))

## ---------------------------------------------------------------- A. INTEGRITY
# A "partial write" is a file that exists but cannot be fully read or has the wrong
# shape — the failure a plain `isfile` count would miss. So every product is OPENED
# and its dataset shapes checked, not merely counted.
struct FiberCheck
    f::Int
    cont_ok::Bool
    cont_islink::Bool
    cont_target_ok::Bool
    faint_ok::Bool
    gspice_ok::Bool
    err::String
end

function check_product(p, nmode)
    isfile(p) || return (false, "missing")
    try
        h5open(p, "r") do fh
            haskey(fh, "Vmat") && haskey(fh, "submsk") && haskey(fh, "λv") ||
                return (false, "datasets missing")
            size(fh["Vmat"]) == (8700, nmode) || return (false, "Vmat shape $(size(fh["Vmat"]))")
            size(fh["submsk"]) == (8700,) || return (false, "submsk shape")
            size(fh["λv"]) == (nmode,) || return (false, "λv shape")
            v = read(fh["λv"])
            all(isfinite, v) || return (false, "λv non-finite")
            return (true, "")
        end
    catch err
        return (false, "unreadable: " * first(sprint(showerror, err), 120))
    end
end

checks = Vector{FiberCheck}(undef, 600)
Threads.@threads for f in 1:600
    pc, pf, pg = skycont_name(new, f), skyfaint_name(new, f), skygspice_name(new, f)
    co, ce = check_product(pc, 30)
    fo, fe = check_product(pf, 120)
    go, ge = check_product(pg, 120)
    islk = islink(pc)
    # an unbroken symlink: the link resolves AND lands in the baseline built/ dir
    tgt_ok = false
    if islk
        try
            t = realpath(pc)
            tgt_ok = isfile(t) && dirname(t) == realpath(base)
        catch
            tgt_ok = false
        end
    end
    checks[f] = FiberCheck(f, co, islk, tgt_ok, fo, go, join(filter(!isempty, [ce, fe, ge]), "; "))
end

n_cont = count(c -> c.cont_ok, checks)
n_link = count(c -> c.cont_islink, checks)
n_tgt = count(c -> c.cont_target_ok, checks)
n_faint = count(c -> c.faint_ok, checks)
n_gspice = count(c -> c.gspice_ok, checks)
say("")
say("== A. INTEGRITY ==")
say(@sprintf("skyCont      present+readable %d/600 | symlinks %d/600 | resolving into baseline built/ %d/600",
    n_cont, n_link, n_tgt))
say(@sprintf("skyLines faint        present+readable %d/600", n_faint))
say(@sprintf("skyLines faint GSPICE present+readable %d/600", n_gspice))
bad = [c for c in checks if !(c.cont_ok && c.faint_ok && c.gspice_ok && c.cont_islink && c.cont_target_ok)]
if isempty(bad)
    say("integrity: PASS — 600/600 complete, no partial writes, all skyCont symlinks intact")
else
    say("integrity: FAIL — $(length(bad)) fiber(s) incomplete:")
    for c in first(bad, 20)
        say(@sprintf("  f%03d cont=%s link=%s tgt=%s faint=%s gspice=%s %s",
            c.f, c.cont_ok, c.cont_islink, c.cont_target_ok, c.faint_ok, c.gspice_ok, c.err))
    end
end
cl = joinpath(new, ".claims")
say(@sprintf("claims: %d (expect 600 if this dir was built in one pass)",
    isdir(cl) ? count(x -> isdir(joinpath(cl, x)), readdir(cl)) : -1))

## ------------------------------------------------------------------- B. MASK
# The whole point of the rebuild: the bright pixels must be GONE from the faint
# support, the removal must be identical on every fiber (fiber-independence), and it
# must be exactly `bright ∩ baseline submsk` — no more, no less.
say("")
say("== B. MASK — did the bright split actually take effect? ==")
n_old = zeros(Int, 600); n_new = zeros(Int, 600)
n_removed = zeros(Int, 600); n_added = zeros(Int, 600)
exact_ok = falses(600); have = falses(600)
Threads.@threads for f in 1:600
    po, pn = skyfaint_name(base, f), skyfaint_name(new, f)
    (isfile(po) && isfile(pn)) || continue
    try
        so = h5read(po, "submsk"); sn = h5read(pn, "submsk")
        n_old[f] = sum(so); n_new[f] = sum(sn)
        n_removed[f] = sum(so .& .!sn)
        n_added[f] = sum(sn .& .!so)
        # fiber-independence: the removed set must equal bright ∩ baseline submsk exactly
        exact_ok[f] = (so .& .!sn) == (so .& bright)
        have[f] = true
    catch
    end
end
nh = count(have)
if nh == 0
    say("no comparable fibers yet")
else
    idx = findall(have)
    fr = n_removed[idx] ./ max.(n_old[idx], 1)
    say(@sprintf("compared %d fibers", nh))
    say(@sprintf("baseline faint support: med %d px   rebuild: med %d px", median(n_old[idx]), median(n_new[idx])))
    say(@sprintf("pixels REMOVED (now bright): med %d  range %d-%d", median(n_removed[idx]), minimum(n_removed[idx]), maximum(n_removed[idx])))
    say(@sprintf("pixels ADDED to faint support: total %d  (MUST be 0)", sum(n_added[idx])))
    say(@sprintf("bright fraction of the faint support: med %.2f%%  range %.2f-%.2f%%",
        100median(fr), 100minimum(fr), 100maximum(fr)))
    say(@sprintf("fiber-independence (removed == bright ∩ baseline submsk, EXACT): %d/%d", count(exact_ok[idx]), nh))
    ok_b = sum(n_added[idx]) == 0 && count(exact_ok[idx]) == nh && minimum(n_removed[idx]) > 0
    say("mask: " * (ok_b ? "PASS — the split is live, fiber-independent, and exactly the delivered mask" :
                    "FAIL — see counts above"))
    # the finding-#35 no-op signature
    if maximum(n_removed[idx]) == 0
        say("*** NO-OP: zero pixels flagged bright on every fiber — this is finding #35 again ***")
    end
end

## --------------------------------------------------------------- C. SUBSPACE
# What the bright mask did to the priors. Leads with ENERGY, not angles.
say("")
say("== C. SUBSPACE — what the new mask changed in the priors ==")

princ_angles(V1, V2, k) = svdvals(Matrix(qr(V1[:, 1:k]).Q)' * Matrix(qr(V2[:, 1:k]).Q))

"""
    energy_capture(V1, V2, k)

Fraction of V1's leading-k energy lying in span(V2[:,1:k]). The meaningful
prior-vs-prior metric; see the methodology note at the top of this file.
"""
function energy_capture(V1, V2, k)
    Q2 = Matrix(qr(V2[:, 1:k]).Q)
    return sum(abs2, Q2' * V1[:, 1:k]) / sum(abs2, V1[:, 1:k])
end

# C1. eigenvalue spectra over all 600 fibers (cheap: λv only)
function lam_matrix(namer, dir, nsub)
    L = fill(NaN, nsub, 600)
    Threads.@threads for f in 1:600
        p = namer(dir, f)
        isfile(p) || continue
        try
            L[:, f] .= h5read(p, "λv")[1:nsub]
        catch
        end
    end
    return L
end
Lf_o = lam_matrix(skyfaint_name, base, 120); Lf_n = lam_matrix(skyfaint_name, new, 120)
Lg_o = lam_matrix(skygspice_name, base, 120); Lg_n = lam_matrix(skygspice_name, new, 120)

for (Lo, Ln, tag) in [(Lf_o, Lf_n, "faint"), (Lg_o, Lg_n, "faint GSPICE")]
    ok = findall(.!isnan.(Lo[1, :]) .& .!isnan.(Ln[1, :]))
    isempty(ok) && continue
    l1o = Lo[1, ok]; l1n = Ln[1, ok]
    tro = vec(sum(Lo[:, ok], dims=1)); trn = vec(sum(Ln[:, ok], dims=1))
    say(@sprintf("%s: lambda1 baseline med %.4g -> rebuild med %.4g  (ratio med %.3f)",
        tag, median(l1o), median(l1n), median(l1n ./ l1o)))
    say(@sprintf("%s: retained variance sum(lambda) baseline med %.4g -> rebuild med %.4g  (ratio med %.3f)",
        tag, median(tro), median(trn), median(trn ./ tro)))
end

# C2. per-mode bright-pixel energy of the BASELINE basis.
# This is the cleanest statement of what the mask removed: for baseline mode j,
# the fraction of its energy that lived on the pixels now declared bright. A
# leading mode with a large value was a sky-line mode contaminating the faint prior.
say("")
say("-- baseline modes: fraction of each mode's energy on the now-bright pixels --")
brightfrac_modes = fill(NaN, 120, length(qa_fibers))
for (i, f) in enumerate(qa_fibers)
    p = skygspice_name(base, f)
    isfile(p) || continue
    try
        V = h5read(p, "Vmat")
        so = h5read(p, "submsk")
        b = so .& bright                       # the pixels this fiber loses
        for j in 1:120
            e = sum(abs2, V[:, j])
            brightfrac_modes[j, i] = e > 0 ? sum(abs2, V[b, j]) / e : NaN
        end
    catch
    end
end
for (i, f) in enumerate(qa_fibers)
    all(isnan, brightfrac_modes[:, i]) && continue
    m = brightfrac_modes[:, i]
    say(@sprintf("  f%03d GSPICE: modes 1-5 %s | med over 120 modes %.3f | modes with >50%% bright energy: %d",
        f, join([@sprintf("%.3f", m[j]) for j in 1:5], " "), median(filter(isfinite, m)), count(x -> isfinite(x) && x > 0.5, m)))
end

# C3. energy capture + LOW-k principal angles on the COMMON support.
# Both bases are restricted to the rebuild's faint pixels (a subset of the baseline's),
# so the comparison is on identical rows. Angles reported only for k<=10.
say("")
say("-- subspace agreement on the common (faint) support --")
for f in qa_fibers
    for (namer, tag) in [(skyfaint_name, "faint"), (skygspice_name, "GSPICE")]
        po, pn = namer(base, f), namer(new, f)
        (isfile(po) && isfile(pn)) || continue
        try
            Vo = h5read(po, "Vmat"); Vn = h5read(pn, "Vmat")
            sn = h5read(pn, "submsk")
            Ao = Vo[sn, :]; An = Vn[sn, :]      # identical rows, common support
            capn = energy_capture(An, Ao, 10)   # new modes explained by baseline span
            capo = energy_capture(Ao, An, 10)   # baseline modes explained by rebuild span
            a = [acosd(clamp(minimum(princ_angles(Ao, An, k)), 0, 1)) for k in (1, 2, 3, 5, 10)]
            say(@sprintf("  f%03d %-6s k=10 capture new-in-old %.4f | old-in-new %.4f | angles k=1,2,3,5,10: %s deg",
                f, tag, capn, capo, join([@sprintf("%.2f", x) for x in a], ", ")))
        catch err
            say(@sprintf("  f%03d %-6s comparison failed: %s", f, tag, first(sprint(showerror, err), 100)))
        end
    end
end

## ----------------------------------------------------------------- FIGURES
try
    ok = findall(.!isnan.(Lg_o[1, :]) .& .!isnan.(Lg_n[1, :]))
    if !isempty(ok)
        fig = Figure(size=(1500, 950))
        ax1 = Axis(fig[1, 1], xlabel="adjfiberindx", ylabel="pixels",
            title="Faint-prior support: baseline vs $(new_tag)")
        lines!(ax1, ok, n_old[ok], color=:gray70, label="baseline submsk")
        lines!(ax1, ok, n_new[ok], color=:cyan, label="rebuild submsk")
        vlines!(ax1, [300.5], color=:orange, linewidth=1)
        axislegend(ax1, position=:rb)

        ax2 = Axis(fig[1, 2], xlabel="adjfiberindx", ylabel="bright px",
            title="Pixels removed from the faint support (now bright)")
        scatter!(ax2, ok, n_removed[ok], markersize=4, color=:magenta)
        vlines!(ax2, [300.5], color=:orange, linewidth=1)

        ax3 = Axis(fig[2, 1], xlabel="mode", ylabel="λ", yscale=log10,
            title="GSPICE eigenvalue spectra (median over fibers)")
        lines!(ax3, 1:120, [median(filter(isfinite, Lg_o[j, ok])) for j in 1:120], color=:gray70, label="baseline")
        lines!(ax3, 1:120, [median(filter(isfinite, Lg_n[j, ok])) for j in 1:120], color=:cyan, label="rebuild")
        axislegend(ax3, position=:rt)

        ax4 = Axis(fig[2, 2], xlabel="mode", ylabel="bright-pixel energy fraction",
            title="Baseline GSPICE modes: energy on the now-bright pixels")
        for (i, f) in enumerate(qa_fibers)
            all(isnan, brightfrac_modes[:, i]) && continue
            lines!(ax4, 1:120, brightfrac_modes[:, i], linewidth=1, label="f$(lpad(f,3,"0"))")
        end
        ylims!(ax4, -0.02, 1.02)
        axislegend(ax4, position=:rt, nbanks=2, labelsize=9)
        save(joinpath(plot_dir, "fig_rebuild_qa_$(variant).png"), fig, px_per_unit=2)
        say("")
        say("figure -> " * joinpath(plot_dir, "fig_rebuild_qa_$(variant).png"))
    end
catch err
    say("figure failed: " * first(sprint(showerror, err), 200))
end

open(joinpath(plot_dir, "rebuild_qa_$(variant).txt"), "w") do io
    for r in report
        println(io, r)
    end
end
println("\nreport -> " * joinpath(plot_dir, "rebuild_qa_$(variant).txt"))
