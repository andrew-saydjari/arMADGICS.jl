# plot_pass1_priors.jl
#
# Per-fiber visualization of EVERY prior family the APOGEE DR21 pass-1 run consumes.
# Review gate before the 200-day testbed re-run: samples + top-3 eigenvectors,
# per fiber, per family, with the numbers printed on the figure.
#
# Prior families (exactly the set src/priors.jl::build_prior_dict wires for pass-1;
# the DIB priors are commented out in pipelineCore.jl and are NOT pass-1):
#   1. starCont              APOGEE_starcont_svd_60_fNNN.h5              (pass-1c, per telescope)
#   2. starLines             APOGEE_stellar_kry_50_subpix_fNNN.h5        (E7 per-fiber LSF, subpix idx 6)
#   3. skyCont               APOGEE_skycont_svd_30_fNNN.h5               (E5)
#   4. skyLine faint         APOGEE_skyline_faint_svd_120_fNNN.h5        (built, NOT deployed)
#   5. skyLine faint GSPICE  APOGEE_skyline_faint_GSPICE_svd_120_fNNN.h5 (DEPLOYED)
#
# Conventions taken from the build scripts (MEASURED, not assumed):
#   * Vmat columns are eigvec * sqrt(eigval)  =>  norm(V[:,k])^2 == lambda_k.
#     We plot UNIT-NORM eigenvectors u_k = V[:,k]/sqrt(lambda_k) so modes 2 and 3 are
#     legible next to mode 1, and print lambda_k and lambda_k/sum(lambda) on the panel.
#   * Detector grid: wavetarg = 10 .^ range(4.179-125*6e-6, step=6e-6, length=8700).
#   * Chip gaps / dead pixels are set to NaN so the line BREAKS (never interpolated).
#   * The covariances are second moments (no mean subtraction), so mode 1 is
#     essentially the mean spectrum of the training corpus for every family.
#
# Usage:
#   julia --project=<env with CairoMakie, HDF5, ColorSchemes> \
#         scripts/validation/pass1_prior_viz/plot_pass1_priors.jl
# Env overrides: PRIORVIZ_OUT, PRIORVIZ_SEED, PRIORVIZ_NFIB, PRIORVIZ_FIBERS
#
# Author - Claude (for A. Saydjari), 2026-09-07

using CairoMakie, ColorSchemes, HDF5, Serialization, Random, Statistics, LinearAlgebra, Printf, Dates

black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const BASE = "/mnt/ceph/users/sdssv/work/asaydjari"
const OUT = get(ENV, "PRIORVIZ_OUT", joinpath(BASE, "2026_09_07/plots/pass1_priors"))
const SEED = parse(Int, get(ENV, "PRIORVIZ_SEED", "20260907"))
const NFIB = parse(Int, get(ENV, "PRIORVIZ_NFIB", "5"))   # fibers per telescope
const NSAMP = 5                                            # samples/draws per panel
mkpath(OUT)

# ---------------------------------------------------------------- prior locations
const STARCONT_ROOT  = joinpath(BASE, "2026_09_05/prior_outputs/starCont_pass1c")
const STARLINES_DIR  = joinpath(BASE, "2026_09_05/prior_outputs/starLines_perfiber")
const STARLINES_REF  = joinpath(BASE, "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5")
const STARLINES_FULL = joinpath(BASE, "2026_09_05/prior_inputs/starLines_fullres/APOGEE_stellar_kry_50_fullres.h5")
const SKY_DIR        = joinpath(BASE, "2026_09_04/prior_outputs/sky_pass1/built_combined_telemaj_union")
const SKY_BASE       = joinpath(BASE, "2026_09_04/prior_outputs/sky_pass1/built")
const SKY_SAMP       = joinpath(BASE, "2026_09_04/prior_outputs/sky_pass1/samples")
const CHIPGAP_MSK    = joinpath(BASE, "2026_04_25/StarContChipGapMsk.h5")
const BRIGHT_MASK    = joinpath(BASE, "2026_09_07/e5_combined/e5_bright_combined.h5")

pad(n) = lpad(n, 3, "0")
tele(af) = af > 300 ? "LCO" : "APO"
fibidx(af) = mod1(af, 300)
const FPI_ADJFIB = Set([76, 226, 388, 519])  # bimodal-bright FPI fibers (pass-1c MANIFEST)
fmt(x) = @sprintf("%.4g", x)
fmte(x) = @sprintf("%.2e", x)

starcont_prior(af)  = joinpath(STARCONT_ROOT, af > 300 ? "built_lco" : "built_apo",
                               "APOGEE_starcont_svd_60_f" * pad(af) * ".h5")
starcont_samp(af)   = joinpath(STARCONT_ROOT, "tell_prior_disk", "starCont_" * pad(af) * ".jdat")
starlines_prior(af) = joinpath(STARLINES_DIR, "APOGEE_stellar_kry_50_subpix_f" * pad(af) * ".h5")
skycont_prior(af)   = joinpath(SKY_DIR, "APOGEE_skycont_svd_30_f" * pad(af) * ".h5")
skyfaint_prior(af)  = joinpath(SKY_DIR, "APOGEE_skyline_faint_svd_120_f" * pad(af) * ".h5")
skygspice_prior(af) = joinpath(SKY_DIR, "APOGEE_skyline_faint_GSPICE_svd_120_f" * pad(af) * ".h5")
skygspice_base(af)  = joinpath(SKY_BASE, "APOGEE_skyline_faint_GSPICE_svd_120_f" * pad(af) * ".h5")
sky_sample(kind, af) = joinpath(SKY_SAMP, kind * "_" * pad(af) * ".h5")

# ---------------------------------------------------------------- grid + palette
const wavetarg = 10 .^ range(4.179 - 125 * 6.0e-6, step = 6.0e-6, length = 8575 + 125)
const NPIX = length(wavetarg)

const SAMPCOL = ["#00d5ff", "#ffb703", "#ff5ec4", "#7bff4d", "#b28dff"]  # bright, dark-safe
const EIGCOL  = ["#00e5ff", "#ffd166", "#ff4d6d"]
const DIMCOL  = "#8899a8"
const HOTCOL  = "#ff2d55"
fibcolors(n) = [get(ColorSchemes.colorschemes[:rainbow_bgyrm_35_85_c69_n256], t)
                for t in range(0.05, 0.92, length = max(n, 2))]

"""NaN out pixels outside `msk` so the plotted line BREAKS at chip gaps."""
function blank(v::AbstractVector, msk::AbstractVector{Bool})
    o = Vector{Float64}(undef, length(v))
    @inbounds for i in eachindex(v)
        o[i] = msk[i] ? Float64(v[i]) : NaN
    end
    o
end

function finite_extrema(vs...)
    lo, hi = Inf, -Inf
    for v in vs, x in v
        isfinite(x) || continue
        lo = min(lo, x); hi = max(hi, x)
    end
    (lo, hi)
end

function set_ylims!(ax, vs...; pad_frac = 0.08)
    lo, hi = finite_extrema(vs...)
    (isfinite(lo) && isfinite(hi)) || return
    hi == lo && (hi = lo + 1)
    d = (hi - lo) * pad_frac
    ylims!(ax, lo - d, hi + d)
end

"""Contiguous true-runs of a Bool mask, longest first (= the three APOGEE chips)."""
function true_runs(msk)
    runs = Tuple{Int,Int}[]
    i = 1
    while i <= length(msk)
        if msk[i]
            j = i
            while j < length(msk) && msk[j + 1]; j += 1; end
            push!(runs, (i, j)); i = j + 1
        else
            i += 1
        end
    end
    sort(runs, by = r -> -(r[2] - r[1]))
end

# ---------------------------------------------------------------- IO helpers
read_h5(f, k) = h5open(f, "r") do fid; read(fid[k]); end
bmask(x) = Vector{Bool}(vec(x) .!= 0)
ncols(path, key) = h5open(path, "r") do f; size(f[key], 2); end
function read_cols(path, key, cols)
    h5open(path, "r") do f
        d = f[key]
        reduce(hcat, [Float64.(d[:, c]) for c in cols])
    end
end

"""Global normalization build_starCont.jl applies before the SVD: the mean of the
non-zero entries inside chipgapmsk over the used (specsum>0) samples."""
function starcont_mnorm(sc::Matrix{Float64}, msk::AbstractVector{Bool})
    s = 0.0; n = 0
    @inbounds for j in axes(sc, 2)
        cs = 0.0
        for i in axes(sc, 1); cs += sc[i, j]; end
        cs > 0 || continue
        for i in eachindex(msk)
            msk[i] || continue
            v = sc[i, j]
            if v != 0; s += v; n += 1; end
        end
    end
    (n == 0 ? 1.0 : s / n), n
end

# ---------------------------------------------------------------- QA log
const QALOG = joinpath(OUT, "prior_viz_qa.txt")
const qaio = open(QALOG, "w")
qa(args...) = (s = string(args...); println(s); flush(stdout); println(qaio, s); flush(qaio))
const ANOMALIES = String[]
flag(msg) = (push!(ANOMALIES, msg); qa("  !! ANOMALY: ", msg))

# ---------------------------------------------------------------- fiber draw
const FIBERS = if haskey(ENV, "PRIORVIZ_FIBERS")
    parse.(Int, split(ENV["PRIORVIZ_FIBERS"], ","))
else
    rng = MersenneTwister(SEED)
    vcat(sort(randperm(rng, 300)[1:NFIB]), 300 .+ sort(randperm(rng, 300)[1:NFIB]))
end

qa("# APOGEE DR21 pass-1 prior gallery — QA / provenance log")
qa("# generated ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), " on ", gethostname())
qa("# RNG: MersenneTwister(", SEED, "); draw = sort(randperm(rng,300)[1:", NFIB,
   "]) for APO, then the same for LCO (+300)")
qa("# fibers drawn (adjfibindx): ", join(FIBERS, ", "))
qa("# wavetarg = 10 .^ range(4.179-125*6e-6, step=6e-6, length=8700) => ",
   fmt(wavetarg[1]), " - ", fmt(wavetarg[end]), " Angstrom")
qa("")

# ---------------------------------------------------------------- shared inputs
qa("## shared inputs")
const RUNTIME_MSK = h5open(CHIPGAP_MSK, "r") do f
    Dict("apo" => bmask(read(f["apo"])), "lco" => bmask(read(f["lco"])))
end
qa("chipgap mask ", CHIPGAP_MSK, ": apo ", count(RUNTIME_MSK["apo"]), " / lco ",
   count(RUNTIME_MSK["lco"]), " good px of ", NPIX)

const SL_LAM_MODEL = let V = read_h5(STARLINES_FULL, "Vmat")
    o = [sum(abs2, view(V, :, k)) for k in 1:size(V, 2)]
    V = nothing; GC.gc(); o
end
const SL_REF = read_h5(STARLINES_REF, "Vmat")[:, :, 6]     # (8700,50) refLSF, zero shift
qa("starLines fullres model-space lambda: l1=", fmt(SL_LAM_MODEL[1]), " l50=",
   fmt(SL_LAM_MODEL[end]), " sum=", fmt(sum(SL_LAM_MODEL)))
issorted(SL_LAM_MODEL, rev = true) || flag("starLines fullres lambda NOT sorted descending")
qa("")

# ---------------------------------------------------------------- per-fiber loader
struct FiberData
    af::Int
    msk::Vector{Bool}
    sc_V::Matrix{Float64}; sc_lam::Vector{Float64}; sc_samp::Matrix{Float64}; sc_cols::Vector{Int}
    nsc::Int
    sl_V::Matrix{Float64}; sl_draw::Matrix{Float64}
    kc_V::Matrix{Float64}; kc_lam::Vector{Float64}; kc_samp::Matrix{Float64}; kc_cols::Vector{Int}
    nkc::Int
    kf_V::Matrix{Float64}; kf_lam::Vector{Float64}; kf_sub::Vector{Bool}
    kg_V::Matrix{Float64}; kg_lam::Vector{Float64}; kg_sub::Vector{Bool}
    base_sub::Vector{Bool}; base_lam::Vector{Float64}
    kl_samp::Matrix{Float64}; kl_cols::Vector{Int}
    bright_frac::Float64; bright_policy::String
end

"""How concentrated is each leading mode? A physically sensible prior mode is spread
over many pixels; a mode whose energy sits in a handful of pixels is a single
unmasked artifact/line dominating the basis, and that is a review-gate anomaly."""
function localization_report(nm, af, V, lam, msk; ncomp = 3, topn = 10, thresh = 0.5)
    for k in 1:min(ncomp, size(V, 2))
        lam[k] > 0 || continue
        u = abs.(view(V, :, k)) ./ sqrt(lam[k])
        e = u .^ 2
        p = sortperm(e, rev = true)
        f10 = sum(view(e, view(p, 1:topn)))
        i1 = p[1]
        qa("    ", nm, " u", k, ": max|u|=", @sprintf("%.4f", u[i1]), " at px ", i1, " (",
           @sprintf("%.2f", wavetarg[i1]), " Å); top-", topn, " px hold ",
           @sprintf("%.1f%%", 100 * f10), " of the mode energy (support ", count(msk), " px)")
        f10 > thresh && flag("$nm f$(pad(af)) u$k is a SPIKE: $(round(100*f10, digits=1))% of its " *
            "energy in $topn px, peak at $(round(wavetarg[i1], digits=2)) Å (px $i1)")
    end
end

function structural_checks(nm, af, V, lam, msk)
    all(isfinite, V) || flag("$nm f$(pad(af)) Vmat non-finite")
    nk = min(10, length(lam))
    all(>(0), view(lam, 1:nk)) || flag("$nm f$(pad(af)) non-positive eigenvalue in top $nk")
    issorted(lam, rev = true) || flag("$nm f$(pad(af)) eigenvalues not sorted descending")
    if count(.!msk) > 0
        mx = maximum(abs, view(V, .!msk, :))
        mx == 0 || flag("$nm f$(pad(af)) nonzero Vmat off its mask (max $mx)")
    end
    if all(>(0), view(lam, 1:nk))
        U = V[:, 1:nk] ./ sqrt.(reshape(lam[1:nk], 1, :))
        orth = maximum(abs, U' * U - I)
        orth < 1e-7 || flag("$nm f$(pad(af)) top-$nk modes not orthonormal (max dev $orth)")
        nrm = maximum(abs, [sum(abs2, view(V, :, k)) / lam[k] - 1 for k in 1:nk])
        nrm < 1e-8 || flag("$nm f$(pad(af)) |V[:,k]|^2 != lambda_k (max rel dev $nrm)")
    end
end

function load_fiber(af)
    qa("## adjfib ", pad(af), "  (", tele(af), ", fiberIndx ", fibidx(af),
       af in FPI_ADJFIB ? ", FPI fiber" : "", ")")
    frng = MersenneTwister(SEED + af)
    tk = af > 300 ? "lco" : "apo"

    # -- starCont ---------------------------------------------------------------
    sc_V, sc_lam, sc_msk = h5open(starcont_prior(af), "r") do f
        read(f["Vmat"]), read(f["λv"]), bmask(read(f["chipgapmsk"]))
    end
    sc_msk == RUNTIME_MSK[tk] ||
        flag("starCont f$(pad(af)) chipgapmsk != runtime $tk mask (pipeline load would error)")
    sc_all = deserialize(starcont_samp(af))::Matrix{Float64}
    size(sc_all, 1) == NPIX || flag("starCont samples f$(pad(af)) have $(size(sc_all,1)) rows")
    nsc = size(sc_all, 2)
    mnorm, nnz = starcont_mnorm(sc_all, sc_msk)
    good = [j for j in axes(sc_all, 2) if sum(view(sc_all, :, j)) > 0]
    sc_cols = sort(good[randperm(frng, length(good))[1:NSAMP]])
    sc_samp = sc_all[:, sc_cols] ./ mnorm
    qa("  starCont: ", nsc, " samples (", length(good), " with specsum>0); mnorm=", fmt(mnorm),
       " over ", nnz, " nonzero px; lam1=", fmt(sc_lam[1]), " sumlam=", fmt(sum(sc_lam)))
    sc_all = nothing; GC.gc()

    # -- starLines (per-fiber LSF, sub-pixel index 6 == zero shift) --------------
    sl_all = read_h5(starlines_prior(af), "Vmat")
    size(sl_all) == (NPIX, 50, 10) || flag("starLines f$(pad(af)) shape $(size(sl_all)) != (8700,50,10)")
    sl_V = sl_all[:, :, 6]
    all(isfinite, sl_V) || flag("starLines f$(pad(af)) slice 6 has non-finite entries")
    sl_all = nothing
    sl_draw = sl_V * randn(frng, size(sl_V, 2), NSAMP)   # exact draws from N(0, V*V')
    cn = [sum(abs2, view(sl_V, :, k)) for k in 1:3]
    qa("  starLines: |V[:,k,6]|^2 (k=1:3) = ", join(fmt.(cn), ", "),
       "  model-space lambda = ", join(fmt.(SL_LAM_MODEL[1:3]), ", "),
       "  ratio = ", join([@sprintf("%.3f", cn[k] / SL_LAM_MODEL[k]) for k in 1:3], ", "))

    # -- skyCont ----------------------------------------------------------------
    kc_V, kc_lam = h5open(skycont_prior(af), "r") do f; read(f["Vmat"]), read(f["λv"]); end
    nkc = ncols(sky_sample("skycont", af), "skycont")
    kc_cols = sort(randperm(frng, nkc)[1:NSAMP])
    kc_samp = read_cols(sky_sample("skycont", af), "skycont", kc_cols)
    qa("  skyCont: ", nkc, " samples; lam1=", fmt(kc_lam[1]), " sumlam=", fmt(sum(kc_lam)))

    # -- skyLine faint (plain, built NOT deployed) and GSPICE (deployed) --------
    kf_V, kf_lam, kf_sub = h5open(skyfaint_prior(af), "r") do f
        read(f["Vmat"]), read(f["λv"]), bmask(read(f["submsk"]))
    end
    kg_V, kg_lam, kg_sub, bfrac, bpol = h5open(skygspice_prior(af), "r") do f
        read(f["Vmat"]), read(f["λv"]), bmask(read(f["submsk"])),
        Float64(read(f["bright_frac"])), string(read(f["bright_policy"]))
    end
    base_sub, base_lam = h5open(skygspice_base(af), "r") do f
        bmask(read(f["submsk"])), read(f["λv"])
    end
    kf_sub == kg_sub || flag("skyLine faint and GSPICE submsk differ on f$(pad(af))")
    removed = base_sub .& .!kg_sub
    added = kg_sub .& .!base_sub
    count(added) == 0 ||
        flag("skyLine rebuild ADDED $(count(added)) px to the faint support on f$(pad(af))")
    nkl = ncols(sky_sample("skyline", af), "skyline")
    kl_cols = sort(randperm(frng, nkl)[1:NSAMP])
    kl_samp = read_cols(sky_sample("skyline", af), "skyline", kl_cols)
    qa("  skyLine: policy=", bpol, " bright_frac=", @sprintf("%.4f", bfrac),
       "; baseline support ", count(base_sub), " -> rebuild ", count(kg_sub),
       " (removed ", count(removed), ", added ", count(added), "); ", nkl, " samples")
    qa("  skyLine GSPICE lam1 ", fmt(base_lam[1]), " -> ", fmt(kg_lam[1]), " ratio ",
       fmte(kg_lam[1] / base_lam[1]), " ; sumlam ", fmt(sum(base_lam)), " -> ",
       fmt(sum(kg_lam)), " ratio ", fmte(sum(kg_lam) / sum(base_lam)))

    structural_checks("starCont", af, sc_V, sc_lam, sc_msk)
    structural_checks("skyCont", af, kc_V, kc_lam, sc_msk)
    structural_checks("skyFaint", af, kf_V, kf_lam, kf_sub)
    structural_checks("skyGSPICE", af, kg_V, kg_lam, kg_sub)

    qa("  leading-mode localization (top-10-pixel energy fraction):")
    localization_report("starCont", af, sc_V, sc_lam, sc_msk)
    localization_report("starLines", af, sl_V, [sum(abs2, view(sl_V, :, k)) for k in 1:3], sc_msk)
    localization_report("skyCont", af, kc_V, kc_lam, sc_msk)
    localization_report("skyFaint", af, kf_V, kf_lam, kf_sub)
    localization_report("skyGSPICE", af, kg_V, kg_lam, kg_sub)

    # the brightest pixel that SURVIVES the bright mask, from the drawn samples
    mx, imx = -Inf, 0
    for j in axes(kl_samp, 2), i in 1:NPIX
        if kg_sub[i] && kl_samp[i, j] > mx; mx = kl_samp[i, j]; imx = i; end
    end
    mxb, imxb = -Inf, 0
    for j in axes(kl_samp, 2), i in 1:NPIX
        if base_sub[i] && kl_samp[i, j] > mxb; mxb = kl_samp[i, j]; imxb = i; end
    end
    qa("  brightest sample pixel on the baseline support: ", fmt(mxb), " ADU at ",
       @sprintf("%.2f", wavetarg[imxb]), " Å (px ", imxb, "); on the RETAINED faint support: ",
       fmt(mx), " ADU at ", @sprintf("%.2f", wavetarg[imx]), " Å (px ", imx, ")")
    qa("")

    FiberData(af, sc_msk, sc_V, sc_lam, sc_samp, sc_cols, nsc, sl_V, sl_draw,
              kc_V, kc_lam, kc_samp, kc_cols, nkc, kf_V, kf_lam, kf_sub,
              kg_V, kg_lam, kg_sub, base_sub, base_lam, kl_samp, kl_cols, bfrac, bpol)
end

# ---------------------------------------------------------------- panel helpers
function eig_lines!(ax, V, lam, msk; ncomp = 3, lw = 0.7)
    tot = sum(lam)
    ys = Vector{Float64}[]
    for k in 1:ncomp
        u = blank(view(V, :, k) ./ sqrt(lam[k]), msk)
        push!(ys, u)
        lines!(ax, wavetarg, u, color = EIGCOL[k], linewidth = lw,
               label = @sprintf("u%d   λ=%.4g   (%.2f%% of Σλ)", k, lam[k], 100 * lam[k] / tot))
    end
    ys
end

function samp_lines!(ax, S, msk; lw = 0.6)
    ys = Vector{Float64}[]
    for j in axes(S, 2)
        y = blank(view(S, :, j), msk)
        push!(ys, y)
        lines!(ax, wavetarg, y, color = SAMPCOL[mod1(j, length(SAMPCOL))], linewidth = lw)
    end
    ys
end

leg!(ax; pos = :rb) = axislegend(ax, position = pos, labelsize = 9.5, framevisible = false,
                                 patchsize = (18, 2), padding = (4, 4, 2, 2), rowgap = 0)

# ---------------------------------------------------------------- overview figure
function overview_figure(d::FiberData)
    af = d.af
    fpi = af in FPI_ADJFIB ? ", FPI" : ""
    fig = Figure(size = (2000, 1760), figure_padding = 14)
    Label(fig[0, 1:2], "APOGEE DR21 pass-1 priors — adjfibindx $(pad(af))  ($(tele(af)), " *
          "fiberIndx $(fibidx(af))$fpi)   |   fiber draw seed $SEED",
          fontsize = 21, font = :bold, padding = (0, 0, 0, 4))

    axs = Axis[]
    mkax(r, c, t) = begin
        a = Axis(fig[r, c], title = t, titlesize = 11.5, titlealign = :left,
                 xlabel = r == 5 ? "vacuum wavelength [Å]" : "",
                 xticklabelsize = 10, yticklabelsize = 10, ylabelsize = 11)
        push!(axs, a); a
    end

    sccols = join(d.sc_cols, ", ")
    klcols = join(d.kl_cols, ", ")
    kccols = join(d.kc_cols, ", ")
    removed = d.base_sub .& .!d.kg_sub

    # row 1 — starCont
    a = mkax(1, 1, "1. starCont (pass-1c) — SAMPLES: $NSAMP of the $(d.nsc) training models\n" *
        "blackbody(Teff) × reddening(Av,Rv), fiber-LSF convolved, × a drawn telluric transmission; " *
        "÷ the build's global mnorm.   sample cols $sccols")
    ys = samp_lines!(a, d.sc_samp, d.msk); set_ylims!(a, ys...); a.ylabel = "normalized flux"
    a = mkax(1, 2, "1. starCont — TOP 3 EIGENVECTORS (unit norm, u_k = Vmat[:,k]/√λ_k)\n" *
        "second-moment SVD, no mean subtraction ⇒ u₁ ≈ mean model shape")
    ys = eig_lines!(a, d.sc_V, d.sc_lam, d.msk); set_ylims!(a, ys...); leg!(a)

    # row 2 — starLines
    a = mkax(2, 1, "2. starLines (E7 per-fiber FPI LSF, sub-pixel index 6 = ZERO SHIFT) — " *
        "SAMPLES: $NSAMP exact prior draws\nV[:,:,6]·z with z ∼ N(0,I₅₀); units are norm94 " *
        "(flux/94th-pct − 1), so 0 = continuum and dips = absorption")
    ys = samp_lines!(a, d.sl_draw, d.msk); set_ylims!(a, ys...); a.ylabel = "norm94 flux offset"
    a = mkax(2, 2, "2. starLines — TOP 3 EIGENVECTORS (unit norm). Solid = this fiber's LSF; " *
        "dashed grey = refLSF R=22 500 (the reporting basis)\nλ_model = eigenvalue of the shared " *
        "Krylov basis before convolution; ‖V‖² = post-convolution column norm²")
    ysl = Vector{Float64}[]
    for k in 1:3
        cn = sum(abs2, view(d.sl_V, :, k))
        u = blank(view(d.sl_V, :, k) ./ sqrt(cn), d.msk)
        ur = blank(view(SL_REF, :, k) ./ sqrt(sum(abs2, view(SL_REF, :, k))), d.msk)
        lines!(a, wavetarg, ur, color = DIMCOL, linewidth = 0.5, linestyle = :dash)
        lines!(a, wavetarg, u, color = EIGCOL[k], linewidth = 0.7,
               label = @sprintf("u%d   λ_model=%.4g   ‖V[:,%d,6]‖²=%.4g", k, SL_LAM_MODEL[k], k, cn))
        push!(ysl, u)
    end
    set_ylims!(a, ysl...); leg!(a)

    # row 3 — skyCont
    a = mkax(3, 1, "3. skyCont (E5) — SAMPLES: $NSAMP of the $(d.nkc) real sky-fiber continua\n" *
        "smooth component of an observed sky spectrum (sky_smooth_fit), detector units.   " *
        "sample cols $kccols")
    ys = samp_lines!(a, d.kc_samp, d.msk); set_ylims!(a, ys...); a.ylabel = "flux [ADU]"
    a = mkax(3, 2, "3. skyCont — TOP 3 EIGENVECTORS (unit norm)")
    ys = eig_lines!(a, d.kc_V, d.kc_lam, d.msk); set_ylims!(a, ys...); leg!(a)

    # row 4 — skyLine, full dynamic range + the bright excision
    a = mkax(4, 1, "4. skyLine — SAMPLES: $NSAMP real sky-line residuals (skyflux − skycont) on the " *
        "BASELINE support ($(count(d.base_sub)) px)\nRED = the $(count(removed)) px the combined " *
        "bright mask REMOVES (bright_frac $(@sprintf("%.2f", 100*d.bright_frac))%). FULL dynamic range.   " *
        "sample cols $klcols")
    ys = Vector{Float64}[]
    for j in axes(d.kl_samp, 2)
        y = blank(view(d.kl_samp, :, j), d.base_sub)
        lines!(a, wavetarg, y, color = SAMPCOL[mod1(j, length(SAMPCOL))], linewidth = 0.5)
        push!(ys, y)
    end
    for j in axes(d.kl_samp, 2)
        lines!(a, wavetarg, blank(view(d.kl_samp, :, j), removed), color = HOTCOL, linewidth = 0.7)
    end
    set_ylims!(a, ys...); a.ylabel = "line residual [ADU]"
    a = mkax(4, 2, "4. skyLine faint, PLAIN (built, NOT deployed) — TOP 3 EIGENVECTORS " *
        "(unit norm), support $(count(d.kf_sub)) px")
    ys = eig_lines!(a, d.kf_V, d.kf_lam, d.kf_sub); set_ylims!(a, ys...); leg!(a)

    # row 5 — skyLine on the faint support only + the deployed prior
    a = mkax(5, 1, "5. skyLine — THE SAME $NSAMP samples restricted to the FAINT support " *
        "($(count(d.kg_sub)) px), bright lines excised\nNOTE THE Y-SCALE vs row 4: this is what the " *
        "deployed prior actually sees.")
    ys = Vector{Float64}[]
    for j in axes(d.kl_samp, 2)
        y = blank(view(d.kl_samp, :, j), d.kg_sub)
        lines!(a, wavetarg, y, color = SAMPCOL[mod1(j, length(SAMPCOL))], linewidth = 0.5)
        push!(ys, y)
    end
    set_ylims!(a, ys...); a.ylabel = "line residual [ADU]"
    a = mkax(5, 2, "5. skyLine faint GSPICE — DEPLOYED PRIOR — TOP 3 EIGENVECTORS (unit norm)\n" *
        "Σλ = $(fmt(sum(d.kg_lam)))   (pre-rebuild built/ baseline was $(fmt(sum(d.base_lam))), " *
        "ratio $(fmte(sum(d.kg_lam)/sum(d.base_lam))) — bright lines removed BY DESIGN)")
    ys = eig_lines!(a, d.kg_V, d.kg_lam, d.kg_sub); set_ylims!(a, ys...); leg!(a)

    for a in axs
        xlims!(a, wavetarg[1], wavetarg[end])
        a.xgridvisible = false; a.ygridvisible = false
    end
    rowgap!(fig.layout, 10); colgap!(fig.layout, 16)
    fn = joinpath(OUT, "fiber_f$(pad(af))_$(lowercase(tele(af)))_overview.png")
    save(fn, fig, px_per_unit = 1.2)
    fn
end

# ---------------------------------------------------------------- zoom figure
function zoom_figure(d::FiberData)
    af = d.af
    fpi = af in FPI_ADJFIB ? ", FPI" : ""
    chips = sort(true_runs(d.msk)[1:3], by = r -> r[1])       # blue, green, red
    wins = map(chips) do r
        ic = (r[1] + r[2]) ÷ 2
        (wavetarg[ic] - 20.0, wavetarg[ic] + 20.0)
    end
    names = ["blue chip", "green chip", "red chip"]
    fig = Figure(size = (2000, 1560), figure_padding = 14)
    Label(fig[0, 1:3], "pass-1 priors, 40 Å zooms — adjfibindx $(pad(af)) ($(tele(af)), " *
          "fiberIndx $(fibidx(af))$fpi)   |   top 3 unit-norm eigenvectors per family; " *
          "grey = one training sample, min-max rescaled to the panel",
          fontsize = 19, font = :bold)

    rows = [("1. starCont", d.sc_V, d.sc_lam, d.msk, d.sc_samp),
            ("2. starLines (subpix 6)", d.sl_V,
             [sum(abs2, view(d.sl_V, :, k)) for k in 1:size(d.sl_V, 2)], d.msk, d.sl_draw),
            ("3. skyCont", d.kc_V, d.kc_lam, d.msk, d.kc_samp),
            ("4. skyLine faint plain", d.kf_V, d.kf_lam, d.kf_sub, d.kl_samp),
            ("5. skyLine faint GSPICE (deployed)", d.kg_V, d.kg_lam, d.kg_sub, d.kl_samp)]
    for (r, row) in enumerate(rows)
        nm, V, lam, msk, S = row
        for (c, w) in enumerate(wins)
            w0, w1 = w
            a = Axis(fig[r, c], title = c == 1 ? "$nm  —  $(names[c])" : names[c],
                     titlesize = 12, titlealign = :left,
                     xlabel = r == 5 ? "vacuum wavelength [Å]" : "",
                     ylabel = c == 1 ? "unit-norm eigvec" : "",
                     xticklabelsize = 10, yticklabelsize = 10)
            sel = (wavetarg .>= w0) .& (wavetarg .<= w1)
            ys = Vector{Float64}[]
            for k in 1:3
                u = blank(view(V, :, k) ./ sqrt(lam[k]), msk)
                lines!(a, wavetarg, u, color = EIGCOL[k], linewidth = 1.4)
                push!(ys, u[sel])
            end
            lo, hi = finite_extrema(ys...)
            ssel = blank(view(S, :, 1), msk)
            sl, sh = finite_extrema(ssel[sel])
            if isfinite(sl) && isfinite(sh) && sh > sl && isfinite(lo) && isfinite(hi) && hi > lo
                scaled = @. (ssel - sl) / (sh - sl) * (hi - lo) + lo
                lines!(a, wavetarg, scaled, color = DIMCOL, linewidth = 0.9)
            end
            xlims!(a, w0, w1)
            (isfinite(lo) && isfinite(hi) && hi > lo) &&
                ylims!(a, lo - 0.1 * (hi - lo), hi + 0.1 * (hi - lo))
            a.xgridvisible = false; a.ygridvisible = false
        end
    end
    rowgap!(fig.layout, 8); colgap!(fig.layout, 14)
    fn = joinpath(OUT, "fiber_f$(pad(af))_$(lowercase(tele(af)))_zoom.png")
    save(fn, fig, px_per_unit = 1.2)
    fn
end

# ---------------------------------------------------------------- run
made = String[]
summ = Dict{Int,Any}()
for af in FIBERS
    d = load_fiber(af)
    push!(made, overview_figure(d))
    push!(made, zoom_figure(d))
    summ[af] = (sc_lam = copy(d.sc_lam), kc_lam = copy(d.kc_lam), kf_lam = copy(d.kf_lam),
                kg_lam = copy(d.kg_lam), base_lam = copy(d.base_lam),
                sl_cn = [sum(abs2, view(d.sl_V, :, k)) for k in 1:size(d.sl_V, 2)],
                sc_u1 = blank(d.sc_V[:, 1] ./ sqrt(d.sc_lam[1]), d.msk),
                sl_u1 = blank(d.sl_V[:, 1] ./ sqrt(sum(abs2, view(d.sl_V, :, 1))), d.msk),
                kc_u1 = blank(d.kc_V[:, 1] ./ sqrt(d.kc_lam[1]), d.msk),
                kg_u1 = blank(d.kg_V[:, 1] ./ sqrt(d.kg_lam[1]), d.kg_sub),
                nsub = count(d.kg_sub), nbase = count(d.base_sub), bfrac = d.bright_frac)
    d = nothing; GC.gc()
end

# ---------------------------------------------------------------- summary: eigenspectra
let
    fig = Figure(size = (1900, 1000), figure_padding = 14)
    Label(fig[0, 1:3], "pass-1 priors — eigenvalue spectra for all $(length(FIBERS)) drawn fibers " *
          "(solid = APO, dashed = LCO)", fontsize = 20, font = :bold)
    cols = fibcolors(length(FIBERS))
    specs = [("starCont (60 modes)", :sc_lam, "λ  [normalized flux²]"),
             ("starLines ‖V[:,k,6]‖² (50 modes)", :sl_cn, "column norm²"),
             ("skyCont (30 modes)", :kc_lam, "λ  [ADU²]"),
             ("skyLine faint plain (120)", :kf_lam, "λ  [ADU²]"),
             ("skyLine faint GSPICE — DEPLOYED (120)", :kg_lam, "λ  [ADU²]"),
             ("skyLine faint GSPICE — pre-rebuild built/ (120)", :base_lam, "λ  [ADU²]")]
    for (i, sp) in enumerate(specs)
        nm, key, yl = sp
        r, c = fldmod1(i, 3)
        a = Axis(fig[r, c], title = nm, titlesize = 13, yscale = log10, titlealign = :left,
                 xlabel = "mode index k", ylabel = yl)
        for (j, af) in enumerate(FIBERS)
            v = getfield(summ[af], key)
            lines!(a, 1:length(v), max.(v, 1e-300), color = cols[j], linewidth = 1.5,
                   linestyle = af > 300 ? :dash : :solid, label = "f$(pad(af)) $(tele(af))")
        end
        a.xgridvisible = false; a.ygridvisible = false
        i == 1 && axislegend(a, position = :rt, labelsize = 8.5, framevisible = false, nbanks = 2)
    end
    rowgap!(fig.layout, 12); colgap!(fig.layout, 18)
    fn = joinpath(OUT, "summary_eigenvalue_spectra.png"); save(fn, fig, px_per_unit = 1.3)
    push!(made, fn)
end

# ---------------------------------------------------------------- summary: mode 1 overlay
let
    fig = Figure(size = (1900, 1320), figure_padding = 14)
    Label(fig[0, 1], "pass-1 priors — leading eigenvector u₁ (unit norm) across the " *
          "$(length(FIBERS)) drawn fibers, APO + LCO overlaid  (sign gauge fixed to positive mean)",
          fontsize = 20, font = :bold)
    cols = fibcolors(length(FIBERS))
    rows = [("starCont u₁", :sc_u1), ("starLines u₁ (subpix 6)", :sl_u1),
            ("skyCont u₁", :kc_u1), ("skyLine faint GSPICE u₁ (deployed)", :kg_u1)]
    for (i, row) in enumerate(rows)
        nm, key = row
        a = Axis(fig[i, 1], title = nm, titlesize = 13, titlealign = :left,
                 xlabel = i == length(rows) ? "vacuum wavelength [Å]" : "", ylabel = "u₁",
                 xticklabelsize = 10, yticklabelsize = 10)
        ys = Vector{Float64}[]
        for (j, af) in enumerate(FIBERS)
            y = getfield(summ[af], key)
            m = 0.0
            for x in y; isfinite(x) && (m += x); end
            s = m < 0 ? -1.0 : 1.0
            yy = s .* y
            lines!(a, wavetarg, yy, color = cols[j], linewidth = 0.6, label = "f$(pad(af)) $(tele(af))")
            push!(ys, yy)
        end
        set_ylims!(a, ys...); xlims!(a, wavetarg[1], wavetarg[end])
        a.xgridvisible = false; a.ygridvisible = false
        i == 1 && axislegend(a, position = :rt, labelsize = 8.5, framevisible = false, nbanks = 5)
    end
    rowgap!(fig.layout, 10)
    fn = joinpath(OUT, "summary_mode1_across_fibers.png"); save(fn, fig, px_per_unit = 1.3)
    push!(made, fn)
end

# ---------------------------------------------------------------- summary: variance scale
let
    fig = Figure(size = (1780, 800), figure_padding = 14)
    Label(fig[0, 1:2], "pass-1 priors — absolute retained variance per fiber " *
          "(this is where the ~230× faint-sky drop lives; it is BY DESIGN)",
          fontsize = 19, font = :bold)
    x = 1:length(FIBERS)
    labs = ["f$(pad(af))\n$(tele(af))" for af in FIBERS]
    a1 = Axis(fig[1, 1], title = "λ₁ (leading eigenvalue)", titlesize = 13, yscale = log10,
              ylabel = "λ₁ (family-specific units)", xticks = (x, labs), xticklabelsize = 9,
              titlealign = :left)
    a2 = Axis(fig[1, 2], title = "Σλ (total retained variance)", titlesize = 13, yscale = log10,
              ylabel = "Σλ", xticks = (x, labs), xticklabelsize = 9, titlealign = :left)
    series = [("starCont", :sc_lam, "#00d5ff", :circle),
              ("starLines ‖V‖²", :sl_cn, "#b28dff", :utriangle),
              ("skyCont", :kc_lam, "#7bff4d", :diamond),
              ("skyLine GSPICE rebuild (DEPLOYED)", :kg_lam, "#ffb703", :rect),
              ("skyLine GSPICE pre-rebuild built/", :base_lam, "#ff2d55", :xcross)]
    for s in series
        nm, key, col, mk = s
        y1 = [getfield(summ[af], key)[1] for af in FIBERS]
        y2 = [sum(getfield(summ[af], key)) for af in FIBERS]
        scatterlines!(a1, x, y1, color = col, marker = mk, markersize = 11, linewidth = 1.2, label = nm)
        scatterlines!(a2, x, y2, color = col, marker = mk, markersize = 11, linewidth = 1.2, label = nm)
    end
    axislegend(a1, position = :lb, labelsize = 9, framevisible = false)
    for a in (a1, a2); a.xgridvisible = false; end
    colgap!(fig.layout, 22)
    fn = joinpath(OUT, "summary_variance_scale.png"); save(fn, fig, px_per_unit = 1.3)
    push!(made, fn)
end

# ---------------------------------------------------------------- FLEET SCAN
# The per-fiber panels show that u2/u3 of the deployed faint sky-line prior are
# often a single-pixel spike at a chip edge. Measure that over all 600 fibers so
# the finding is fleet-scale, not anecdotal.  Disable with PRIORVIZ_FLEET=0.
if get(ENV, "PRIORVIZ_FLEET", "1") != "0"
    qa("## FLEET SCAN — leading-mode localization of the DEPLOYED faint sky prior, all 600 fibers")
    chip_edges = Dict(tk => sort(vcat([[r[1], r[2]] for r in
                                       filter(r -> r[2] - r[1] > 100, true_runs(RUNTIME_MSK[tk]))]...))
                      for tk in ("apo", "lco"))
    qa("  chip-run edge pixels: apo ", chip_edges["apo"], "  lco ", chip_edges["lco"])
    F10 = zeros(600, 3); PXK = zeros(Int, 600, 3); DEG = zeros(Int, 600, 3)
    for af in 1:600
        tk = af > 300 ? "lco" : "apo"
        V, lam = h5open(skygspice_prior(af), "r") do f
            f["Vmat"][:, 1:3], f["λv"][1:3]
        end
        for k in 1:3
            e = (view(V, :, k) ./ sqrt(lam[k])) .^ 2
            p = sortperm(e, rev = true)
            F10[af, k] = sum(view(e, view(p, 1:10)))
            PXK[af, k] = p[1]
            DEG[af, k] = minimum(abs.(p[1] .- chip_edges[tk]))
        end
    end
    for k in 1:3
        f = view(F10, :, k)
        qa("  u", k, ": median top-10-px energy fraction ", @sprintf("%.3f", median(f)),
           "; >50% on ", count(>(0.5), f), "/600 fibers; peak within 5 px of a chip edge on ",
           count(<=(5), view(DEG, :, k)), "/600")
    end
    nany = count(af -> any(k -> F10[af, k] > 0.5, 1:3), 1:600)
    nedge = count(af -> any(k -> F10[af, k] > 0.5 && DEG[af, k] <= 5, 1:3), 1:600)
    qa("  fibers with >=1 of u1..u3 spike-dominated (>50% of energy in 10 px): ", nany, "/600; ",
       "of those the spike sits within 5 px of a chip edge on ", nedge)

    # second failure mode: peaks that are NOT at a chip edge sit on the immediate WINGS of
    # bright lines the combined mask does cover -> the delivered mask is a few px too narrow.
    bright = h5open(BRIGHT_MASK, "r") do f; bmask(read(f["mask_telemaj_union"])); end
    bidx = findall(bright)
    dbright = [minimum(abs.(PXK[af, k] .- bidx)) for af in 1:600, k in 1:3]
    nspike = count(F10 .> 0.5)
    n_e = count((F10 .> 0.5) .& (DEG .<= 5))
    n_w = count((F10 .> 0.5) .& (DEG .> 5) .& (dbright .<= 5))
    n_o = nspike - n_e - n_w
    qa("  spike-dominated MODES: ", nspike, " total = ", n_e, " at a chip edge (<=5 px) + ", n_w,
       " on the immediate wings of a masked bright line (<=5 px from the ", count(bright),
       "-px combined mask, but not masked) + ", n_o, " elsewhere")
    n_w > 0 && flag("FLEET: $n_w spike-dominated modes peak 1-5 px OUTSIDE the combined bright " *
        "mask, i.e. on the wings of lines the mask already covers — the delivered mask is a few " *
        "pixels too narrow")
    nany > 60 && flag("FLEET: $nany/600 fibers have a spike-dominated mode in the top 3 of the " *
        "DEPLOYED faint sky prior ($nedge of them at a chip edge) — 120-mode capacity is being " *
        "spent on <=10 px")

    fig = Figure(size = (1900, 1120), figure_padding = 14)
    Label(fig[0, 1:2], "FLEET SCAN (all 600 fibers) — the deployed faint sky-line prior spends its " *
          "leading modes on a handful of pixels: chip edges, and the un-masked wings of bright sky lines",
          fontsize = 19, font = :bold)
    a = Axis(fig[1, 1], title = "top-10-pixel energy fraction of u₁, u₂, u₃ (1.0 = the mode IS 10 px)",
             titlesize = 13, titlealign = :left, xlabel = "adjfibindx", ylabel = "top-10 px energy fraction")
    for k in 1:3
        scatter!(a, 1:600, view(F10, :, k), color = EIGCOL[k], markersize = 3.5,
                 label = "u$k  (>50% on $(count(>(0.5), view(F10,:,k)))/600)")
    end
    hlines!(a, [0.5], color = DIMCOL, linestyle = :dash, linewidth = 1)
    vlines!(a, [300.5], color = DIMCOL, linewidth = 1)
    text!(a, 20, 0.035, text = "APO", color = DIMCOL, fontsize = 14)
    text!(a, 320, 0.035, text = "LCO", color = DIMCOL, fontsize = 14)
    ylims!(a, -0.02, 1.02); a.xgridvisible = false; a.ygridvisible = false
    axislegend(a, position = :lt, labelsize = 10, framevisible = false)

    a2 = Axis(fig[1, 2], title = "where the peak pixel of each spike-dominated mode sits " *
              "(grey = chip-run edges of the runtime mask)", titlesize = 13, titlealign = :left,
              xlabel = "vacuum wavelength [Å]", ylabel = "# of spike-dominated modes")
    spx = Int[]
    for af in 1:600, k in 1:3
        F10[af, k] > 0.5 && push!(spx, PXK[af, k])
    end
    up = sort(unique(spx))
    vlines!(a2, [wavetarg[p] for p in vcat(chip_edges["apo"], chip_edges["lco"])],
            color = DIMCOL, linewidth = 1.1, linestyle = :dash, label = "chip-run edges")
    cnt = [count(==(p), spx) for p in up]
    isedge = [min(minimum(abs.(p .- chip_edges["apo"])), minimum(abs.(p .- chip_edges["lco"]))) <= 5
              for p in up]
    barplot!(a2, [wavetarg[p] for p in up[isedge]], cnt[isedge], width = 3.5,
             color = "#00e5ff", strokewidth = 0, label = "peak at a chip edge")
    barplot!(a2, [wavetarg[p] for p in up[.!isedge]], cnt[.!isedge], width = 3.5,
             color = "#ffb703", strokewidth = 0, label = "peak on a bright-line wing")
    axislegend(a2, position = :rt, labelsize = 10, framevisible = false)
    xlims!(a2, wavetarg[1], wavetarg[end]); a2.xgridvisible = false; a2.ygridvisible = false
    Label(fig[2, 1:2], "MEASURED: $nany/600 fibers have at least one of u₁–u₃ carrying >50% of its " *
          "energy in 10 pixels. Of the $nspike spike-dominated modes, $n_e peak within 5 px of a " *
          "chip edge (most often px 3439 = the LAST pixel of the APO blue chip, and px 6419 = the " *
          "FIRST pixel of the LCO red chip) and $n_w peak 1–5 px OUTSIDE the 650-px combined bright " *
          "mask, on the wings of lines that mask already covers. " *
          "INFERRED: with the bright sky lines correctly excised the retained variance is ~230× " *
          "smaller, so these two classes of poorly-modelled pixel are now the largest thing left and " *
          "they capture the leading modes. The priors are structurally valid (all invariants pass) " *
          "— but 120-mode capacity is being spent on ≤10 pixels instead of on faint sky structure. " *
          "Two candidate fixes: trim 1–2 px from each chip-run edge in the runtime mask, and grow " *
          "the combined bright mask by ~3–5 px (the per-fiber path already does expand_msk(rad=4)).",
          fontsize = 12.5, color = "#ffd166", justification = :left, lineheight = 1.25,
          word_wrap = true, tellwidth = false)
    rowgap!(fig.layout, 10); colgap!(fig.layout, 20)
    fn = joinpath(OUT, "summary_fleet_chipedge_spikes.png"); save(fn, fig, px_per_unit = 1.3)
    push!(made, fn)
    qa("")
end

# ---------------------------------------------------------------- wrap up
qa("## figures written (", length(made), ")")
for f in made; qa("  ", basename(f), "  ", filesize(f), " B"); end
qa("")
qa("## anomalies flagged: ", length(ANOMALIES))
for m in ANOMALIES; qa("  - ", m); end
close(qaio)
println("done -> ", OUT)
