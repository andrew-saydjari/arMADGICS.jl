# build.jl -- orchestration: run the stages in order, assemble the columns, write one file.

"Everything the builder needs to know, from the CLI."
Base.@kwdef struct BuildConfig
    full_list::String
    almanac::String
    apred::String
    workup::String
    lists::String                 # directory with bad_nights.csv, bad_fibers.csv, rv_railed.csv
    fiber_offsets::String = ""    # fiber_offsets.csv; "" => omit the three v_rad_corr columns
    offtarget::Bool = true
    out::String
    report::String = ""
    cache::String                 # cache directory for identity / times / kernel outputs
    iers::String                  # pinned IERS-A table (finals2000A format); required, never downloaded
    iers_b::String = ""           # pinned IERS-B table (eopc04) substituted as astropy IERS_Auto does
    bary::Symbol = :julia         # :julia (bary.jl) or :astropy (Python kernel, cross-check only)
    python::String = "python3"    # only for numerics legacy040 or bary astropy
    kernels::String = joinpath(@__DIR__, "..", "python", "legacy_kernels.py")
    numerics::Numerics = JuliaNumerics()
end

"Path recorded in the product: relative to scripts/release_tables when the file is vendored there."
function _recorded_path(p::AbstractString)
    root = rstrip(normpath(joinpath(@__DIR__, "..", "..")), '/')
    a = abspath(p)
    startswith(a, root * "/") ? relpath(a, root) : a
end

"Refuse a file whose sha256 differs from a SHA256SUMS entry beside it (vendored tables)."
function _check_sum(path::AbstractString, sha::AbstractString)
    sums = joinpath(dirname(abspath(path)), "SHA256SUMS")
    isfile(sums) || return nothing
    for l in eachline(sums)
        f = split(l)
        length(f) == 2 && f[2] == basename(path) && f[1] != sha &&
            error("IERS file $path: sha256 $sha does not match $sums ($(f[1]))")
    end
    return nothing
end

"Load the pinned IERS tables (A, and B unless empty) and check the exposure times; loud warnings, never silent."
function check_iers(path::AbstractString, path_b::AbstractString, mjd_mid::AbstractVector{Float64})
    for p in (path, path_b)
        isempty(p) || isfile(p) || error("IERS file not found: $p (pass --iers / --iers-b; nothing is downloaded)")
    end
    isempty(path_b) && @warn "IERS: no IERS-B file given; Bulletin B columns of $path are used as is (NOT what 0.4.0 did)"
    T = read_finals2000A(path; iers_b = path_b)
    _check_sum(path, T.sha256)
    isempty(path_b) || _check_sum(path_b, T.sha256_b)
    cov = eop_coverage(T, mjd_mid)
    warn = String[]
    cov.n_before_table + cov.n_after_table > 0 &&
        push!(warn, "$(cov.n_before_table + cov.n_after_table) exposure times are OUTSIDE the IERS table " *
                    "($(T.first_mjd)-$(T.last_mjd)); astropy would hold the end values")
    cov.n_predicted > 0 &&
        push!(warn, "$(cov.n_predicted) exposure times are at or after MJD $(T.predictive_mjd), where the IERS " *
                    "table holds PREDICTED (not measured) Earth orientation; use a newer table")
    for w in warn
        @warn "IERS: " * w
    end
    say(@sprintf("IERS-A %s sha256 %s: MJD %.0f-%.0f, measured (non-predicted) before %.0f, Bulletin B to %.0f",
                 path, T.sha256, T.first_mjd, T.last_mjd, T.predictive_mjd, T.last_bulletin_b_mjd))
    say(@sprintf("IERS-B %s sha256 %s: substituted up to MJD %.0f; exposure times %.3f-%.3f; warnings: %d",
                 path_b, T.sha256_b, T.last_iers_b_mjd, cov.min_mjd, cov.max_mjd, length(warn)))
    rep = Dict{String,Any}("path" => abspath(path), "recorded_path" => _recorded_path(path), "sha256" => T.sha256,
        "bytes" => T.nbytes, "path_b" => isempty(path_b) ? "" : abspath(path_b),
        "recorded_path_b" => isempty(path_b) ? "" : _recorded_path(path_b), "sha256_b" => T.sha256_b,
        "first_mjd" => T.first_mjd, "last_mjd" => T.last_mjd,
        "predictive_mjd" => T.predictive_mjd, "last_bulletin_b_mjd" => T.last_bulletin_b_mjd,
        "last_iers_b_mjd" => T.last_iers_b_mjd,
        "coverage" => Dict(string(k) => v for (k, v) in pairs(cov)), "warnings" => warn)
    return T, rep
end

erfa_version() = unsafe_string(ccall((:eraVersion, LIBERFA), Cstring, ()))

say(s) = (println(stdout, "[$(Dates.format(now(), "HH:MM:SS"))] ", s); flush(stdout))

"Run a stage, report wall time, return its value."
macro stage(name, ex)
    quote
        local t0 = time()
        say("stage " * $(esc(name)) * " ...")
        local v = $(esc(ex))
        say(@sprintf("stage %s done in %.1f s", $(esc(name)), time() - t0))
        v
    end
end

function _run_kernel(cfg::BuildConfig, args::Vector{String}, out::String)
    isfile(out) && return out
    tmp = out * ".tmp"
    run(pipeline(`$(cfg.python) -I $(cfg.kernels) $args --out $tmp`; stdout = stdout, stderr = stderr))
    mv(tmp, out; force = true)
    return out
end

read_workup(dir, key) = h5open(f -> read(f[key]), joinpath(dir, "arMADGICS_out_$(key).h5"), "r")

function build(cfg::BuildConfig)
    T0 = time()
    mkpath(cfg.cache)
    timing = Dict{String,Float64}()
    tick(k, t) = (timing[k] = time() - t)

    t = time(); id = @stage "identity" load_identity(cfg.full_list, joinpath(cfg.cache, "identity.h5")); tick("identity", t)
    N = length(id); say("rows: $N")

    tpath = joinpath(cfg.cache, "exposure_times.h5")
    t = time()
    et = @stage "times" begin
        if isfile(tpath)
            load_times(tpath)
        else
            x = harvest_times(id, cfg.apred); save_times(tpath, x); x
        end
    end
    tick("times", t)

    # barycentric basis per exposure, from the pinned IERS table: Julia (default) or the astropy
    # kernel (cross-check); numpy pow10 kernel only for legacy040. Python kernels run while the
    # almanac join reads in parallel; their outputs are cached under names that carry the table hash.
    T_eop, iers_rep = check_iers(cfg.iers, cfg.iers_b, et.mjd_mid_precise[et.status .== "ok"])
    pow_path = joinpath(cfg.cache, "pow10.h5")
    basis_path = joinpath(cfg.cache, "bary_basis_astropy_$(first(T_eop.sha256, 12))_$(first(T_eop.sha256_b, 12)).h5")
    t = time()
    ktask = @async begin
        cfg.bary === :astropy &&
            _run_kernel(cfg, vcat(["bary-basis", "--times", tpath, "--iers-a", cfg.iers],
                                  isempty(cfg.iers_b) ? String[] : ["--iers-b", cfg.iers_b]), basis_path)
        cfg.numerics isa Legacy040Numerics && _run_kernel(cfg, ["pow10", "--workup", cfg.workup], pow_path)
        nothing
    end

    t = time(); C, fia, nnights = @stage "almanac join" join_almanac(id, cfg.almanac); tick("almanac", t)
    say("telescope-nights: $nnights")

    # mid-exposure time per row (any harvest status, as 0.4.0)
    t = time()
    tk = expkey.(et.tele, et.mjd, et.expnum)
    issorted(tk) || error("exposure_times not sorted")
    mjd_mid = fill(NaN, N)
    Threads.@threads for i in 1:N
        rk = expkey(id.tele[i], id.mjd[i], id.expnum[i])
        p = clamp(searchsortedfirst(tk, rk), 1, length(tk))
        tk[p] == rk && (mjd_mid[i] = et.mjd_mid_precise[p])
    end
    tick("mid_time", t)

    t = time(); wait(ktask); tick("kernels_wait", t)
    t = time()
    B = if cfg.bary === :astropy
        load_basis(basis_path)
    elseif cfg.bary === :julia
        ek, U = @stage "barycentric basis (Julia/ERFA)" bary_basis_julia(et, T_eop)
        BaryBasis(ek, U, "", Dict{String,Any}(tn => Dict("n_exposures" => length(ek[tc])) for (tc, tn) in TELENAME))
    else
        error("bary must be :julia or :astropy")
    end
    tick("bary_basis", t)
    ra = C.f64["ra"]; dec = C.f64["dec"]
    t = time(); v_bc = @stage "barycentric" barycentric(cfg.numerics, B, id, ra, dec); tick("bary", t)

    t = time()
    wk = Dict(k => read_workup(cfg.workup, k) for k in ("RV_pixoff_final", "RV_pix_var", "RV_flag", "ingestBit", "snr"))
    tick("workup_read", t)
    pix = Float64.(wk["RV_pixoff_final"]); var = Float64.(wk["RV_pix_var"])
    t = time()
    pw = pow10(cfg.numerics, pix, pow_path)
    v_rel, v_rad, e_v = @stage "velocities" velocities(pw, var, v_bc)
    tick("velocities", t)

    t = time(); bflag, nbad = @stage "bary_flag" bary_flag(id, et, v_bc, C.sdss_id_agrees, wk["RV_pixoff_final"]); tick("bary_flag", t)
    t = time(); glon, glat = @stage "galactic" galactic(ra, dec); tick("galactic", t)
    t = time(); qa, qrep = @stage "qa_bitmask" qa_bitmask(id, cfg.lists); tick("qa", t)

    cols = Dict{String,AbstractVector}(
        "tele" => id.tele, "mjd" => id.mjd, "expnum" => id.expnum, "adjfiberindx" => id.adjfiberindx,
        "fiber_index_almanac" => fia, "sdss_id_workup" => id.sdss_id,
        "join_status" => C.join_status, "era" => C.era, "sdss_id_agrees" => C.sdss_id_agrees,
        "mjd_mid_exposure" => mjd_mid, "v_barycentric_correction" => v_bc,
        "v_rel" => v_rel, "v_rad" => v_rad, "e_v_rad" => e_v,
        "bary_flag" => bflag, "glon" => glon, "glat" => glat, "qa_bitmask" => qa)
    for (k, v) in C.i64; cols[k] = v; end
    for (k, v) in C.f64; cols[k] = v; end
    for (k, v) in C.i32; cols[k] = v; end
    for (k, v) in C.f32; cols[k] = v; end
    for (k, v) in C.cats; cols[k] = v; end
    for (k, v) in wk; cols[k] = v; end
    with_off = !isempty(cfg.fiber_offsets)
    if with_off
        t = time()
        vc, ov, fl = @stage "fiber_offsets" fiber_offsets(cfg.fiber_offsets, id, C.era, qa, wk["ingestBit"], wk["RV_flag"], bflag, v_rad)
        cols["v_rad_corr"] = vc; cols["v_rad_offset_applied"] = ov; cols["v_rad_corr_flag"] = fl
        tick("fiber_offsets", t)
    end
    cfg.offtarget && (cols["offtarget"] = C.offtarget)

    rattrs = if cfg.numerics isa Legacy040Numerics
        # 0.4.0 regression: exactly the 0.4.0 root attributes (the IERS record goes to the log and report)
        root_attrs(isempty(B.astropy_version) ? "8.0.1" : B.astropy_version; with_fiber_offsets = with_off)
    else
        root_attrs_041(; with_fiber_offsets = with_off, bary = cfg.bary, astropy_version = B.astropy_version,
                       erfa = erfa_version(), iers = iers_rep)
    end
    t = time(); @stage "write" write_table(cfg.out, cols, C.codes, rattrs, N); tick("write", t)
    timing["total"] = time() - T0

    rep = Dict(
        "n_rows" => N,
        "join_status" => Dict(string(k) => count(==(UInt8(k)), C.join_status) for k in 0:5),
        "era" => Dict(string(k) => count(==(UInt8(k)), C.era) for k in 0:2),
        "sdss_id_agrees" => Dict(string(k) => count(==(UInt8(k)), C.sdss_id_agrees) for k in 0:3),
        "n_finite_vbary" => count(isfinite, v_bc), "n_finite_mid" => count(isfinite, mjd_mid),
        "bary" => B.report, "bary_source" => string(cfg.bary), "erfa_version" => erfa_version(),
        "iers" => iers_rep, "bary_flag_bad_time_exposures" => nbad,
        "bary_flag" => Dict(n => count(x -> x & b != 0, bflag) for (n, b) in
                            (("NO_BARY", NO_BARY), ("BAD_TIME", BAD_TIME), ("ID_MISMATCH", ID_MISMATCH), ("NO_RV", NO_RV))),
        "qa_bits" => Dict(n => count(x -> x & b != 0, qa) for (n, b) in
                          (("BAD_NIGHT", BAD_NIGHT), ("CHRONIC_FIBER", CHRONIC_FIBER), ("BAD_FIBER_NIGHT", BAD_FIBER_NIGHT),
                           ("BAD_FIBER_RANGE", BAD_FIBER_RANGE), ("RV_RAILED", RV_RAILED))),
        "qa_lists" => Dict(string(k) => v for (k, v) in pairs(qrep)),
        "offtarget" => count(==(Int8(1)), C.offtarget),
        "numerics" => string(typeof(cfg.numerics)),
        "threads" => Threads.nthreads(), "workers" => nworkers(),
        "timing_s" => timing)
    isempty(cfg.report) || open(io -> JSON.print(io, rep, 2), cfg.report, "w")
    say(@sprintf("TOTAL %.1f s", timing["total"]))
    return rep
end

"Write all present schema columns, in schema order, with attributes; temp file then rename."
function write_table(out::AbstractString, cols, codes, rattrs, N::Int)
    tmp = out * ".tmp-$(getpid())"
    h5open(tmp, "w") do f
        for spec in SCHEMA
            haskey(cols, spec.name) || continue
            v = cols[spec.name]
            eltype(v) == spec.T || error("column $(spec.name): $(eltype(v)) != schema $(spec.T)")
            length(v) == N || error("column $(spec.name): length $(length(v)) != $N")
            d = write_column(f, spec.name, v, spec.layout)
            spec.name in CATEGORICAL && put_attr!(d, "categories", VlenStrVec(codes[spec.name].order))
            for (k, a) in spec.attrs; put_attr!(d, k, a); end
            close(d)
        end
        for (k, a) in rattrs
            put_attr!(f, k, a === nothing ? N : a)
        end
    end
    mv(tmp, out; force = true)
    return out
end
