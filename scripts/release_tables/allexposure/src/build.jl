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
    python::String = "python3"
    iers_a::String = ""           # pin an IERS-A table for the barycentric kernel ("" = astropy IERS_Auto)
    kernels::String = joinpath(@__DIR__, "..", "python", "legacy_kernels.py")
    numerics::Numerics = LegacyNumerics()
end

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

    # Python kernels (legacy numerics): started now, run while the almanac join reads in parallel
    basis_path = joinpath(cfg.cache, "bary_basis.h5")
    pow_path = joinpath(cfg.cache, "pow10.h5")
    t = time()
    ktask = @async begin
        _run_kernel(cfg, vcat(["bary-basis", "--times", tpath], isempty(cfg.iers_a) ? String[] : ["--iers-a", cfg.iers_a]), basis_path)
        cfg.numerics isa LegacyNumerics && _run_kernel(cfg, ["pow10", "--workup", cfg.workup], pow_path)
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
    B = load_basis(basis_path)
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

    t = time(); @stage "write" write_table(cfg.out, cols, C.codes, root_attrs(B.astropy_version; with_fiber_offsets = with_off), N); tick("write", t)
    timing["total"] = time() - T0

    rep = Dict(
        "n_rows" => N,
        "join_status" => Dict(string(k) => count(==(UInt8(k)), C.join_status) for k in 0:5),
        "era" => Dict(string(k) => count(==(UInt8(k)), C.era) for k in 0:2),
        "sdss_id_agrees" => Dict(string(k) => count(==(UInt8(k)), C.sdss_id_agrees) for k in 0:3),
        "n_finite_vbary" => count(isfinite, v_bc), "n_finite_mid" => count(isfinite, mjd_mid),
        "bary" => B.report, "bary_flag_bad_time_exposures" => nbad,
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
