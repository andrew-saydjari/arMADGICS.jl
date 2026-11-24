# derived.jl -- the row-wise stages: barycentric correction, velocities, bary_flag, Galactic
# coordinates, qa_bitmask, fiber zero points. All threaded over rows (no HDF5 inside the loops).
#
# NUMERICS. Bit-identity with 0.4.0 pins three implementation choices (README, "Numerics"):
#   * cos/sin of the per-row target unit vector: numpy float64 cos/sin match glibc's libm bit for
#     bit, Julia's Base.cos/sin do not always; so `LegacyNumerics` calls libm through ccall.
#   * 10^(p*6e-6): numpy uses Intel SVML on AVX-512 hardware (1-ulp differences from correctly
#     rounded pow in ~7% of rows); `LegacyNumerics` takes that array from the Python kernel.
#   * Galactic: astropy goes ICRS -> FK5(J2000) -> Galactic with ERFA s2c / rxp / c2s and its own
#     unit scalings and longitude wrapping; reproduced here with ERFA's C library (ERFA_jll).
# `NativeNumerics` uses Base.cos/sin and Base.exp10 (pure Julia, CPU-independent); it is NOT
# bit-identical to 0.4.0 (differences are reported by the regression in ulps).

abstract type Numerics end
struct LegacyNumerics <: Numerics end
struct NativeNumerics <: Numerics end

@inline _cos(::LegacyNumerics, x::Float64) = ccall(:cos, Float64, (Float64,), x)
@inline _sin(::LegacyNumerics, x::Float64) = ccall(:sin, Float64, (Float64,), x)
@inline _cos(::NativeNumerics, x::Float64) = cos(x)
@inline _sin(::NativeNumerics, x::Float64) = sin(x)

const C_KM_S = 299792.458
const DELTA_LOG = 6e-6
const DEG2RAD_NP = 0.017453292519943295      # numpy.radians factor (== pi/180 in double)

# ------------------------------------------------------------------------------ barycentric

"Per-exposure basis coefficients from the Python kernel (legacy_kernels.py bary-basis)."
struct BaryBasis
    exp_key::Dict{Int,Vector{Int64}}       # tele => exposure keys (exposure_times.h5 order)
    U::Dict{Int,Matrix{Float64}}           # tele => 4 x ne
    astropy_version::String
    report::Dict{String,Any}
end
function load_basis(path)
    h5open(path, "r") do f
        ek = Dict{Int,Vector{Int64}}(); U = Dict{Int,Matrix{Float64}}(); rep = Dict{String,Any}()
        for (t, tn) in TELENAME
            g = f[tn]
            ek[t] = read(g["exp_key"]); U[t] = read(g["U"])
            # HDF5.jl reads the (4, ne) numpy array as (ne, 4): transpose to 4 x ne
            size(U[t], 2) == 4 && size(U[t], 1) == length(ek[t]) && (U[t] = permutedims(U[t]))
            haskey(attributes(g), "n_exposures") &&
                (rep[tn] = Dict("n_exposures" => read(attributes(g)["n_exposures"]),
                                "basis_residual_m_s" => read(attributes(g)["basis_residual_m_s"])))
        end
        for k in ("iers_table", "iers_table_sha256")
            haskey(attributes(f), k) && (rep[k] = read(attributes(f)[k]))
        end
        BaryBasis(ek, U, String(read(attributes(f)["astropy_version"])), rep)
    end
end

"v_bc per row = c*(u0 + ux nx + uy ny + uz nz - 1), same operation order as 0.4.0 s6."
function barycentric(nm::Numerics, B::BaryBasis, id::Identity, ra::Vector{Float64}, dec::Vector{Float64})
    N = length(id)
    v = fill(NaN, N)
    for t in (1, 2)
        ek = B.exp_key[t]; U = B.U[t]; ne = length(ek)
        ne == 0 && continue
        so = sortperm(ek); eks = ek[so]
        Threads.@threads for i in 1:N
            id.tele[i] == t || continue
            (isfinite(ra[i]) && isfinite(dec[i])) || continue
            rk = expkey(id.tele[i], id.mjd[i], id.expnum[i])
            p = clamp(searchsortedfirst(eks, rk), 1, ne)
            eks[p] == rk || continue
            e = so[p]
            rr = ra[i] * DEG2RAD_NP; dd = dec[i] * DEG2RAD_NP
            cd = _cos(nm, dd)
            nx = cd * _cos(nm, rr); ny = cd * _sin(nm, rr); nz = _sin(nm, dd)
            v[i] = C_KM_S * (U[1, e] + U[2, e] * nx + U[3, e] * ny + U[4, e] * nz - 1.0)
        end
    end
    return v
end

# ------------------------------------------------------------------------------ velocities

"10^(p*6e-6): legacy = numpy/SVML array from the Python kernel; native = Base.exp10."
pow10(::LegacyNumerics, p::Vector{Float64}, legacy_path::AbstractString) =
    h5open(f -> read(f["pow10"]), legacy_path, "r")
pow10(::NativeNumerics, p::Vector{Float64}, ::AbstractString) = exp10.(p .* DELTA_LOG)

"v_rel, v_rad, e_v_rad (almanac finalize_radial_velocities convention, 0.4.0 s6 op order)."
function velocities(pw::Vector{Float64}, var::Vector{Float64}, v_bc::Vector{Float64})
    N = length(pw)
    v_rel = Vector{Float64}(undef, N); v_rad = similar(v_rel); e_v = similar(v_rel)
    k = DELTA_LOG * log(10.0)
    Threads.@threads for i in 1:N
        z = pw[i] - 1
        zc = v_bc[i] / C_KM_S
        zp = z + 1
        v_rel[i] = ((zp^2 - 1) / (zp^2 + 1)) * C_KM_S
        zt = z + zc + zc * z
        ztp = zt + 1
        v_rad[i] = ((ztp^2 - 1) / (ztp^2 + 1)) * C_KM_S
        e_z = k * pw[i] * Base.sqrt_llvm(var[i])     # IEEE sqrt: NaN (not an error) for var < 0, as numpy
        e_v[i] = abs(4 * zp / ((zp^2 + 1)^2)) * C_KM_S * e_z
    end
    return v_rel, v_rad, e_v
end

# ------------------------------------------------------------------------------ bary_flag

const NO_BARY, BAD_TIME, ID_MISMATCH, NO_RV = 0x01, 0x02, 0x04, 0x08

"bary_flag (0.4.0 s9): NO_BARY | BAD_TIME | ID_MISMATCH | NO_RV."
function bary_flag(id::Identity, et::ExposureTimes, v_bc, agrees, pixoff)
    bad = Set{Int64}()
    for j in eachindex(et.tele)
        off = et.mjd_mid_precise[j] - et.mjd[j]
        (off < -0.6 || off > 1.6 || !isfinite(et.mjd_mid_precise[j])) &&
            push!(bad, expkey(et.tele[j], et.mjd[j], et.expnum[j]))
    end
    N = length(id); fl = zeros(UInt8, N)
    Threads.@threads for i in 1:N
        f = 0x00
        isfinite(v_bc[i]) || (f |= NO_BARY)
        expkey(id.tele[i], id.mjd[i], id.expnum[i]) in bad && (f |= BAD_TIME)
        agrees[i] == 0x02 && (f |= ID_MISMATCH)
        isfinite(pixoff[i]) || (f |= NO_RV)
        fl[i] = f
    end
    return fl, length(bad)
end

# ------------------------------------------------------------------------------ Galactic

const LIBERFA = ERFA_jll.liberfa
# astropy 8.0.1 ICRS->FK5(J2000) and FK5(J2000)->Galactic matrices (row-major), as hex floats so
# the bits are exact. Extracted from astropy.coordinates.builtin_frames (icrs_to_fk5, fk5_to_gal).
const M_ICRS_FK5 = (0x1.fffffffffffb9p-1, -0x1.dcd6589697144p-24, -0x1.7af89b5df8a29p-25,
                    0x1.dcd659c94c1dap-24, 0x1.fffffffffff9cp-1, 0x1.9e5e98fa40aeep-24,
                    0x1.7af8955652f0dp-25, -0x1.9e5e9a66d39d8p-24, 0x1.fffffffffffcdp-1)
const M_FK5_GAL = (-0x1.c18a58e1c22ddp-5, -0x1.bf3325efd5a22p-1, -0x1.ef7270a9417fcp-2,
                   0x1.f9f7d41fd3bccp-2, -0x1.c78168bf39b7ep-2, 0x1.7e747511bc55fp-1,
                   -0x1.bc3ebc593c29fp-1, -0x1.95a912fbe9079p-3, 0x1.d2ed6a73d5f56p-2)
const AP_DEG2RAD = 0x1.1df46a2529d39p-6      # astropy u.deg.to(u.rad)
const AP_RAD2DEG = 0x1.ca5dc1a63c1f8p+5      # astropy u.rad.to(u.deg)
const AP_A360    = 0x1.921fb54442d18p+2      # astropy u.degree.to(u.rad, 360.0)

"astropy Angle._wrap_at(360 deg) for one value in radians (wrap_angle_floor = 0)."
@inline function _wrap360(x::Float64)
    (x < 0.0 || x >= AP_A360) || return x
    w = fld(x - 0.0, AP_A360)
    x -= w * AP_A360
    x >= AP_A360 && (x -= AP_A360)
    x < 0.0 && (x += AP_A360)
    return x
end

@inline function _rotate(M::NTuple{9,Float64}, lon::Float64, lat::Float64)
    p = Ref{NTuple{3,Float64}}(); q = Ref{NTuple{3,Float64}}(); m = Ref(M)
    th = Ref(0.0); ph = Ref(0.0)
    ccall((:eraS2c, LIBERFA), Cvoid, (Cdouble, Cdouble, Ptr{NTuple{3,Float64}}), lon, lat, p)
    ccall((:eraRxp, LIBERFA), Cvoid, (Ptr{NTuple{9,Float64}}, Ptr{NTuple{3,Float64}}, Ptr{NTuple{3,Float64}}), m, p, q)
    ccall((:eraC2s, LIBERFA), Cvoid, (Ptr{NTuple{3,Float64}}, Ref{Cdouble}, Ref{Cdouble}), q, th, ph)
    return _wrap360(th[]), ph[]
end

"(glon, glat) in degrees, NaN where ra/dec not finite (0.4.0 s10 via astropy)."
function galactic(ra::Vector{Float64}, dec::Vector{Float64})
    N = length(ra); gl = fill(NaN, N); gb = fill(NaN, N)
    Threads.@threads for i in 1:N
        (isfinite(ra[i]) && isfinite(dec[i])) || continue
        lon, lat = _rotate(M_ICRS_FK5, ra[i] * AP_DEG2RAD, dec[i] * AP_DEG2RAD)
        lon, lat = _rotate(M_FK5_GAL, lon, lat)
        gl[i] = lon * AP_RAD2DEG; gb[i] = lat * AP_RAD2DEG
    end
    return gl, gb
end

# ------------------------------------------------------------------------------ qa_bitmask

const BAD_NIGHT, CHRONIC_FIBER, BAD_FIBER_NIGHT, BAD_FIBER_RANGE, RV_RAILED = 1, 2, 4, 8, 16
_tcode(s) = s == "apo" ? 1 : s == "lco" ? 2 : error("unknown tele $s")
_int(x) = x isa Integer ? Int(x) : parse(Int, string(x))

"qa_bitmask from bad_nights.csv, bad_fibers.csv, rv_railed.csv (0.4.0 s10). fiber = adjfiberindx - 300(tele-1)."
function qa_bitmask(id::Identity, lists::AbstractString)
    bn = CSV.File(joinpath(lists, "bad_nights.csv"); types = Dict(:tele => String))
    bf = CSV.File(joinpath(lists, "bad_fibers.csv"); types = Dict(:tele => String, :scope => String))
    rv = CSV.File(joinpath(lists, "rv_railed.csv"); types = Dict(:tele => String))
    nights = Set((_tcode(r.tele), _int(r.mjd)) for r in bn)
    chronic = Set((_tcode(r.tele), _int(r.fiber)) for r in bf if r.scope == "chronic")
    fnight = Set((_tcode(r.tele), _int(r.mjd), _int(r.fiber)) for r in bf if r.scope == "fiber_night")
    franges = [(_tcode(r.tele), _int(r.fiber), _int(r.mjd_start), _int(r.mjd_end)) for r in bf if r.scope == "fiber_mjd_range"]
    other = [r.scope for r in bf if !(r.scope in ("chronic", "fiber_night", "fiber_mjd_range"))]
    isempty(other) || error("unexpected scope values: $(unique(other))")
    railed = Set((_tcode(r.tele), _int(r.mjd), _int(r.expnum), _int(r.adjfiberindx)) for r in rv)
    N = length(id); qa = zeros(Int32, N)
    Threads.@threads for i in 1:N
        t = Int(id.tele[i]); m = Int(id.mjd[i]); fb = Int(id.adjfiberindx[i]) - 300 * (t - 1)
        q = Int32(0)
        (t, m) in nights && (q |= BAD_NIGHT)
        (t, fb) in chronic && (q |= CHRONIC_FIBER)
        (t, m, fb) in fnight && (q |= BAD_FIBER_NIGHT)
        for (tt, ff, a, b) in franges
            (tt == t && ff == fb && a <= m <= b) && (q |= BAD_FIBER_RANGE; break)
        end
        (t, m, Int(id.expnum[i]), Int(id.adjfiberindx[i])) in railed && (q |= RV_RAILED)
        qa[i] = q
    end
    return qa, (n_bad_nights = length(nights), n_chronic = length(chronic), n_fiber_night = length(fnight),
                n_fiber_range = length(franges), n_railed = length(railed))
end

# ------------------------------------------------------------------------------ fiber offsets

"Release default mask (RECOMMENDED_CUTS.md, without offtarget): the mask v_rad_corr used."
@inline default_mask(qa, ing, rvf, bf, v) = qa == 0 && ing == 0 && rvf == 0 && bf == 0 && isfinite(v)

"v_rad_corr, v_rad_offset_applied, v_rad_corr_flag (0.4.0 fiber_offsets/s10_final_build.jl)."
function fiber_offsets(csv::AbstractString, id::Identity, era, qa, ing, rvf, bf, v)
    off = fill(NaN, 1200); cal = falses(1200)
    d = CSV.File(csv; comment = "#")
    for r in d
        b = r.fiber + (r.tele == "lco" ? 300 : 0) + 600 * (r.era - 1)
        off[b] = r.offset_kms; cal[b] = r.calibratable
    end
    N = length(v); vc = fill(NaN, N); ov = fill(NaN, N); fl = zeros(UInt8, N)
    Threads.@threads for i in 1:N
        b = id.adjfiberindx[i] + 600 * (Int(era[i]) - 1)
        if !default_mask(qa[i], ing[i], rvf[i], bf[i], v[i])
            fl[i] = 0x02
        elseif !cal[b]
            fl[i] = 0x01
        else
            ov[i] = off[b]; vc[i] = v[i] - ov[i]
        end
    end
    return vc, ov, fl
end
