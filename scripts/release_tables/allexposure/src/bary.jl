# bary.jl -- the barycentric correction in Julia: a port of astropy 8.0.1
#   SkyCoord.radial_velocity_correction(kind="barycentric", obstime=Time(mjd, scale="utc", location=site))
# with astropy's builtin ephemeris, evaluated with ERFA's C library (ERFA_jll, the library pyerfa
# wraps) and an explicit IERS-A table (finals2000A format). Nothing is downloaded at run time.
#
# The correction is linear in the target unit vector n (no distance, no proper motion):
#     v_bc(n) = c * (zb(n) - 1),   zb(n) = gamma * (1 + beta . n) / (1 + gr / c)
# so per exposure it is fully described by U = (u0, ux, uy, uz) = (gamma, gamma*beta) / (1 + gr/c),
# the same 4 numbers the 0.4.0 build fitted from 6 astropy evaluations (legacy_kernels.py bary-basis).
# Here U is computed directly. Step by step, as astropy does it:
#   times      Time(mjd, utc): astropy day_frac; UTC->TAI->TT (eraUtctai, eraTaitt); TT->TDB with
#              eraDtdb at the site (astropy Time._get_delta_tdb_tt); UTC->UT1 with UT1-UTC linearly
#              interpolated in the IERS table (astropy IERS._interpolate).
#   observer   EarthLocation.get_gcrs_posvel: CIRS->ITRS = eraC2tcio(I, eraEra00(UT1),
#              eraPom00(xp, yp, eraSp00(TT))), GCRS->CIRS = eraC2i06a(TT); position = R * site,
#              velocity = (OMEGA_EARTH * CIP axis in GCRS) x position.
#   Earth      eraEpv00(TDB) barycentric velocity.
#   gr         EarthLocation.gravitational_redshift: -GM/(c r) summed over Earth (geocentric site
#              distance), Moon (eraMoon98), Jupiter (eraPlan94) and Sun (eraEpv00), positions at TDB.
# The unit handling (AU/d vs m/s, astropy's float conversion factors) follows astropy's Quantity
# arithmetic; the conversion factors are taken from astropy 8.0.1 as exact hex floats.

const MJD_ZERO = 2400000.5

# ------------------------------------------------------------------------------ IERS table

"""
The Earth-orientation table astropy 8.0.1 builds with IERS_Auto (what 0.4.0 used): an IERS-A file
(finals2000A format; rapid service + predictions, with Bulletin B columns) whose Bulletin B values are
overwritten by an IERS-B file (EOP 20 C04, `eopc04.1962-now`, IERS_Auto._substitute_iers_b), then
combined with B values where present, else A (IERS_A._combine_a_b_columns).
"""
struct EOPTable
    path::String                 # IERS-A file
    sha256::String
    nbytes::Int
    path_b::String               # IERS-B file ("" = none: plain astropy IERS_A)
    sha256_b::String
    mjd::Vector{Float64}
    ut1_utc::Vector{Float64}     # s
    pm_x::Vector{Float64}        # arcsec
    pm_y::Vector{Float64}        # arcsec
    first_mjd::Float64
    last_mjd::Float64
    predictive_mjd::Float64      # first row whose Bulletin A UT1 or PM value is a prediction ('P')
    last_bulletin_b_mjd::Float64 # last row with Bulletin B (final) UT1 and PM values
    last_iers_b_mjd::Float64     # last row whose values came from the IERS-B file
end

_fld(l, a, b) = (s = strip(l[min(a, end+1):min(b, end)]); isempty(s) ? NaN : parse(Float64, s))
_flag(l, a) = length(l) >= a ? l[a] : ' '

"Read an eopc04 (IERS-B, EOP 20 C04) file: MJD, PM_x, PM_y (arcsec), UT1_UTC (s); columns per its ReadMe."
function read_eopc04(path::AbstractString)
    mjd = Float64[]; px = Float64[]; py = Float64[]; ut = Float64[]
    for l in eachline(path)
        (isempty(strip(l)) || startswith(l, "#")) && continue
        push!(mjd, _fld(l, 17, 26)); push!(px, _fld(l, 27, 38)); push!(py, _fld(l, 39, 50)); push!(ut, _fld(l, 51, 62))
    end
    return mjd, px, py, ut
end

"""
    read_finals2000A(path_a; iers_b = "")

Read a finals2000A.all file (columns per its ReadMe; astropy IERS_A.read). With `iers_b` (an eopc04
file) its values replace the Bulletin B columns first, as astropy IERS_Auto does.
"""
function read_finals2000A(path::AbstractString; iers_b::AbstractString = "")
    bytes = read(path)
    mjd = Float64[]; ua = Float64[]; pxa = Float64[]; pya = Float64[]; ub = Float64[]; pxb = Float64[]; pyb = Float64[]
    pmflag = Char[]; utflag = Char[]
    for l in eachline(IOBuffer(bytes))
        isempty(strip(l)) && continue
        u = _fld(l, 59, 68); pf = _flag(l, 17)
        (isfinite(u) && pf != ' ') || continue          # astropy drops rows without a UT1/PM value
        push!(mjd, _fld(l, 8, 15)); push!(ua, u); push!(pxa, _fld(l, 19, 27)); push!(pya, _fld(l, 38, 46))
        push!(pxb, _fld(l, 135, 144)); push!(pyb, _fld(l, 145, 154)); push!(ub, _fld(l, 155, 165))
        push!(pmflag, pf); push!(utflag, _flag(l, 58))
    end
    issorted(mjd) || error("$path: MJD column not sorted")
    sha_b = ""; nb = 0
    if !isempty(iers_b)
        bb = read(iers_b); sha_b = bytes2hex(SHA.sha256(bb))
        bm, bx, by, bu = read_eopc04(iers_b)
        hasb = findall(isfinite, ub)
        i0 = searchsortedfirst(bm, mjd[hasb[1]]); i1 = searchsortedlast(bm, mjd[hasb[end]])
        nb = max(0, i1 - i0 + 1)
        if nb > 0
            all(isapprox.(mjd[1:nb], bm[i0:i1])) || error("IERS-B $iers_b does not line up with $path")
            ub[1:nb] .= bu[i0:i1]; pxb[1:nb] .= bx[i0:i1]; pyb[1:nb] .= by[i0:i1]
        end
    end
    ut1 = [isnan(ub[i]) ? ua[i] : ub[i] for i in eachindex(mjd)]
    bpm = [isnan(pxb[i]) || isnan(pyb[i]) for i in eachindex(mjd)]
    px = [bpm[i] ? pxa[i] : pxb[i] for i in eachindex(mjd)]
    py = [bpm[i] ? pya[i] : pyb[i] for i in eachindex(mjd)]
    p_ut = something(findfirst(==('P'), utflag), length(mjd))
    p_pm = something(findfirst(==('P'), pmflag), length(mjd))
    p = min(p_ut, p_pm)                                   # astropy meta["predictive_mjd"]
    lb = findlast(i -> !isnan(ub[i]) && !bpm[i], eachindex(mjd))
    EOPTable(String(path), bytes2hex(SHA.sha256(bytes)), length(bytes), String(iers_b), sha_b, mjd, ut1, px, py,
             first(mjd), last(mjd), mjd[p], lb === nothing ? NaN : mjd[lb], nb > 0 ? mjd[nb] : NaN)
end

"astropy IERS._interpolate at a two-part UTC JD; returns (value(s), status) with status 0 = in table."
@inline function _eop_index(T::EOPTable, jd1::Float64, jd2::Float64)
    mjd = floor(jd1 - MJD_ZERO + jd2)
    utc = jd1 - (MJD_ZERO + mjd) + jd2
    n = length(T.mjd)
    i = searchsortedlast(T.mjd, mjd)                     # == numpy searchsorted(side="right")
    i1 = clamp(i, 1, n - 1)                              # 0-based i1; rows i1 (i0) and i1+1 (i1), 1-based
    return mjd, utc, i, i1
end
function eop_ut1_utc(T::EOPTable, jd1::Float64, jd2::Float64)
    mjd, utc, i, i1 = _eop_index(T, jd1, jd2)
    n = length(T.mjd)
    i == 0 && return T.ut1_utc[1]
    i == n && return T.ut1_utc[n]
    m0 = T.mjd[i1]; m1 = T.mjd[i1+1]; v0 = T.ut1_utc[i1]; d = T.ut1_utc[i1+1] - v0
    d -= round(d)                                        # leap-second jump (numpy round: half to even)
    return v0 + (mjd - m0 + utc) / (m1 - m0) * d
end
function eop_pm(T::EOPTable, jd1::Float64, jd2::Float64)
    mjd, utc, i, i1 = _eop_index(T, jd1, jd2)
    n = length(T.mjd)
    i == 0 && return T.pm_x[1], T.pm_y[1]
    i == n && return T.pm_x[n], T.pm_y[n]
    m0 = T.mjd[i1]; m1 = T.mjd[i1+1]; f = (mjd - m0 + utc) / (m1 - m0)
    return T.pm_x[i1] + f * (T.pm_x[i1+1] - T.pm_x[i1]), T.pm_y[i1] + f * (T.pm_y[i1+1] - T.pm_y[i1])
end

# ------------------------------------------------------------------------------ sites, constants

"An observatory as astropy EarthLocation.of_site gives it (geocentric ITRS, m; dtdb arguments)."
struct Site
    name::String
    xyz::NTuple{3,Float64}   # m
    lon_rad::Float64         # location.lon.to_value(rad)
    rxy_km::Float64          # hypot(x, y) in km
    z_km::Float64
end
# astropy 8.0.1 EarthLocation.of_site (sites registry; WGS84):
#   apo: lon -105.82 deg, lat 32.78 deg, height 2798 m;  lco: lon -70.70167, lat -29.00333, height 2282 m
const SITE_APO = Site("apo", (-0x1.656a14d462780p+20, -0x1.3b59455e72acep+22, 0x1.a34f4db244fdep+21),
                      -0x1.d8ceec9fd8eafp+0, 0x1.4fa1388d9a40fp+12, 0x1.ad5f8af3fceb6p+11)
const SITE_LCO = Site("lco", (0x1.c29977fc1f6cdp+20, -0x1.41b5212dca735p+22, -0x1.7768163889ecdp+21),
                      -0x1.3be5e87c96960p+0, 0x1.5d0a7a739a2c6p+12, -0x1.806a95baeb718p+11)
const SITES = (SITE_APO, SITE_LCO)                     # index = tele code

const C_M_S = 299792458.0
const ARCSEC2RAD = 0x1.455a5b2ff8f9dp-18                # astropy u.arcsec.to(u.rad)
const AUD_TO_MS  = 0x1.a6b80d638e38ep+20                # (u.au/u.day).to(u.m/u.s)
const MS_TO_AUD  = 0x1.3611a33fcf53bp-21               # (u.m/u.s).to(u.au/u.day)
const AUD_MS2    = 0x1.5d01a83bbc59bp+41                # ((u.au/u.day)/(u.m/u.s))**2 -> dimensionless
const AU_M       = 0x1.16a5d2d360000p+37                # u.au.to(u.m)
const OMEGA_EARTH = 0x1.31da7d7cba179p-14               # astropy.coordinates.earth.OMEGA_EARTH, 1/s
const GM_SUN   = 0x1.cc70534284d5ep+66                  # astropy.constants, m^3/s^2
const GM_JUP   = 0x1.c214ba527a140p+56
const GM_EARTH = 0x1.6a86633fd4000p+48
const GM_MOON  = 0x1.1d74268bf6accp+42                  # G * 7.34767309e22 kg (astropy default)

# ------------------------------------------------------------------------------ ERFA calls

const M3 = NTuple{9,Float64}       # row-major 3x3, as ERFA's double[3][3]
const V3 = NTuple{3,Float64}
const PV = NTuple{6,Float64}       # double[2][3]: p then v

macro erfa(f, args...)
    esc(:(ccall(($(QuoteNode(f)), LIBERFA), $(args...))))
end
function _jd2(f::Symbol, a::Float64, b::Float64)
    x = Ref(0.0); y = Ref(0.0)
    st = if f === :utctai
        @erfa eraUtctai Cint (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}) a b x y
    elseif f === :taitt
        @erfa eraTaitt Cint (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}) a b x y
    elseif f === :tttai
        @erfa eraTttai Cint (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}) a b x y
    elseif f === :taiutc
        @erfa eraTaiutc Cint (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}) a b x y
    else
        error("unknown $f")
    end
    st < 0 && error("ERFA $f failed ($st)")
    return x[], y[]
end
function _tttdb(a, b, dtr)
    x = Ref(0.0); y = Ref(0.0)
    @erfa eraTttdb Cint (Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}) a b dtr x y
    return x[], y[]
end
function _utcut1(a, b, dut1)
    x = Ref(0.0); y = Ref(0.0)
    st = @erfa eraUtcut1 Cint (Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}) a b dut1 x y
    st < 0 && error("ERFA utcut1 failed ($st)")
    return x[], y[]
end
_dtdb(a, b, ut, el, u, v) = @erfa eraDtdb Cdouble (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble) a b ut el u v
_sp00(a, b) = @erfa eraSp00 Cdouble (Cdouble, Cdouble) a b
_era00(a, b) = @erfa eraEra00 Cdouble (Cdouble, Cdouble) a b
function _pom00(xp, yp, sp)
    r = Ref{M3}(); @erfa eraPom00 Cvoid (Cdouble, Cdouble, Cdouble, Ref{M3}) xp yp sp r; r[]
end
function _c2i06a(a, b)
    r = Ref{M3}(); @erfa eraC2i06a Cvoid (Cdouble, Cdouble, Ref{M3}) a b r; r[]
end
function _c2tcio(c2i::M3, era, pom::M3)
    r = Ref{M3}(); @erfa eraC2tcio Cvoid (Ref{M3}, Cdouble, Ref{M3}, Ref{M3}) Ref(c2i) era Ref(pom) r; r[]
end
function _epv00(a, b)
    h = Ref{PV}(); bb = Ref{PV}(); @erfa eraEpv00 Cint (Cdouble, Cdouble, Ref{PV}, Ref{PV}) a b h bb; h[], bb[]
end
function _moon98(a, b)
    r = Ref{PV}(); @erfa eraMoon98 Cvoid (Cdouble, Cdouble, Ref{PV}) a b r; r[]
end
function _plan94(a, b, np)
    r = Ref{PV}(); st = @erfa eraPlan94 Cint (Cdouble, Cdouble, Cint, Ref{PV}) a b np r
    st < 0 && error("ERFA plan94 failed ($st)"); r[]
end
_pm(p::V3) = @erfa eraPm Cdouble (Ref{V3},) Ref(p)
_pdp(a::V3, b::V3) = @erfa eraPdp Cdouble (Ref{V3}, Ref{V3}) Ref(a) Ref(b)
function _pxp(a::V3, b::V3)
    r = Ref{V3}(); @erfa eraPxp Cvoid (Ref{V3}, Ref{V3}, Ref{V3}) Ref(a) Ref(b) r; r[]
end
function _rxp(m::M3, p::V3)
    r = Ref{V3}(); @erfa eraRxp Cvoid (Ref{M3}, Ref{V3}, Ref{V3}) Ref(m) Ref(p) r; r[]
end

@inline _T(m::M3) = (m[1], m[4], m[7], m[2], m[5], m[8], m[3], m[6], m[9])
"""
3x3 matrix product (row-major) as astropy's numpy `@` computes it through OpenBLAS dgemm: per
element, products accumulated in k order with fused multiply-adds. Julia's `fma` is correctly rounded
on every CPU (hardware or software), so this is deterministic; it reproduces astropy's matrices on
the AVX-512 head node for all 144,381 exposures (plain mul+add differs by a few ulp in ~25 %).
"""
@inline function _mm(a::M3, b::M3)
    ntuple(9) do q
        i = (q - 1) ÷ 3; j = (q - 1) % 3
        fma(a[3i+3], b[j+7], fma(a[3i+2], b[j+4], a[3i+1] * b[j+1]))
    end
end
@inline _p(pv::PV) = (pv[1], pv[2], pv[3])
@inline _v(pv::PV) = (pv[4], pv[5], pv[6])
@inline _sub(a::V3, b::V3) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _add(a::V3, b::V3) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])

# ------------------------------------------------------------------------------ astropy Time

@inline function _two_sum(a::Float64, b::Float64)
    x = a + b; eb = x - a; ea = x - eb; eb = b - eb; ea = a - ea
    return x, ea + eb
end
"astropy.time.utils.day_frac(val1, val2) (no factor/divisor)."
@inline function day_frac(v1::Float64, v2::Float64)
    s, e = _two_sum(v1, v2)
    day = round(s)
    frac, check = _two_sum(s - day, e)
    excess = frac * sign(check) != 0.5 ? round(frac) : round(frac + 2 * check)
    day += excess
    frac = s - day
    frac += e
    return day, frac
end
"Time(mjd, format='mjd', scale='utc') -> (jd1, jd2)  (astropy TimeMJD.set_jds)."
@inline function mjd_to_jd(mjd::Float64)
    j1, j2 = day_frac(mjd, 0.0)
    j1 += MJD_ZERO
    return day_frac(j1, j2)
end

"UTC two-part JD -> (TT, TDB, UT1) two-part JDs and UT1-UTC, exactly as astropy Time does it."
function time_scales(site::Site, T::EOPTable, u1::Float64, u2::Float64)
    a1, a2 = _jd2(:utctai, u1, u2)
    t1, t2 = _jd2(:taitt, a1, a2)
    tt = day_frac(t1, t2)
    # TT -> TDB: delta from eraDtdb at the (un-normalised) TT JD, astropy _get_delta_tdb_tt
    n1, n2 = _jd2(:tttai, t1, t2)
    n1, n2 = _jd2(:taiutc, n1, n2)
    ut = day_frac(n1 - 0.5, n2)[2]
    dtr = _dtdb(t1, t2, ut, site.lon_rad, site.rxy_km, site.z_km)
    tdb = day_frac(_tttdb(t1, t2, dtr)...)
    dut1 = eop_ut1_utc(T, u1, u2)
    ut1 = day_frac(_utcut1(u1, u2, dut1)...)
    return tt, tdb, ut1, dut1
end

# ------------------------------------------------------------------------------ the correction

# The 0.4.0 build evaluated astropy for 6 directions per exposure and fitted U = pinv(A) @ Y
# (s6_build_allexposure.py; legacy_kernels.py bary-basis). bary_coeffs does the same in Julia, with
# astropy's direction vectors and numpy's pinv(A) as exact hex floats, and the 4x6 product
# accumulated as OpenBLAS dgemm does (fused multiply-adds in k order). Result: the U of the astropy
# kernel bit for bit (README, "Barycentric correction").
const BASIS_DIRS = (   # SkyCoord(ra, dec).icrs.cartesian for (ra, dec) = (0,0) (180,0) (90,0) (270,0) (0,90) (0,-90) deg
    (0x1.0000000000000p+0, 0.0, 0.0),
    (-0x1.0000000000000p+0, 0x1.1a62633145c07p-53, 0.0),
    (0x1.1a62633145c07p-54, 0x1.0000000000000p+0, 0.0),
    (-0x1.a79394c9e8a0ap-53, -0x1.0000000000000p+0, 0.0),
    (0x1.1a62633145c07p-54, 0.0, 0x1.0000000000000p+0),
    (0x1.1a62633145c07p-54, 0.0, -0x1.0000000000000p+0))
const APINV = (        # numpy.linalg.pinv([ones(6) BAS]), 4 x 6
    (0x1.5555555555558p-3, 0x1.5555555555557p-3, 0x1.5555555555557p-3, 0x1.5555555555557p-3, 0x1.5555555555557p-3, 0x1.5555555555557p-3),
    (0x1.ffffffffffffep-2, -0x1.ffffffffffffcp-2, 0x1.644f8dd8c2de4p-55, 0x1.644f8dd8c2de4p-55, 0x1.644f8dd8c2de4p-55, 0x1.644f8dd8c2de4p-55),
    (0x1.7a97d1c75937ep-54, 0x1.8891e06cca878p-55, 0x1.ffffffffffffcp-2, -0x1.ffffffffffffep-2, 0x1.dc1ecda59b30ap-56, 0x1.dc1ecda59b30ap-56),
    (-0x1.85fef5c60c875p-58, -0x1.c43c7056a0e1ep-60, -0x1.c43c7056a0e1ep-60, -0x1.c43c7056a0e1ep-60, 0x1.ffffffffffffep-2, -0x1.ffffffffffffep-2))

"""
    bary_coeffs(site, eop, mjd_utc) -> (u0, ux, uy, uz)

Per-exposure coefficients of v_bc(n) = c * (u0 + ux nx + uy ny + uz nz - 1) (c in km/s), by the
0.4.0 method: astropy's radial_velocity_correction for the 6 basis directions, then pinv.
"""
function bary_coeffs(site::Site, T::EOPTable, mjd::Float64)
    st = bary_state(site, T, mjd)
    Y = ntuple(d -> _vkms(st, BASIS_DIRS[d]) / C_KM_S + 1.0, 6)
    return ntuple(4) do r
        a = APINV[r]
        acc = a[1] * Y[1]
        for k in 2:6
            acc = fma(a[k], Y[k], acc)
        end
        acc
    end
end

"""
    bary_coeffs_direct(site, eop, mjd_utc) -> (u0, ux, uy, uz)

The same coefficients computed analytically, U = (gamma, gamma*beta) / (1 + gr/c), without the
6-direction fit (whose rounding is up to ~2.7e-10 km/s). Not used by the build; for comparisons.
"""
function bary_coeffs_direct(site::Site, T::EOPTable, mjd::Float64)
    st = bary_state(site, T, mjd)
    u0 = st.gamma / st.den
    k = st.gamma * AUD_TO_MS / st.den
    return (u0, k * st.beta[1], k * st.beta[2], k * st.beta[3])
end

"""
    bary_state(site, eop, mjd_utc)

The intermediate quantities of astropy's radial_velocity_correction at one time: `beta` in astropy's
units (AU/d per m/s; multiply by AUD_TO_MS for dimensionless), `gamma`, `den` = 1 + gr/c, the
gravitational redshift `gr` (m/s), the observer's GCRS position (m) and velocity (m/s).
"""
function bary_state(site::Site, T::EOPTable, mjd::Float64)
    u1, u2 = mjd_to_jd(mjd)
    (tt1, tt2), (tb1, tb2), (ut11, ut12), _ = time_scales(site, T, u1, u2)
    # observer GCRS position and velocity (EarthLocation._get_gcrs_posvel)
    xpa, ypa = eop_pm(T, u1, u2)
    xp = xpa * ARCSEC2RAD; yp = ypa * ARCSEC2RAD
    pom = _pom00(xp, yp, _sp00(tt1, tt2))
    ref_to_itrs = _c2tcio((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), _era00(ut11, ut12), pom)
    gcrs_to_ref = _c2i06a(tt1, tt2)
    ref_to_gcrs = _T(gcrs_to_ref)
    itrs_to_gcrs = _mm(ref_to_gcrs, _T(ref_to_itrs))
    rot = (ref_to_gcrs[3] * OMEGA_EARTH, ref_to_gcrs[6] * OMEGA_EARTH, ref_to_gcrs[9] * OMEGA_EARTH)
    gpos = _rxp(itrs_to_gcrs, site.xyz)                  # m
    gvel = _pxp(rot, gpos)                               # m/s
    # Earth barycentric velocity (AU/d), plus bodies for the gravitational redshift (AU)
    eh, eb = _epv00(tb1, tb2)
    ve = _v(eb)
    pe = _p(eb)
    psun = _p(eb) .- _p(eh)                              # pvmpv(earth_bary, earth_helio)
    pmoon = _add(_p(_moon98(tb1, tb2)), pe)              # pvppv(moon_geo, earth_bary)
    pjup = _add(_p(_plan94(tb1, tb2, 5)), psun)          # pvppv(jupiter_helio, sun_bary)
    # gravitational redshift (m/s): summed Earth, Moon, Jupiter, Sun (astropy's order)
    rk = 0x1.d662d0f88171ap-38                           # (m^2/s/AU).to(m/s) = 1/AU in m
    gr = 0.0 + (-GM_EARTH / C_M_S / _pm(site.xyz))
    gr += (-GM_MOON / C_M_S / _pm(_sub(pmoon, pe))) * rk
    gr += (-GM_JUP / C_M_S / _pm(_sub(pjup, pe))) * rk
    gr += (-GM_SUN / C_M_S / _pm(_sub(psun, pe))) * rk
    # beta in AU/d/(m/s) as astropy carries it; S converts it to dimensionless
    v = _add(ve, gvel .* MS_TO_AUD)
    b = v ./ C_M_S
    nb = _pm(b)
    gamma = 1 / sqrt(1 - (nb * nb) * AUD_MS2)
    den = 1 + gr / C_M_S
    return (; gamma, beta = b, den, gr, gpos, gvel, tdb = (tb1, tb2))
end

"""
    bary_kms(site, eop, mjd_utc, n)

One value of astropy's radial_velocity_correction (km/s) for the unit vector n, with astropy's
own operation order: zb = gamma * (1 + beta.n) / (1 + gr/c);  v = (zb - 1) c.
"""
bary_kms(site::Site, T::EOPTable, mjd::Float64, n::V3) = _vkms(bary_state(site, T, mjd), n)
@inline function _vkms(st, n::V3)
    zb = st.gamma * (1 + _pdp(st.beta, n) * AUD_TO_MS) / st.den
    return ((zb - 1) * C_M_S) * 0.001                    # m/s -> km/s: astropy scale 0.001
end

"""
    eop_coverage(eop, mjds) -> NamedTuple

Where the exposure times fall in the IERS table: counts before/after the table, in its predicted
part (>= predictive_mjd) and after the last Bulletin B row.
"""
function eop_coverage(T::EOPTable, mjds::AbstractVector{Float64})
    f = filter(isfinite, mjds)
    (min_mjd = minimum(f), max_mjd = maximum(f),
     n_before_table = count(<(T.first_mjd), f), n_after_table = count(>(T.last_mjd), f),
     n_predicted = count(>=(T.predictive_mjd), f), n_after_bulletin_b = count(>(T.last_bulletin_b_mjd), f))
end

"Per-exposure U for all harvested exposures (status ok, finite time), threaded; BaryBasis layout."
function bary_basis_julia(et::ExposureTimes, T::EOPTable)
    ek = Dict{Int,Vector{Int64}}(); U = Dict{Int,Matrix{Float64}}()
    for (t, _) in TELENAME
        sel = findall(j -> et.tele[j] == t && et.status[j] == "ok" && isfinite(et.mjd_mid_precise[j]), eachindex(et.tele))
        ek[t] = [expkey(et.tele[j], et.mjd[j], et.expnum[j]) for j in sel]
        Ut = Matrix{Float64}(undef, 4, length(sel))
        site = SITES[t]
        Threads.@threads for q in eachindex(sel)
            c = bary_coeffs(site, T, et.mjd_mid_precise[sel[q]])
            @inbounds for r in 1:4; Ut[r, q] = c[r]; end
        end
        U[t] = Ut
    end
    return ek, U
end
