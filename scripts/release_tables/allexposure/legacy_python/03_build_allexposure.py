#!/usr/bin/env python
"""03_build_allexposure.py -- build the pass-1 DR21 allExposure summary.

CHAIN POSITION: allExposure step 3 of 5 (after 02_harvest_times.py, before 04_add_bary_flag.py).
INPUTS:  --identity  identity.h5 from step 1
         --almanac   the almanac file the workup was run on (raw/<obs>/<mjd>/{exposures,fibers})
         --workup    arMADGICS workup directory (arMADGICS_out_<key>.h5 for RV_pixoff_final,
                     RV_pix_var, RV_flag, ingestBit, snr)
         --times     exposure_times.h5 from step 2
OUTPUTS: --out       allExposure h5 with 44 datasets (the 0.4.0 table before bary_flag, glon, glat,
                     qa_bitmask, v_rad_corr*, offtarget were added)
         --report    build_report.json (counts and barycentric validation)
SOURCE:  work/2026_09_13/allexposure/s6_build_allexposure.py (0.4.0). Logic unchanged; the
         hard-coded paths became the arguments above.

Product = the almanac's per-exposure and per-fiber METADATA columns, joined row-for-row
onto the 26,715,540 arMADGICS pass-1 workup rows, plus BOTH velocity frames:

    v_rel                      observed / topocentric velocity, exactly as pass-1 produces it
    v_barycentric_correction   the barycentric correction
    v_rad                      barycentric-corrected velocity

The observed frame is first-class and is NOT overwritten: sky emission-line work needs it.

JOIN (the part that has to be right -- see s1/s3 for the verification):
    workup row identity  = (tele, mjd, expnum, adjfiberindx)   from identity.h5
    almanac exposure row = raw/<obs>/<mjd>/exposures  where exposure == expnum
    almanac fiber table  = raw/<obs>/<mjd>/fibers/<key>, key = config_id (fps era)
                                                          or plate_id  (plate era)
    fiber index          = 301 - (adjfiberindx - 300*(tele-1))      <-- REVERSED
The reversal was established, not assumed: it reproduces the almanac sdss_id on 100.00%
of rows where both sides carry one (s3_join_verify.jl).

BARYCENTRIC METHOD: astropy radial_velocity_correction(kind='barycentric'), following
almanac/src/almanac/postprocess/radial_velocity.py, evaluated at mjd_mid_exposure_precise
(harvested from the pass-1 apred 1D outputs by s4, since the almanac carries only the
integer night MJD).  It is evaluated EXACTLY but cheaply by exploiting the fact that

    v_bc(site, t, nhat) = c * (u0(site,t) + uvec(site,t) . nhat - 1)

is exactly linear in the target unit vector nhat -- verified out-of-sample to 4e-7 m/s.
So astropy is called on 6 basis directions per exposure (144,381 exposures) and the
26.7M per-row values follow from one numpy matvec.  A direct astropy recomputation on a
random row sample validates the result.

Velocity composition follows almanac's finalize_radial_velocities() exactly:
    z      = 10**(RV_pixoff_final * 6e-6) - 1
    z_corr = v_bary / c
    v_rel  = z_to_v(z)
    v_rad  = z_to_v(z + z_corr + z_corr*z)        i.e. (1+z_tot) = (1+z)(1+z_corr)
    e_v_rad= |dv/dz| * dz/dp * sqrt(RV_pix_var)

READ-ONLY on every pipeline input.  Writes only --out and --report.
"""
import argparse
import json
import os
import sys
import time
import warnings

import h5py
import numpy as np

warnings.filterwarnings("ignore")
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

_ap = argparse.ArgumentParser(description="build the allExposure summary (44 datasets)")
_ap.add_argument("--identity", required=True, help="identity.h5 from 01_build_identity.jl")
_ap.add_argument("--almanac", required=True, help="almanac allobs_<start>_<end>.h5 the workup used")
_ap.add_argument("--workup", required=True, help="arMADGICS workup directory (arMADGICS_out_*.h5)")
_ap.add_argument("--times", required=True, help="exposure_times.h5 from 02_harvest_times.py")
_ap.add_argument("--out", required=True, help="output allExposure h5")
_ap.add_argument("--report", required=True, help="output build_report.json")
_args = _ap.parse_args()
IDENT = _args.identity
ALM = _args.almanac
WK = _args.workup
TIMES = _args.times
OUT = _args.out
REPORT = _args.report

C_KM_S = 299792.458
DELTA_LOG = 6e-6
TELENAME = {1: "apo", 2: "lco"}
# The recovered reference pattern uses EarthLocation.of_site('apo') / ('lco'), and so does
# almanac postprocess/radial_velocity.py.  Use the same two registry keys verbatim.
ASTROPY_SITE = {1: "apo", 2: "lco"}

# per-fiber columns common to BOTH the plate-era and FPS-era almanac fiber tables
FIB_NUM = ["sdss_id", "catalogid", "ra", "dec", "lambda_eff", "x_focal", "y_focal"]
FIB_CAT = ["category", "hole_type"]
# per-exposure columns
EXP_NUM = ["config_id", "plate_id", "design_id", "field_id", "cart_id", "map_id",
           "n_read", "fps", "flagged_bad", "chip_flags", "seeing", "focus",
           "collpist", "colpitch", "dithered_pixels"]
EXP_CAT = ["image_type"]
# workup scalars that are the provenance of the velocities
WK_KEYS = ["RV_pixoff_final", "RV_pix_var", "RV_flag", "ingestBit", "snr"]


def say(s):
    print(s, flush=True)


def as_str(a):
    """bytes/other -> python str array."""
    a = np.asarray(a)
    if a.dtype.kind == "S":
        return np.char.decode(a, "utf-8", errors="replace")
    return a.astype(str)


def main():
    t_start = time.time()

    say("=== loading identity ===")
    with h5py.File(IDENT, "r") as f:
        tele = f["tele"][:].astype(np.int8)
        mjd = f["mjd"][:].astype(np.int32)
        expnum = f["expnum"][:].astype(np.int32)
        adjf = f["adjfiberindx"][:].astype(np.int16)
        sdid_wk = f["sdss_id"][:].astype(np.int64)
    N = tele.size
    say(f"N = {N}")

    say("=== loading harvested mid-exposure times ===")
    with h5py.File(TIMES, "r") as f:
        et_tele = f["tele"][:].astype(np.int8)
        et_mjd = f["mjd"][:].astype(np.int32)
        et_exp = f["expnum"][:].astype(np.int32)
        et_mid = f["mjd_mid_exposure_precise"][:]
        et_status = as_str(f["status"][:])
    say(f"exposures with a harvested time: {np.sum(et_status == 'ok')} / {et_tele.size}")
    tkey = (et_tele.astype(np.int64) * 10**12 + et_mjd.astype(np.int64) * 10**6
            + et_exp.astype(np.int64))
    order = np.argsort(tkey)
    tkey_s, et_mid_s = tkey[order], et_mid[order]

    # ---------------------------------------------------------------- allocate
    say("=== allocating output columns ===")
    out = {}
    out["tele"] = tele
    out["mjd"] = mjd
    out["expnum"] = expnum
    out["adjfiberindx"] = adjf
    out["fiber_index_almanac"] = (301 - (adjf.astype(np.int16) - 300 * (tele.astype(np.int16) - 1))).astype(np.int16)
    out["sdss_id_workup"] = sdid_wk
    for c in FIB_NUM:
        out[c] = np.full(N, np.nan) if c in ("ra", "dec", "lambda_eff", "x_focal", "y_focal") \
            else np.full(N, -1, dtype=np.int64)
    for c in EXP_NUM:
        out[c] = np.full(N, np.nan, dtype=np.float32) if c in ("seeing", "focus", "collpist", "colpitch", "dithered_pixels") \
            else np.full(N, -1, dtype=np.int32)
    cat_codes = {c: np.zeros(N, dtype=np.uint8) for c in FIB_CAT + EXP_CAT}
    cat_maps = {c: {"": 0} for c in FIB_CAT + EXP_CAT}
    out["join_status"] = np.zeros(N, dtype=np.uint8)
    # 0 unset, 1 ok, 2 no exposure row, 3 no fiber table, 4 wrong length, 5 idx range
    out["era"] = np.zeros(N, dtype=np.uint8)   # 1 plate, 2 fps
    out["sdss_id_agrees"] = np.zeros(N, dtype=np.uint8)  # 1 agree, 2 disagree, 3 untestable

    # ---------------------------------------------------------------- the join
    say("=== joining almanac metadata (3,657 telescope-nights) ===")
    idx_by_night = {}
    okey = tele.astype(np.int64) * 10**6 + mjd.astype(np.int64)
    o = np.argsort(okey, kind="stable")
    ok_s = okey[o]
    bounds = np.searchsorted(ok_s, np.unique(ok_s), side="left")
    uniq = np.unique(ok_s)
    ends = np.append(bounds[1:], len(ok_s))
    for uq, b, e in zip(uniq, bounds, ends):
        idx_by_night[(int(uq // 10**6), int(uq % 10**6))] = o[b:e]
    say(f"telescope-nights: {len(idx_by_night)}")

    t0 = time.time()
    with h5py.File(ALM, "r") as A:
        for ni, (t, night) in enumerate(sorted(idx_by_night)):
            rows = idx_by_night[(t, night)]
            g = A["raw"][TELENAME[t]][str(night)]
            ex = g["exposures"]
            e_exp = ex["exposure"][:].astype(np.int64)
            e_fps = ex["fps"][:].astype(np.int64)
            e_cfg = ex["config_id"][:].astype(np.int64)
            e_plt = ex["plate_id"][:].astype(np.int64)
            e_num = {c: (ex[c][:] if c in ex else None) for c in EXP_NUM}
            e_cat = {c: (as_str(ex[c][:]) if c in ex else None) for c in EXP_CAT}
            pos = {int(v): i for i, v in enumerate(e_exp)}
            fg = g["fibers"] if "fibers" in g else None
            fkeys = set(fg.keys()) if fg is not None else set()
            cache = {}

            erow = np.array([pos.get(int(v), -1) for v in expnum[rows]])
            miss = erow < 0
            if miss.any():
                out["join_status"][rows[miss]] = 2
            have = ~miss
            if not have.any():
                continue
            r_ok = rows[have]
            er = erow[have]
            era = np.where(e_fps[er] == 1, 2, 1).astype(np.uint8)
            out["era"][r_ok] = era
            for c in EXP_NUM:
                if e_num[c] is not None:
                    out[c][r_ok] = e_num[c][er]
            for c in EXP_CAT:
                if e_cat[c] is not None:
                    vals = e_cat[c][er]
                    m = cat_maps[c]
                    for vv in np.unique(vals):
                        if vv not in m:
                            m[vv] = len(m)
                    cat_codes[c][r_ok] = np.array([m[v] for v in vals], dtype=np.uint8)

            keystr = np.where(e_fps[er] == 1, e_cfg[er].astype(str), e_plt[er].astype(str))
            fidx = out["fiber_index_almanac"][r_ok].astype(np.int64) - 1   # 0-based
            for kk in np.unique(keystr):
                sel = keystr == kk
                rr = r_ok[sel]
                if fg is None or kk not in fkeys:
                    out["join_status"][rr] = 3
                    continue
                if kk not in cache:
                    grp = fg[kk]
                    d = {}
                    for c in FIB_NUM:
                        d[c] = grp[c][:] if c in grp else None
                    for c in FIB_CAT:
                        d[c] = as_str(grp[c][:]) if c in grp else None
                    d["_n"] = grp["ra"].shape[0] if "ra" in grp else 0
                    cache[kk] = d
                d = cache[kk]
                if d["_n"] != 300:
                    out["join_status"][rr] = 4
                    continue
                ii = fidx[sel]
                bad = (ii < 0) | (ii > 299)
                if bad.any():
                    out["join_status"][rr[bad]] = 5
                gd = ~bad
                if not gd.any():
                    continue
                rr2, ii2 = rr[gd], ii[gd]
                for c in FIB_NUM:
                    if d[c] is not None:
                        out[c][rr2] = d[c][ii2]
                for c in FIB_CAT:
                    if d[c] is not None:
                        vals = d[c][ii2]
                        m = cat_maps[c]
                        for vv in np.unique(vals):
                            if vv not in m:
                                m[vv] = len(m)
                        cat_codes[c][rr2] = np.array([m[v] for v in vals], dtype=np.uint8)
                out["join_status"][rr2] = 1
                # sdss_id cross-check
                if d["sdss_id"] is not None:
                    a_s = d["sdss_id"][ii2].astype(np.int64)
                    w_s = sdid_wk[rr2]
                    testable = (a_s > 0) & (w_s > 0)
                    st = np.where(~testable, 3, np.where(a_s == w_s, 1, 2)).astype(np.uint8)
                    out["sdss_id_agrees"][rr2] = st
            if (ni + 1) % 400 == 0:
                say(f"  {ni+1}/{len(idx_by_night)} nights  {time.time()-t0:.0f}s")
    say(f"join done in {time.time()-t0:.0f}s")

    # attach mid-exposure time per row
    say("=== attaching mid-exposure times ===")
    rk = (tele.astype(np.int64) * 10**12 + mjd.astype(np.int64) * 10**6
          + expnum.astype(np.int64))
    p = np.searchsorted(tkey_s, rk)
    p = np.clip(p, 0, len(tkey_s) - 1)
    hit = tkey_s[p] == rk
    mjd_mid = np.full(N, np.nan)
    mjd_mid[hit] = et_mid_s[p[hit]]
    out["mjd_mid_exposure"] = mjd_mid
    say(f"rows with a mid-exposure time: {np.sum(np.isfinite(mjd_mid))} "
        f"({100*np.sum(np.isfinite(mjd_mid))/N:.4f}%)")

    # ---------------------------------------------------------------- barycentric
    say("=== barycentric correction (factorized exact astropy) ===")
    ra = out["ra"]; dec = out["dec"]
    v_bc = np.full(N, np.nan)
    # basis directions: octahedron -> design matrix [1, nx, ny, nz] is perfectly conditioned
    BAS = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0], [0, -1., 0], [0, 0, 1.], [0, 0, -1.]])
    bas_ra = np.degrees(np.arctan2(BAS[:, 1], BAS[:, 0])) % 360.0
    bas_dec = np.degrees(np.arcsin(np.clip(BAS[:, 2], -1, 1)))
    Amat = np.hstack([np.ones((6, 1)), BAS])
    Apinv = np.linalg.pinv(Amat)
    say(f"design matrix cond = {np.linalg.cond(Amat):.3f}")

    fit_report = {}
    # exposures that actually carry rows AND have a time
    for t in (1, 2):
        m_et = (et_tele == t) & (et_status == "ok") & np.isfinite(et_mid)
        exp_mid = et_mid[m_et]
        exp_key = tkey[m_et]
        ne = exp_mid.size
        if ne == 0:
            continue
        say(f"  {TELENAME[t]}: {ne} exposures x 6 directions = {6*ne} astropy evaluations")
        loc = EarthLocation.of_site(ASTROPY_SITE[t])
        Y = np.empty((6, ne))
        t1 = time.time()
        CH = 200000
        for d in range(6):
            vals = np.empty(ne)
            for s0 in range(0, ne, CH):
                s1 = min(s0 + CH, ne)
                sc = SkyCoord(ra=np.full(s1 - s0, bas_ra[d]) * u.deg,
                              dec=np.full(s1 - s0, bas_dec[d]) * u.deg)
                tt = Time(exp_mid[s0:s1], format="mjd", scale="utc", location=loc)
                vals[s0:s1] = sc.radial_velocity_correction(
                    kind="barycentric", obstime=tt).to(u.km / u.s).value
            Y[d] = vals / C_KM_S + 1.0
            say(f"     dir {d+1}/6  {time.time()-t1:.0f}s")
        U = Apinv @ Y                      # (4, ne)
        resid = Amat @ U - Y
        rmax = np.abs(resid).max() * C_KM_S * 1000.0
        say(f"     4-parameter model residual on the 6 basis directions: max {rmax:.3e} m/s")
        fit_report[TELENAME[t]] = dict(n_exposures=int(ne), basis_residual_m_s=float(rmax))

        # map each row to its exposure's U
        so = np.argsort(exp_key)
        eks = exp_key[so]
        ep = np.clip(np.searchsorted(eks, rk), 0, ne - 1)
        rowhit = (eks[ep] == rk) & (tele == t) & np.isfinite(ra) & np.isfinite(dec)
        ecol = so[ep[rowhit]]
        rr = np.radians(ra[rowhit]); dd = np.radians(dec[rowhit])
        nx = np.cos(dd) * np.cos(rr); ny = np.cos(dd) * np.sin(rr); nz = np.sin(dd)
        v_bc[rowhit] = C_KM_S * (U[0, ecol] + U[1, ecol] * nx + U[2, ecol] * ny
                                 + U[3, ecol] * nz - 1.0)
        say(f"     {TELENAME[t]}: filled {np.sum(rowhit)} rows")

    out["v_barycentric_correction"] = v_bc
    say(f"v_bary: finite {np.sum(np.isfinite(v_bc))} ({100*np.sum(np.isfinite(v_bc))/N:.4f}%)  "
        f"range {np.nanmin(v_bc):.3f} .. {np.nanmax(v_bc):.3f} km/s")

    # ---------------------------------------------------------------- velocities
    say("=== velocities (almanac finalize_radial_velocities convention) ===")
    wkv = {}
    for k in WK_KEYS:
        with h5py.File(f"{WK}/arMADGICS_out_{k}.h5", "r") as f:
            wkv[k] = f[k][:]
        say(f"  read {k}")
    p_off = wkv["RV_pixoff_final"].astype(np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        z = 10 ** (p_off * DELTA_LOG) - 1
        z_corr = v_bc / C_KM_S
        v_rel = ((z + 1) ** 2 - 1) / ((z + 1) ** 2 + 1) * C_KM_S
        z_tot = z + z_corr + z_corr * z
        v_rad = ((z_tot + 1) ** 2 - 1) / ((z_tot + 1) ** 2 + 1) * C_KM_S
        e_z = DELTA_LOG * np.log(10) * 10 ** (p_off * DELTA_LOG) * np.sqrt(wkv["RV_pix_var"].astype(np.float64))
        e_v = np.abs(4 * (z + 1) / (((z + 1) ** 2 + 1) ** 2)) * C_KM_S * e_z
    out["v_rel"] = v_rel
    out["v_rad"] = v_rad
    out["e_v_rad"] = e_v
    for k in WK_KEYS:
        out[k] = wkv[k]
    say(f"v_rel  finite {np.sum(np.isfinite(v_rel))}   v_rad finite {np.sum(np.isfinite(v_rad))}")

    # ---------------------------------------------------------------- validation
    say("=== validation: direct astropy on a random row sample ===")
    rng = np.random.default_rng(0)
    cand = np.flatnonzero(np.isfinite(v_bc))
    samp = rng.choice(cand, size=min(20000, cand.size), replace=False)
    dmax = 0.0; dall = []
    for t in (1, 2):
        s = samp[tele[samp] == t]
        if s.size == 0:
            continue
        loc = EarthLocation.of_site(ASTROPY_SITE[t])
        sc = SkyCoord(ra=ra[s] * u.deg, dec=dec[s] * u.deg)
        tt = Time(mjd_mid[s], format="mjd", scale="utc", location=loc)
        ref = sc.radial_velocity_correction(kind="barycentric", obstime=tt).to(u.km / u.s).value
        d = (v_bc[s] - ref) * 1000.0
        dall.append(d)
        say(f"  {TELENAME[t]}: n={s.size}  max|factorized - direct astropy| = {np.abs(d).max():.3e} m/s  "
            f"rms {np.sqrt((d**2).mean()):.3e} m/s")
        dmax = max(dmax, float(np.abs(d).max()))
    fit_report["validation_max_m_s"] = dmax

    # ---------------------------------------------------------------- write
    say("=== writing ===")
    if os.path.exists(OUT):
        os.remove(OUT)
    with h5py.File(OUT, "w") as f:
        for k, v in out.items():
            f.create_dataset(k, data=v, compression="gzip", compression_opts=1, shuffle=True)
        for c, codes in cat_codes.items():
            assert len(cat_maps[c]) <= 255, f"{c}: {len(cat_maps[c])} categories overflows uint8"
            ds = f.create_dataset(c, data=codes, compression="gzip", compression_opts=1)
            inv = [""] * len(cat_maps[c])
            for kk, vv in cat_maps[c].items():
                inv[vv] = kk
            ds.attrs["categories"] = np.array(inv, dtype=h5py.string_dtype())
        f.attrs["description"] = (
            "APOGEE DR21 pass-1 allExposure summary: almanac per-exposure and per-fiber "
            "metadata joined to the arMADGICS workup rows, with both the observed "
            "(topocentric) and barycentric velocity frames.")
        f.attrs["n_rows"] = N
        f.attrs["bary_method"] = (
            "astropy SkyCoord.radial_velocity_correction(kind='barycentric'), evaluated at "
            "mjd_mid_exposure_precise (UTC) with EarthLocation.of_site; computed via an exact "
            "linear factorization in the target unit vector (6 basis directions per exposure).")
        f.attrs["astropy_version"] = __import__("astropy").__version__
        f.attrs["velocity_convention"] = (
            "v_rel = observed/topocentric (NO barycentric correction applied); "
            "v_barycentric_correction = astropy barycentric correction, km/s; "
            "v_rad = barycentric-corrected, composed in redshift as (1+z_tot)=(1+z)(1+z_corr), "
            "matching almanac postprocess/radial_velocity.py finalize_radial_velocities().")
        f.attrs["fiber_index_convention"] = (
            "fiber_index_almanac = 301 - (adjfiberindx - 300*(tele-1)); REVERSED relative to "
            "adjfiberindx. Verified against almanac sdss_id on 100.00% of testable rows.")
        f.attrs["join_status_codes"] = "1=ok 2=no_exposure_row 3=no_fiber_table 4=wrong_length 5=idx_out_of_range"
        f.attrs["era_codes"] = "1=plate 2=fps"
        f.attrs["sdss_id_agrees_codes"] = "1=agree 2=disagree 3=untestable(one side <=0)"
        f.attrs["tele_codes"] = "1=apo 2=lco"
    with open(REPORT, "w") as f:
        json.dump(dict(
            n_rows=int(N),
            join_status={str(k): int(np.sum(out["join_status"] == k)) for k in range(6)},
            era={str(k): int(np.sum(out["era"] == k)) for k in range(3)},
            sdss_id_agrees={str(k): int(np.sum(out["sdss_id_agrees"] == k)) for k in range(4)},
            n_finite_vbary=int(np.sum(np.isfinite(v_bc))),
            n_finite_mid=int(np.sum(np.isfinite(mjd_mid))),
            bary=fit_report,
        ), f, indent=2)
    say(f"wrote {OUT}  ({os.path.getsize(OUT)/1e9:.2f} GB)")
    say(f"TOTAL {time.time()-t_start:.0f}s")
    say("DONE")


if __name__ == "__main__":
    main()
