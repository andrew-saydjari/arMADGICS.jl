#!/usr/bin/env python
"""legacy_kernels.py -- the two numpy/astropy kernels the Julia allExposure builder delegates.

Used only with `--numerics legacy` (the default, which reproduces 0.4.0 bit for bit). Everything
else in the allExposure build is Julia. Two kernels stay in Python because their exact bits are
defined by library implementations that have no bit-identical Julia equivalent:

  bary-basis  Per exposure, the 4 coefficients U = (u0, ux, uy, uz) of the exact linear
              factorization  v_bc(site, t, nhat) = c * (u0 + u . nhat - 1)  of astropy's
              SkyCoord.radial_velocity_correction(kind='barycentric'). astropy evaluates it with
              ERFA plus its own IERS interpolation; U = pinv(A) @ Y uses numpy/BLAS dgemm. This is
              the code of 0.4.0 s6_build_allexposure.py lines 273-302, unchanged, run per telescope.
              Cost: 6 astropy evaluations per exposure (~150 s for 144k exposures).
  pow10       Per row, 10 ** (RV_pixoff_final * 6e-6). On AVX-512 CPUs numpy evaluates float64
              power with Intel SVML, which differs from correctly rounded pow (glibc, Julia) in
              ~7% of rows by 1 ulp. 0.4.0 was built on an AVX-512 node, so its v_rel / v_rad /
              e_v_rad bits carry SVML's rounding. NOTE: on a CPU without AVX-512 (e.g. AMD rome
              nodes) numpy itself takes a different path and the 0.4.0 Python chain would NOT
              reproduce 0.4.0; the same holds for this kernel.

INPUTS / OUTPUTS:
  bary-basis --times exposure_times.h5 --out bary_basis.h5 [--iers-a finals2000A.all]
      out: per telescope group "<tele>/" with exp_key (int64, tele*1e12+mjd*1e6+expnum, in
      exposure_times.h5 order), U (4 x ne float64), plus attrs astropy_version and the basis residual.
  pow10 --workup DIR --out pow10.h5
      out: dataset pow10 (float64, one per workup row)

SOURCE: kernels lifted verbatim from work/2026_09_13/allexposure/s6_build_allexposure.py (0.4.0).
Run with the Python environment in ../../requirements.txt.
"""
import argparse
import time
import warnings

import h5py
import numpy as np

warnings.filterwarnings("ignore")

C_KM_S = 299792.458
DELTA_LOG = 6e-6
TELENAME = {1: "apo", 2: "lco"}
ASTROPY_SITE = {1: "apo", 2: "lco"}


def as_str(a):
    a = np.asarray(a)
    if a.dtype.kind == "S":
        return np.char.decode(a, "utf-8", errors="replace")
    return a.astype(str)


def bary_basis(times, out, iers_a=""):
    import hashlib
    import astropy
    from astropy.utils import iers
    from astropy import units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    with h5py.File(times, "r") as f:
        et_tele = f["tele"][:].astype(np.int8)
        et_mjd = f["mjd"][:].astype(np.int32)
        et_exp = f["expnum"][:].astype(np.int32)
        et_mid = f["mjd_mid_exposure_precise"][:]
        et_status = as_str(f["status"][:])
    tkey = (et_tele.astype(np.int64) * 10**12 + et_mjd.astype(np.int64) * 10**6
            + et_exp.astype(np.int64))

    BAS = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0], [0, -1., 0], [0, 0, 1.], [0, 0, -1.]])
    bas_ra = np.degrees(np.arctan2(BAS[:, 1], BAS[:, 0])) % 360.0
    bas_dec = np.degrees(np.arcsin(np.clip(BAS[:, 2], -1, 1)))
    Amat = np.hstack([np.ones((6, 1)), BAS])
    Apinv = np.linalg.pinv(Amat)

    # Earth-orientation table: astropy's IERS_Auto by default (what 0.4.0 used; its table state at
    # run time is NOT recorded in 0.4.0), or a pinned IERS-A file. Either way the table used is
    # recorded so the barycentric bits can be reproduced later.
    if iers_a:
        iers.earth_orientation_table.set(iers.IERS_A.open(iers_a))
    tab = iers.earth_orientation_table.get()
    src = str(tab.meta.get("data_path", "") or iers_a or "unknown")
    try:
        sha = hashlib.sha256(open(src, "rb").read()).hexdigest()
    except OSError:
        sha = ""
    with h5py.File(out, "w") as g:
        g.attrs["astropy_version"] = astropy.__version__
        g.attrs["iers_table"] = src
        g.attrs["iers_table_sha256"] = sha
        for t in (1, 2):
            m_et = (et_tele == t) & (et_status == "ok") & np.isfinite(et_mid)
            exp_mid = et_mid[m_et]
            exp_key = tkey[m_et]
            ne = exp_mid.size
            grp = g.create_group(TELENAME[t])
            grp["exp_key"] = exp_key
            if ne == 0:
                grp["U"] = np.zeros((4, 0))
                continue
            print(f"  {TELENAME[t]}: {ne} exposures x 6 directions", flush=True)
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
                print(f"     dir {d+1}/6  {time.time()-t1:.0f}s", flush=True)
            U = Apinv @ Y
            resid = Amat @ U - Y
            rmax = np.abs(resid).max() * C_KM_S * 1000.0
            grp["U"] = U
            grp.attrs["n_exposures"] = int(ne)
            grp.attrs["basis_residual_m_s"] = float(rmax)
            print(f"     basis residual max {rmax:.3e} m/s", flush=True)


def pow10(workup, out):
    with h5py.File(f"{workup}/arMADGICS_out_RV_pixoff_final.h5", "r") as f:
        p_off = f["RV_pixoff_final"][:].astype(np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        y = 10 ** (p_off * DELTA_LOG)
    with h5py.File(out, "w") as g:
        g["pow10"] = y


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("bary-basis")
    a.add_argument("--times", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--iers-a", default="", help="pin an IERS-A finals2000A file (default: astropy IERS_Auto)")
    b = sub.add_parser("pow10")
    b.add_argument("--workup", required=True)
    b.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()
    if args.cmd == "bary-basis":
        bary_basis(args.times, args.out, args.iers_a)
    else:
        pow10(args.workup, args.out)
    print(f"{args.cmd} done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
