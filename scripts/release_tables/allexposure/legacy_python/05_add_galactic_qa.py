#!/usr/bin/env python
"""05_add_galactic_qa.py -- add Galactic coordinates and the exclusion bitmask in place.

CHAIN POSITION: allExposure step 5 of 5 (legacy Python chain; after 04 and after
                ../../exclusions/build_exclusion_lists.jl).
INPUTS:  --prod     the allExposure h5 from step 4, MODIFIED IN PLACE (adds glon, glat, qa_bitmask)
         --lists    directory holding bad_nights.csv, bad_fibers.csv, rv_railed.csv
OUTPUTS: --report   qa_report.json
OPTIONS: --expect-rows / --expect-railed   the 0.4.0 sanity asserts (26,715,540 rows, 393 railed);
                    pass -1 to skip for a different release
SOURCE:  work/2026_09_13/allexposure/s10_add_galactic_qa.py (0.4.0). Logic unchanged; hard-coded
         paths and the two release-specific assert constants became the arguments above.

Two approved additions to the allExposure product, and nothing else:

  1. glon, glat  (float64, degrees)  ICRS (ra, dec) -> Galactic (l, b), computed with
     astropy.  Rows whose ra or dec is not finite get NaN.

  2. qa_bitmask  (int32)  built from the three public exclusion lists under
     --lists (read-only):

       bit 0 (1)  BAD_NIGHT        (tele, mjd) in bad_nights.csv
       bit 1 (2)  CHRONIC_FIBER    (tele, fiber) in a scope=chronic row of bad_fibers.csv
       bit 2 (4)  BAD_FIBER_NIGHT  (tele, mjd, fiber) in a scope=fiber_night row
       bit 3 (8)  BAD_FIBER_RANGE  (tele, fiber) in a scope=fiber_mjd_range row AND
                                   mjd_start <= mjd <= mjd_end
       bit 4 (16) RV_RAILED        (tele, mjd, expnum, adjfiberindx) in rv_railed.csv

FIBER CONVENTION.  bad_fibers.csv `fiber` is the 1-based fiber index *as it appears in the
product arrays*, 1-300 per telescope (public/README.md conventions: "For LCO rows,
fiber = adjfiberindx - 300").  That is `adjfiberindx - 300*(tele-1)`, NOT
`fiber_index_almanac`, which the product attrs document as REVERSED
(fiber_index_almanac = 301 - (adjfiberindx - 300*(tele-1))).  Verified empirically: the
three lowest-median-S/N APO product-index fibers (28, 90, 276) and the lowest-median-S/N
LCO fiber (170) are all in the WAVECAL_CHRONIC list under this convention, and are not
under the reversed one.

TELESCOPE CODING.  Product `tele` is int8 with 1=apo, 2=lco (product attr `tele_codes`);
the CSVs use the strings `apo`/`lco`.  Mapped explicitly.

READ-ONLY on the exclusion lists; modifies only the --prod file, and only by ADDING the three new datasets.  Nothing is deleted or rewritten.
"""
import argparse
import csv
import json
import sys

import h5py
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

_ap = argparse.ArgumentParser(description="add glon, glat, qa_bitmask to allExposure in place")
_ap.add_argument("--prod", required=True)
_ap.add_argument("--lists", required=True)
_ap.add_argument("--report", required=True)
_ap.add_argument("--expect-rows", type=int, default=26_715_540)
_ap.add_argument("--expect-railed", type=int, default=393)
_args = _ap.parse_args()
PROD = _args.prod
PUB = _args.lists

TELE_CODE = {"apo": 1, "lco": 2}

BAD_NIGHT, CHRONIC_FIBER, BAD_FIBER_NIGHT, BAD_FIBER_RANGE, RV_RAILED = 1, 2, 4, 8, 16

CHUNK = 2_000_000


def log(*a):
    print(*a)
    sys.stdout.flush()


def main():
    report = {}

    # ---------------------------------------------------------------- read keys
    log("=== reading identity columns ===")
    with h5py.File(PROD, "r") as f:
        for name in ("glon", "glat", "qa_bitmask"):
            if name in f:
                raise SystemExit(f"refusing to clobber existing dataset {name!r}")
        n = f["ra"].shape[0]
        ra = f["ra"][:]
        dec = f["dec"][:]
        tele = f["tele"][:].astype(np.int64)
        mjd = f["mjd"][:].astype(np.int64)
        expnum = f["expnum"][:].astype(np.int64)
        adjfib = f["adjfiberindx"][:].astype(np.int64)
    log(f"N = {n}")
    assert _args.expect_rows < 0 or n == _args.expect_rows, n

    # product-array fiber index, 1-300 per telescope
    fiber = adjfib - 300 * (tele - 1)
    assert fiber.min() == 1 and fiber.max() == 300

    # ------------------------------------------------------------ glon / glat
    log("=== ICRS -> Galactic (astropy) ===")
    glon = np.full(n, np.nan, dtype=np.float64)
    glat = np.full(n, np.nan, dtype=np.float64)
    good = np.isfinite(ra) & np.isfinite(dec)
    log(f"rows with finite ra/dec: {int(good.sum())} / {n}  (NaN rows: {int((~good).sum())})")
    idx = np.flatnonzero(good)
    for s in range(0, idx.size, CHUNK):
        sl = idx[s:s + CHUNK]
        c = SkyCoord(ra=ra[sl] * u.deg, dec=dec[sl] * u.deg, frame="icrs").galactic
        glon[sl] = c.l.deg
        glat[sl] = c.b.deg
        log(f"  {min(s + CHUNK, idx.size)}/{idx.size}")
    report["glon_range"] = [float(np.nanmin(glon)), float(np.nanmax(glon))]
    report["glat_range"] = [float(np.nanmin(glat)), float(np.nanmax(glat))]
    report["n_radec_nan"] = int((~good).sum())

    # ------------------------------------------------------------ qa_bitmask
    log("=== building qa_bitmask ===")
    qa = np.zeros(n, dtype=np.int32)

    # key helpers.  mjd <= 61230 < 1e5, expnum <= 243 < 1e3, fiber <= 300 < 1e3.
    def k_tn(t, m):            # (tele, mjd)
        return t * 100_000 + m

    def k_tf(t, fb):           # (tele, fiber)
        return t * 1_000 + fb

    def k_tnf(t, m, fb):       # (tele, mjd, fiber)
        return (t * 100_000 + m) * 1_000 + fb

    def k_tnef(t, m, e, af):   # (tele, mjd, expnum, adjfiberindx)
        return ((t * 100_000 + m) * 1_000 + e) * 1_000 + af

    # -- bit 0: bad nights
    bn = list(csv.DictReader(open(f"{PUB}/bad_nights.csv")))
    bn_keys = np.array(sorted({k_tn(TELE_CODE[r["tele"]], int(r["mjd"])) for r in bn}),
                       dtype=np.int64)
    log(f"bad_nights.csv rows={len(bn)} distinct (tele,mjd)={bn_keys.size}")
    hit0 = np.isin(k_tn(tele, mjd), bn_keys)
    qa[hit0] |= BAD_NIGHT

    # -- bits 1-3: bad fibers
    bf = list(csv.DictReader(open(f"{PUB}/bad_fibers.csv")))
    chronic = [r for r in bf if r["scope"] == "chronic"]
    fnight = [r for r in bf if r["scope"] == "fiber_night"]
    frange = [r for r in bf if r["scope"] == "fiber_mjd_range"]
    other = [r for r in bf if r["scope"] not in ("chronic", "fiber_night", "fiber_mjd_range")]
    assert not other, f"unexpected scope values: {{r['scope'] for r in other}}"
    log(f"bad_fibers.csv rows={len(bf)}  chronic={len(chronic)} "
        f"fiber_night={len(fnight)} fiber_mjd_range={len(frange)}")

    ch_keys = np.array(sorted({k_tf(TELE_CODE[r["tele"]], int(r["fiber"])) for r in chronic}),
                       dtype=np.int64)
    hit1 = np.isin(k_tf(tele, fiber), ch_keys)
    qa[hit1] |= CHRONIC_FIBER

    fn_keys = np.array(sorted({k_tnf(TELE_CODE[r["tele"]], int(r["mjd"]), int(r["fiber"]))
                               for r in fnight}), dtype=np.int64)
    hit2 = np.isin(k_tnf(tele, mjd, fiber), fn_keys)
    qa[hit2] |= BAD_FIBER_NIGHT

    hit3 = np.zeros(n, dtype=bool)
    tf_all = k_tf(tele, fiber)
    for r in frange:
        key = k_tf(TELE_CODE[r["tele"]], int(r["fiber"]))
        hit3 |= (tf_all == key) & (mjd >= int(r["mjd_start"])) & (mjd <= int(r["mjd_end"]))
    qa[hit3] |= BAD_FIBER_RANGE

    # -- bit 4: rv railed
    rv = list(csv.DictReader(open(f"{PUB}/rv_railed.csv")))
    rv_keys = np.array(sorted({k_tnef(TELE_CODE[r["tele"]], int(r["mjd"]), int(r["expnum"]),
                                      int(r["adjfiberindx"])) for r in rv}), dtype=np.int64)
    log(f"rv_railed.csv rows={len(rv)} distinct keys={rv_keys.size}")
    hit4 = np.isin(k_tnef(tele, mjd, expnum, adjfib), rv_keys)
    qa[hit4] |= RV_RAILED

    counts = {
        "bit0_BAD_NIGHT": int(hit0.sum()),
        "bit1_CHRONIC_FIBER": int(hit1.sum()),
        "bit2_BAD_FIBER_NIGHT": int(hit2.sum()),
        "bit3_BAD_FIBER_RANGE": int(hit3.sum()),
        "bit4_RV_RAILED": int(hit4.sum()),
        "any_nonzero": int((qa != 0).sum()),
    }
    for k, v in counts.items():
        log(f"  {k:24s} {v:,}")
    report["bit_counts"] = counts

    # -- independent cross-checks
    # bit 0: count rows on the 141 bad tele-nights with a plain per-night loop
    indep0 = 0
    for r in bn:
        indep0 += int(((tele == TELE_CODE[r["tele"]]) & (mjd == int(r["mjd"]))).sum())
    log(f"cross-check bit0: per-night loop = {indep0:,}  (vectorized = {counts['bit0_BAD_NIGHT']:,})")
    assert indep0 == counts["bit0_BAD_NIGHT"]
    report["bit0_independent_count"] = indep0
    assert _args.expect_railed < 0 or counts["bit4_RV_RAILED"] == _args.expect_railed, counts["bit4_RV_RAILED"]

    # ------------------------------------------------------------ write
    log("=== writing datasets in place ===")
    with h5py.File(PROD, "r+") as f:
        f.create_dataset("glon", data=glon)
        f.create_dataset("glat", data=glat)
        f.create_dataset("qa_bitmask", data=qa)
        f.attrs["galactic_convention"] = (
            "glon/glat are Galactic longitude/latitude in degrees, from ICRS (ra, dec) via "
            "astropy SkyCoord(...).galactic; NaN where ra or dec is not finite"
        )
        f.attrs["qa_bitmask_codes"] = (
            "bit0(1)=BAD_NIGHT (tele,mjd) in public/bad_nights.csv; "
            "bit1(2)=CHRONIC_FIBER (tele,fiber) scope=chronic in public/bad_fibers.csv; "
            "bit2(4)=BAD_FIBER_NIGHT (tele,mjd,fiber) scope=fiber_night; "
            "bit3(8)=BAD_FIBER_RANGE (tele,fiber) scope=fiber_mjd_range with "
            "mjd_start<=mjd<=mjd_end; bit4(16)=RV_RAILED (tele,mjd,expnum,adjfiberindx) in "
            "public/rv_railed.csv. fiber = adjfiberindx - 300*(tele-1) (product-array "
            "1-based fiber index, NOT fiber_index_almanac)."
        )
        f.attrs["qa_bitmask_source"] = "pass1_dr21/arMADGICS/public/ exclusion lists, 2026-09-16"

    with open(_args.report, "w") as fh:
        json.dump(report, fh, indent=2)
    log("done")


if __name__ == "__main__":
    main()
