#!/usr/bin/env python
"""04_add_bary_flag.py -- add the bary_flag quality column to the allExposure product.

CHAIN POSITION: allExposure step 4 of 5 (legacy Python chain; after 03, before 05).
INPUTS:  --times   exposure_times.h5 from step 2
         --prod    the allExposure h5 from step 3, MODIFIED IN PLACE (adds bary_flag + a root attr)
OUTPUTS: --report  flag_report.json
SOURCE:  work/2026_09_13/allexposure/s9_add_flags.py (0.4.0). Logic unchanged; hard-coded paths
         became the arguments above.

Nothing here is silently dropped: every row keeps its velocities, and the problems found
during the build are surfaced as an explicit bitmask so a consumer can cut on them.

bit 0 (1)  NO_BARY      v_barycentric_correction is not finite (no RA/Dec or no time)
bit 1 (2)  BAD_TIME     mjd_mid_exposure is implausible for its night MJD (outside
                        -0.6 .. +1.6 d).  Affects LCO 57800 exposures 1-8, whose
                        harvested mid-exposure time is frozen at 57580.122532, i.e.
                        219.88 d before the night -- a bad DATE-OBS in the raw headers
                        for LCO's very first night.  The barycentric correction on those
                        rows is evaluated at the wrong epoch and is WRONG BY UP TO THE
                        FULL +/-30 km/s AMPLITUDE.  Cut on this bit.
bit 2 (4)  ID_MISMATCH  the almanac sdss_id at the joined fiber disagrees with the workup
                        sdss_id, so the RA/Dec provenance for that row is uncertain.
                        1,012 rows, all APO plate-era, confined to MJD 58036-58049.
                        The targets involved sit on the same plate, so the induced
                        barycentric error is bounded by the plate radius: roughly
                        30 km/s * 1.5 deg = 0.8 km/s worst case.
bit 3 (8)  NO_RV        RV_pixoff_final is not finite, so v_rel/v_rad are not defined.

READ-ONLY on pipeline inputs; modifies only the --prod file.
"""
import argparse
import json

import h5py
import numpy as np

_ap = argparse.ArgumentParser(description="add bary_flag to allExposure in place")
_ap.add_argument("--times", required=True)
_ap.add_argument("--prod", required=True)
_ap.add_argument("--report", required=True)
_args = _ap.parse_args()
PROD = _args.prod
TIMES = _args.times
REPORT = _args.report

NO_BARY, BAD_TIME, ID_MISMATCH, NO_RV = 1, 2, 4, 8


def main():
    with h5py.File(TIMES, "r") as f:
        et = f["tele"][:]
        em = f["mjd"][:]
        ee = f["expnum"][:]
        emid = f["mjd_mid_exposure_precise"][:]
    off = emid - em
    bad_exp = (off < -0.6) | (off > 1.6) | ~np.isfinite(emid)
    bad_keys = (et[bad_exp].astype(np.int64) * 10**12
                + em[bad_exp].astype(np.int64) * 10**6
                + ee[bad_exp].astype(np.int64))
    print(f"exposures with implausible mid-exposure time: {bad_exp.sum()}")

    with h5py.File(PROD, "r+") as f:
        n = f.attrs["n_rows"]
        tele = f["tele"][:]
        mjd = f["mjd"][:]
        expnum = f["expnum"][:]
        vbc = f["v_barycentric_correction"][:]
        agree = f["sdss_id_agrees"][:]
        pix = f["RV_pixoff_final"][:]

        flag = np.zeros(n, dtype=np.uint8)
        flag[~np.isfinite(vbc)] |= NO_BARY
        rk = (tele.astype(np.int64) * 10**12 + mjd.astype(np.int64) * 10**6
              + expnum.astype(np.int64))
        flag[np.isin(rk, bad_keys)] |= BAD_TIME
        flag[agree == 2] |= ID_MISMATCH
        flag[~np.isfinite(pix)] |= NO_RV

        if "bary_flag" in f:
            del f["bary_flag"]
        ds = f.create_dataset("bary_flag", data=flag, compression="gzip",
                              compression_opts=1)
        ds.attrs["bits"] = ("1=NO_BARY (no barycentric correction available); "
                            "2=BAD_TIME (mid-exposure time implausible for the night; "
                            "LCO 57800 exp 1-8, correction evaluated at the wrong epoch, "
                            "wrong by up to the full +/-30 km/s); "
                            "4=ID_MISMATCH (almanac vs workup sdss_id disagree, RA/Dec "
                            "provenance uncertain, error bounded by plate radius ~0.8 km/s); "
                            "8=NO_RV (RV_pixoff_final not finite)")
        f.attrs["recommended_cut"] = "bary_flag & 0b1111 == 0 for barycentric science"

        rep = {
            "n_rows": int(n),
            "flag_clean": int(np.sum(flag == 0)),
            "NO_BARY": int(np.sum(flag & NO_BARY != 0)),
            "BAD_TIME": int(np.sum(flag & BAD_TIME != 0)),
            "ID_MISMATCH": int(np.sum(flag & ID_MISMATCH != 0)),
            "NO_RV": int(np.sum(flag & NO_RV != 0)),
        }
    print(json.dumps(rep, indent=2))
    for k, v in rep.items():
        if k != "n_rows":
            print(f"  {k:14s} {v:10d}  {100*v/rep['n_rows']:8.5f}%")
    with open(REPORT, "w") as f:
        json.dump(rep, f, indent=2)
    print("DONE")


if __name__ == "__main__":
    main()
