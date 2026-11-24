#!/usr/bin/env python
"""02_harvest_times.py -- harvest the PRECISE MID-EXPOSURE TIME for every workup exposure.

CHAIN POSITION: allExposure step 2 of 5 (after 01_build_identity.jl, before 03_build_allexposure.py).
INPUTS:  --identity  identity.h5 from step 1
         --apred     pass-1 apred 1D output root (contains <mjd>/ar1Dunical_<tele>_<mjd>_<expnum>_object.h5)
OUTPUTS: --out       exposure_times.h5 (tele, mjd, expnum, mjd_mid_exposure, mjd_mid_exposure_precise,
                     nread_total, source, status), one row per distinct workup exposure
OPTIONS: --nproc     worker processes (0.4.0 used 12)
SOURCE:  work/2026_09_13/allexposure/s4_harvest_times.py (0.4.0). Logic unchanged; the hard-coded
         paths and process count became the arguments above.

The almanac carries only the integer night MJD.  A barycentric correction evaluated at a
night boundary would be wrong by up to ~12 h of Earth rotation plus orbital motion, i.e.
O(km/s) -- four orders above the 60 m/s pass-1 RV floor.  The precise time DOES exist in
the pass-1 apred 1D outputs:

    ar1Dunical_<tele>_<mjd>_<expnum:04d>_object.h5 : metadata/mjd_mid_exposure_precise

which is exactly the file and key almanac's own postprocess loader uses
(almanac/src/almanac/postprocess/loaders.py:19,53-65).  This script harvests it for each
distinct (tele, mjd, expnum) workup exposure.

Accuracy budget: d(v_bary)/dt is bounded by the diurnal term, ~3.3e-5 km/s per second, so
even a 30 s timing error is ~1 m/s -- far below the 60 m/s floor.  A few-second-accurate
mid-exposure time is therefore ample.

Why Python and not Julia: HDF5.jl serializes all file access on a global lock, so
Threads.@threads gives no parallelism on 144k small file opens.  A process pool does.

READ-ONLY on all inputs.  Writes only the --out file.
"""
import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import h5py
import numpy as np

_ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
_ap.add_argument("--identity", required=True, help="identity.h5 from 01_build_identity.jl")
_ap.add_argument("--apred", required=True, help="pass-1 apred 1D output root (per-MJD subdirectories)")
_ap.add_argument("--out", required=True, help="output exposure_times.h5")
_ap.add_argument("--nproc", type=int, default=12, help="worker processes (0.4.0 used 12)")
_args = _ap.parse_args()
IDENT = _args.identity
APRED = _args.apred
OUT = _args.out
TELENAME = {1: "apo", 2: "lco"}
# ar1Dunical is one file per exposure (no chip axis) and is what almanac's loader reads.
# ar1Dcal is per chip; keep it as a fallback so a missing unical file is not fatal.
CHIPS = ["a", "b", "c", "R", "G", "B"]
NPROC = _args.nproc


def do_night(arg):
    """Harvest every requested exposure of one (tele, night).  Runs in a worker process."""
    t, night, exps = arg
    tn = TELENAME[t]
    d = os.path.join(APRED, str(night))
    res = []
    if not os.path.isdir(d):
        return [(t, night, e, np.nan, np.nan, 0, "", "no_night_dir") for e in exps]
    for e in exps:
        mid = midp = np.nan
        nread = 0
        src = ""
        status = "no_object_file"
        cands = [(os.path.join(d, f"ar1Dunical_{tn}_{night}_{e:04d}_object.h5"), "unical")]
        cands += [(os.path.join(d, f"ar1Dcal_{tn}_{night}_{e:04d}_{c}_object.h5"), f"cal_{c}")
                  for c in CHIPS]
        for path, tag in cands:
            if not os.path.isfile(path):
                continue
            try:
                with h5py.File(path, "r") as h:
                    g = h["metadata"]
                    if "mjd_mid_exposure" in g:
                        mid = float(g["mjd_mid_exposure"][()])
                    if "mjd_mid_exposure_precise" in g:
                        midp = float(g["mjd_mid_exposure_precise"][()])
                    if "nread_total" in g:
                        nread = int(g["nread_total"][()])
                src, status = tag, "ok"
            except Exception:
                status = "read_error"
            if status == "ok":
                break
        res.append((t, night, e, mid, midp, nread, src, status))
    return res


def main():
    t_start = time.time()
    print("loading identity ...", flush=True)
    with h5py.File(IDENT, "r") as f:
        tele = f["tele"][:].astype(np.int8)
        mjd = f["mjd"][:].astype(np.int32)
        expnum = f["expnum"][:].astype(np.int32)
    print(f"N = {tele.size}", flush=True)

    key = (tele.astype(np.int64) * 10**12 + mjd.astype(np.int64) * 10**6
           + expnum.astype(np.int64))
    uk = np.unique(key)
    ut = (uk // 10**12).astype(np.int8)
    um = ((uk // 10**6) % 10**6).astype(np.int32)
    ue = (uk % 10**6).astype(np.int32)
    NE = uk.size
    print(f"distinct exposures: {NE}", flush=True)

    # group by (tele, night) so each task is one directory -> good locality on ceph
    tasks = {}
    for t, m, e in zip(ut, um, ue):
        tasks.setdefault((int(t), int(m)), []).append(int(e))
    tasklist = [(t, n, sorted(v)) for (t, n), v in sorted(tasks.items())]
    print(f"tasks (telescope-nights): {len(tasklist)}  with {NPROC} processes", flush=True)

    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=NPROC) as ex:
        for i, r in enumerate(ex.map(do_night, tasklist, chunksize=4)):
            rows.extend(r)
            if (i + 1) % 300 == 0:
                el = time.time() - t0
                print(f"  {i+1}/{len(tasklist)} nights  {len(rows)} exposures  {el:.0f}s  "
                      f"(eta {el/(i+1)*(len(tasklist)-i-1):.0f}s)", flush=True)
    print(f"harvest done in {time.time()-t0:.0f}s", flush=True)

    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    a_t = np.array([r[0] for r in rows], dtype=np.int8)
    a_m = np.array([r[1] for r in rows], dtype=np.int32)
    a_e = np.array([r[2] for r in rows], dtype=np.int32)
    a_mid = np.array([r[3] for r in rows], dtype=np.float64)
    a_midp = np.array([r[4] for r in rows], dtype=np.float64)
    a_nr = np.array([r[5] for r in rows], dtype=np.int32)
    a_src = np.array([r[6] for r in rows], dtype=object)
    a_st = np.array([r[7] for r in rows], dtype=object)
    assert a_t.size == NE, f"{a_t.size} != {NE}"

    print("\n=== harvest status ===", flush=True)
    for s in sorted(set(a_st)):
        n = int(np.sum(a_st == s))
        print(f"  {s:16s} {n:8d}  {100*n/NE:7.4f}%")
    print("sources: " + "  ".join(f"{s}={int(np.sum(a_src==s))}" for s in sorted(set(a_src)) if s))

    good = a_st == "ok"
    fin = good & np.isfinite(a_midp)
    print(f"\nmjd_mid_exposure_precise finite: {int(fin.sum())} / {NE}")
    if fin.any():
        print(f"  range {a_midp[fin].min():.6f} .. {a_midp[fin].max():.6f}")
        off = a_midp[fin] - a_m[fin]
        print(f"  mid-exposure minus night MJD: min {off.min():.4f}  max {off.max():.4f}  "
              f"median {np.median(off):.4f} d")
        nb = int(np.sum((off < -0.6) | (off > 1.6)))
        print(f"  implausible night offsets (<-0.6 or >1.6 d): {nb}")
        d = np.abs(a_midp[fin] - a_mid[fin])
        d = d[np.isfinite(d)]
        if d.size:
            print(f"  |precise - default mjd_mid| max = {d.max():.3e} d = {86400*d.max():.2f} s")

    with h5py.File(OUT, "w") as f:
        f["tele"] = a_t
        f["mjd"] = a_m
        f["expnum"] = a_e
        f["mjd_mid_exposure"] = a_mid
        f["mjd_mid_exposure_precise"] = a_midp
        f["nread_total"] = a_nr
        f.create_dataset("source", data=np.array([s.encode() for s in a_src]))
        f.create_dataset("status", data=np.array([s.encode() for s in a_st]))
    print(f"\nwrote {OUT}  ({NE} exposures, {int(good.sum())} ok)")
    print(f"TOTAL {time.time()-t_start:.0f}s")
    print("DONE")


if __name__ == "__main__":
    main()
