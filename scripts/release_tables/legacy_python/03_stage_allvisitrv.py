#!/usr/bin/env python
"""03_stage_allvisitrv.py -- copy-and-swap: add `offtarget` to allVisitRV.

CHAIN POSITION: offtarget step 3 of 3; produces the final allVisitRV (25 datasets).
SOURCE: work/2026_09_16/offtarget_sidecar/s4_stage_columns.py (0.4.0). Logic unchanged for the
allVisitRV job; paths became arguments; the allExposure job was removed because the Julia
allExposure builder writes `offtarget` directly (same dtype, chunking, filters, description).

Copies the input file to --dst and ADDS exactly one dataset,
`offtarget` (int8, 0/1), matching the host file's own chunking and filter settings.
Nothing existing is touched.  The released files are opened read-only.
"""
import argparse, os, shutil, subprocess, sys
import h5py
import numpy as np


_ap = argparse.ArgumentParser()
_ap.add_argument("--src", required=True, help="allVisitRV with 24 columns")
_ap.add_argument("--flag", required=True, help="offtarget_visit_flag.h5 from step 2")
_ap.add_argument("--dst", required=True, help="output allVisitRV (25 columns)")
_args = _ap.parse_args()
with h5py.File(_args.src, "r") as _f:
    _n = _f["sdss_id"].shape[0]

JOBS = [
    dict(src=_args.src,
         dst=_args.dst,
         flag=_args.flag,
         n=_n, chunks=(65536,), shuffle=True, gzip=9,
         desc=("FPS off-target flag, constant within a visit (verified).  1 where the almanac fiber "
               "table for this visit's (telescope, mjd, config_id, fiber) has assigned == 1 and "
               "on_target == 0, i.e. the robotic positioner was assigned this target but did not "
               "reach it, so the spectra are of whatever sky the fiber did reach while still "
               "carrying the assigned target's identity.  Task 94's sdss_id 85604486 is the worked "
               "example: a field star 9 arcmin away at S/N 96, clean on every other flag.  0 for all "
               "plate-era visits -- the plate era has no assigned/on_target columns and cannot be "
               "flagged this way.  Source: pass1_dr21/almanac/allobs_57618_61230.h5.")),
]


def log(*a):
    print(*a); sys.stdout.flush()


def main():
    for j in JOBS:
        log(f"\n=== {os.path.basename(j['dst'])} ===")
        with h5py.File(j["flag"], "r") as f:
            flag = f["offtarget"][:].astype(np.int8)
        assert flag.shape == (j["n"],), flag.shape
        assert set(np.unique(flag)) <= {0, 1}
        log(f"flag loaded: {flag.shape[0]} rows, {int(flag.sum())} flagged")

        if os.path.exists(j["dst"]):
            os.remove(j["dst"])
        log("copying...")
        shutil.copyfile(j["src"], j["dst"])
        log(f"copied, {os.path.getsize(j['dst'])} bytes")

        with h5py.File(j["dst"], "r+") as f:
            if "offtarget" in f:
                raise SystemExit("refusing to clobber an existing offtarget dataset")
            d = f.create_dataset("offtarget", data=flag, dtype="i1",
                                 chunks=j["chunks"], shuffle=j["shuffle"],
                                 compression="gzip", compression_opts=j["gzip"])
            d.attrs["description"] = j["desc"]
        with h5py.File(j["dst"], "r") as f:
            d = f["offtarget"]
            log(f"added: dtype={d.dtype} shape={d.shape} chunks={d.chunks} "
                f"comp={d.compression}{d.compression_opts} shuffle={d.shuffle} sum={int(d[:].sum())}")
        log(f"final size {os.path.getsize(j['dst'])} bytes")
        log("md5: " + subprocess.run(["md5sum", j["dst"]], capture_output=True,
                                     text=True).stdout.split()[0])


if __name__ == "__main__":
    main()
