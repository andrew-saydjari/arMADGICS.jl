#!/usr/bin/env python
"""recommended_cuts.py -- the default quality mask of the released allExposure table.

Definition = data/work/mwm/0.4.0/public/RECOMMENDED_CUTS.md ("The recipe"):
    qa_bitmask == 0, ingestBit == 0, RV_flag == 0, bary_flag == 0, isfinite(v_rad), offtarget == 0
`include_offtarget=False` gives the 5-cut mask that allVisitRV membership and the fiber-offset fit
used (0.4.0 built those before the offtarget column existed).

SOURCE: work/2026_09_16/public_exclusions/archive_pre_aks_edit/recommended_cuts.py (0.4.0 draft,
5 cuts, hard-coded path), extended with the offtarget cut that the released RECOMMENDED_CUTS.md
added on 2026-09-16 (offtarget_sidecar/recommended_cuts_patch.md). The path is an argument.

Usage: python recommended_cuts.py ALLEXPOSURE.h5 [--expect N5 N6]
  prints the kept counts; with --expect, exits non-zero unless they equal N5 (5 cuts) and N6 (6 cuts).
  0.4.0: N5 = 25516821, N6 = 25207025.
"""
import argparse
import sys

import h5py
import numpy as np


def recommended_mask(path, include_offtarget=True):
    """Boolean keep-mask (True = keep), one entry per allExposure row."""
    with h5py.File(path, "r") as f:
        keep = ((f["qa_bitmask"][:] == 0) & (f["ingestBit"][:] == 0) & (f["RV_flag"][:] == 0)
                & (f["bary_flag"][:] == 0) & np.isfinite(f["v_rad"][:]))
        if include_offtarget:
            keep &= f["offtarget"][:] == 0
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("allexposure")
    ap.add_argument("--expect", nargs=2, type=int)
    a = ap.parse_args()
    k5 = int(recommended_mask(a.allexposure, False).sum())
    k6 = int(recommended_mask(a.allexposure, True).sum())
    print(f"kept without offtarget cut: {k5}\nkept with offtarget cut   : {k6}")
    if a.expect and (k5, k6) != tuple(a.expect):
        print(f"MISMATCH: expected {a.expect}")
        sys.exit(1)


if __name__ == "__main__":
    main()
