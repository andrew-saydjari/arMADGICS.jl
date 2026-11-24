#!/usr/bin/env python
"""make_restrict_ids.py -- the sdss_id set that `almanac add metadata --restrict-to` queries.

CHAIN POSITION: metadata step 1 of 2 (after the allExposure build, before build_metadata.sh).
INPUT:  argv[1]  allExposure h5 (reads `sdss_id`)
OUTPUT: argv[2]  .npy of the unique sdss_id > 0, sorted ascending, int64
For 0.4.0 this was 4,019,650 ids (RUN_INFO.txt: "0.4.0 allExposure sdss_id > 0"). The 0.4.0 file
was made inline, without a kept script; this script reproduces it byte for byte (checked against
work/2026_11_18/targets_040/sdss_id_allexposure_040.npy in the regression).
NOTE for 0.4.1: sdss_id 0 (the solar-system target, Vesta) is excluded by `> 0` here.
"""
import sys

import h5py
import numpy as np

with h5py.File(sys.argv[1], "r") as f:
    s = f["sdss_id"][:].astype(np.int64)
ids = np.unique(s[s > 0])
np.save(sys.argv[2], ids)
print(f"{ids.size} unique sdss_id > 0 -> {sys.argv[2]}")
