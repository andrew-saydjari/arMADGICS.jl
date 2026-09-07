"""Build the FULL E3 domeflat ar1Dunical input lists (both telescopes) with a
minimum-flux hygiene cut.

Extends 2026_09_01/telluric_refit/build_list.py (subset builder) to the full census:

  1. Bulk almanac census: 2026_05_01/outdir/almanac/allobs_57600_61160.h5,
     raw/<tele>/<mjd>/exposures with image_type == "domeflat".
  2. Existence check: 2026_05_01/outdir/apred/<mjd>/ar1Dunical_<tele>_<mjd>_<nnnn>_domeflat.h5
  3. Per-file screen (parallel, 4 workers):
       - schema: flux_1d/ivar_1d/mask_1d present, flux_1d shape == (300, 8700)
       - MIN-FLUX cut: masked median flux (every-10th-fiber subsample, mask_1d & ivar>0)
         must be >= MIN_MEDIAN_FLUX counts.  E3_REPORT.md finding: lco_60855_0011 is a
         near-zero-flux exposure (median ~1 count) labeled domeflat in the almanac; both
         fit implementations produce garbage on it.  Real domeflats have median flux
         O(10^3-10^4) counts; junk shutter/dark-like exposures sit at O(1).  Threshold
         100 counts sits >1 dex from both populations (verified in the printed histogram).

Outputs (in CWD):
  list_apo_full.txt / list_lco_full.txt / list_full.txt   accepted file paths (APO then LCO,
                                                          sorted by (mjd, exp) within tele)
  excluded_files.txt                                      excluded path, reason, median flux
  flux_medians.tsv                                        path <tab> masked-median-flux for
                                                          every existing candidate (record)
  list_stats.txt                                          summary counts

Run:  nice -n 10 uv run python build_full_list.py   (ccalin051; IO + light CPU only)
"""
import re
import sys
import numpy as np
import h5py
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

ALMANAC = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/almanac/allobs_57600_61160.h5"
APRED = Path("/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/apred")
MIN_MEDIAN_FLUX = 100.0     # counts; see docstring
EXPECTED_SHAPE = (300, 8700)
N_WORKERS = 4               # ccalin051 shared-node etiquette


def screen_file(fp):
    """Return (fp, ok_schema, median_flux). median_flux is np.nan if unreadable."""
    try:
        with h5py.File(fp, "r") as f:
            if not all(k in f for k in ("flux_1d", "ivar_1d", "mask_1d")):
                return fp, False, np.nan
            if f["flux_1d"].shape != EXPECTED_SHAPE:
                return fp, False, np.nan
            flux = f["flux_1d"][::10, :]
            ivar = f["ivar_1d"][::10, :]
            mask = f["mask_1d"][::10, :].astype(bool) & (ivar > 0)
        med = float(np.median(flux[mask])) if mask.any() else 0.0
        return fp, True, med
    except Exception:
        return fp, False, np.nan


def main():
    # 1) almanac census
    alm = {}
    with h5py.File(ALMANAC, "r") as f:
        for tele in ("apo", "lco"):
            g = f[f"raw/{tele}"]
            for mjd in g.keys():
                e = g[mjd]["exposures"]
                it = e["image_type"][:].astype(str)
                expstr = e["exposure_string"][:].astype(str)
                sel = it == "domeflat"
                if sel.any():
                    alm[(tele, int(mjd))] = [s[-4:] for s in expstr[sel]]
    n_alm = sum(len(v) for v in alm.values())
    print(f"almanac: {n_alm} domeflat exposures across {len(alm)} (tele,mjd) groups", flush=True)

    # 2) existence
    cands = {"apo": [], "lco": []}
    for (tele, mjd), exps in sorted(alm.items()):
        d = APRED / str(mjd)
        if not d.is_dir():
            continue
        for e4 in exps:
            fp = d / f"ar1Dunical_{tele}_{mjd}_{e4}_domeflat.h5"
            if fp.exists():
                cands[tele].append(str(fp))
    n_exist = {t: len(v) for t, v in cands.items()}
    print(f"existing ar1Dunical: apo={n_exist['apo']}  lco={n_exist['lco']}", flush=True)

    # 3) parallel schema + flux screen (incremental: appends to flux_medians.tsv
    #    as results arrive; on restart, resumes past already-screened paths)
    results = {}
    tsv = Path("flux_medians.tsv")
    if tsv.exists():
        for line in tsv.read_text().splitlines():
            try:
                fp, med, tag = line.split("\t")
                results[fp] = (tag == "ok", float(med))
            except ValueError:
                pass  # partial last line from an interrupted run
        print(f"resuming: {len(results)} already screened", flush=True)

    all_paths = cands["apo"] + cands["lco"]
    todo = [fp for fp in all_paths if fp not in results]
    with open(tsv, "a") as fm, ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for i, (fp, ok, med) in enumerate(ex.map(screen_file, todo, chunksize=8)):
            results[fp] = (ok, med)
            fm.write(f"{fp}\t{med:.3f}\t{'ok' if ok else 'bad_schema'}\n")
            if (i + 1) % 500 == 0:
                fm.flush()
                print(f"  screened {i+1}/{len(todo)}", flush=True)

    # 4) apply cuts, write lists
    kept = {"apo": [], "lco": []}
    excluded = []
    for tele in ("apo", "lco"):
        for fp in cands[tele]:
            ok, med = results[fp]
            if not ok:
                excluded.append((fp, "bad_schema", med))
            elif not np.isfinite(med) or med < MIN_MEDIAN_FLUX:
                excluded.append((fp, f"low_flux(<{MIN_MEDIAN_FLUX:g})", med))
            else:
                kept[tele].append(fp)

    for tele in ("apo", "lco"):
        with open(f"list_{tele}_full.txt", "w") as out:
            out.write("\n".join(kept[tele]) + "\n")
    with open("list_full.txt", "w") as out:
        out.write("\n".join(kept["apo"] + kept["lco"]) + "\n")
    with open("excluded_files.txt", "w") as out:
        for fp, reason, med in excluded:
            out.write(f"{fp}\t{reason}\tmedian_flux={med:.3f}\n")

    meds = np.array([results[fp][1] for fp in all_paths if results[fp][0]])
    lo = meds[np.isfinite(meds)]
    hist_edges = [0, 1, 10, 100, 1000, 10000, 100000, 1e9]
    hist, _ = np.histogram(lo, bins=hist_edges)
    junk = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/apred/60855/ar1Dunical_lco_60855_0011_domeflat.h5"
    junk_excluded = any(fp == junk for fp, _, _ in excluded)

    lines = [
        f"almanac domeflat exposures: {n_alm}",
        f"existing ar1Dunical: apo={n_exist['apo']} lco={n_exist['lco']} total={sum(n_exist.values())}",
        f"MIN_MEDIAN_FLUX threshold: {MIN_MEDIAN_FLUX:g} counts (masked median, every-10th-fiber subsample)",
        f"kept: apo={len(kept['apo'])} lco={len(kept['lco'])} total={len(kept['apo'])+len(kept['lco'])}",
        f"excluded: {len(excluded)} "
        f"(low_flux={sum(1 for _, r, _ in excluded if r.startswith('low_flux'))}, "
        f"bad_schema={sum(1 for _, r, _ in excluded if r == 'bad_schema')})",
        f"known junk lco_60855_0011 excluded: {junk_excluded}",
        "median-flux histogram (bins 0/1/10/100/1e3/1e4/1e5/inf): " + str(hist.tolist()),
    ]
    with open("list_stats.txt", "w") as out:
        out.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if not junk_excluded:
        print("ERROR: lco_60855_0011 was NOT excluded — check threshold!", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
