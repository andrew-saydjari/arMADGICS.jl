"""Merge per-shard E3 refit outputs into per-telescope product files + verify.

Run AFTER all status/shard_<k>.STATUS say state=DONE:

    nice -n 10 uv run python merge_shards.py

Produces (in out/):
    tellurics_refit_20260902_apo.h5     delivered-style layout (see 20260323.py
    tellurics_refit_20260902_lco.h5     docstring) + chi_sq_fiber, concatenated
                                        along the file axis in full-list order
    MERGE_REPORT.txt                    counts, T-update spot check, Tfunsample
                                        population quantiles vs delivered refs
    ALL_DONE                            sentinel, written only if every check passes

Screening note: chi_sq_fiber is RECORDED here, not applied — the per-fiber chi2
screen (top ~0.1%) is applied at consumption when drawing prior training samples.
"""
import argparse
import pickle
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import h5py

_ap = argparse.ArgumentParser()
_ap.add_argument("--base", default=None, help="run directory (default: this file's)")
_ap.add_argument("--date-tag", default="20260902")
_args = _ap.parse_args()

BASE = Path(_args.base).resolve() if _args.base else Path(__file__).resolve().parent
OUT = BASE / "out"
DATE_TAG = _args.date_tag
STACK_KEYS = ["theta", "gamma", "chi_sq", "median_resid", "chi_sq_fiber",
              "stage", "paths", "success", "T_final"]
# Per-file (not per-shard) datasets, copied once from the first shard and then
# asserted identical across every shard of the telescope.
SHARED_KEYS = ["design_matrix", "support"]
QUANTS = [1, 5, 25, 50, 75, 95, 99]
# E3_REPORT.md delivered-production reference quantiles (chip-pixel Tfunsample).
# NOTE: these are PRE-support-fix references.  A corrected support changes the
# canonical basis coverage and therefore theta everywhere, so these quantiles
# are EXPECTED to move; they are printed for comparison, not gated on.
REF = {
    "apo": [0.298, 0.414, 0.720, 1.039, 1.196, 1.408, 1.730],
    "lco": [0.343, 0.465, 0.786, 1.036, 1.150, 1.269, 1.292],
}
N_SAMPLE_FILES = 200      # per telescope, for quantile + T-change spot checks
rng = np.random.default_rng(20260902)

# ── 1. read shard map, require all shards DONE ────────────────────────────────
shards = []   # (sid, tele, n_files, list_path, out_path)
for line in (BASE / "shards" / "shard_map.tsv").read_text().splitlines():
    sid, tele, n, lp, op = line.split("\t")
    shards.append((int(sid), tele, int(n), lp, op))

not_done = []
for sid, tele, n, lp, op in shards:
    st = BASE / "status" / f"shard_{sid}.STATUS"
    state = ""
    if st.exists():
        for l in st.read_text().splitlines():
            if l.startswith("state="):
                state = l.split("=", 1)[1]
    if state != "DONE":
        not_done.append((sid, state or "NO_STATUS"))
if not_done:
    print(f"REFUSING to merge: shards not DONE: {not_done}", file=sys.stderr)
    sys.exit(1)

with open(BASE / "T_init.pkl", "rb") as fp:
    T_init = np.asarray(pickle.load(fp), dtype=np.float64)


report = [f"E3 full refit merge — {datetime.now().isoformat()}", ""]
ok_all = True

# ── 2. merge per telescope ────────────────────────────────────────────────────
for tele in ("apo", "lco"):
    tele_shards = [s for s in shards if s[1] == tele]
    tele_shards.sort(key=lambda s: s[0])
    full_list = [l for l in (BASE / f"list_{tele}_full.txt").read_text().splitlines() if l.strip()]

    merged_path = OUT / f"tellurics_refit_{DATE_TAG}_{tele}.h5"
    if merged_path.exists():
        print(f"{merged_path} exists — remove it to re-merge", file=sys.stderr)
        sys.exit(1)

    with h5py.File(merged_path, "w") as fo:
        offset = 0
        attrs_done = False
        for sid, _, n, lp, op in tele_shards:
            with h5py.File(op, "r") as fi:
                ns = fi["success"].shape[0]
                assert ns == n, f"shard {sid}: {ns} rows != map {n}"
                if "support" not in fi:
                    print(f"shard {sid}: no `support` dataset -- refit with "
                          "fit_domeflats.py before merging", file=sys.stderr)
                    sys.exit(1)
                if not attrs_done:
                    total = sum(s[2] for s in tele_shards)
                    for k in SHARED_KEYS:
                        fo.create_dataset(k, data=fi[k][:])
                    for k in STACK_KEYS:
                        src = fi[k]
                        fo.create_dataset(k, (total,) + src.shape[1:], dtype=src.dtype)
                    for a, v in fi.attrs.items():
                        fo.attrs[a] = v
                    fo.attrs["merged_from_shards"] = [s[0] for s in tele_shards]
                    fo.attrs["merge_date"] = DATE_TAG
                    attrs_done = True
                else:
                    # A shard fit against a different support would make theta
                    # mean something different -- refuse rather than merge.
                    for k in SHARED_KEYS:
                        if not np.array_equal(fo[k][:], fi[k][:]):
                            print(f"shard {sid}: `{k}` differs from shard "
                                  f"{tele_shards[0][0]}; refusing to merge",
                                  file=sys.stderr)
                            sys.exit(1)
                for k in STACK_KEYS:
                    fo[k][offset:offset + n] = fi[k][:]
            offset += n

        succ = fo["success"][:]
        paths = [p.decode() for p in fo["paths"][:]]
        n_total, n_succ = len(succ), int(succ.sum())
        list_match = paths == full_list

        report += [f"[{tele}] merged {len(tele_shards)} shards -> {merged_path.name}",
                   f"[{tele}] rows={n_total} success={n_succ} "
                   f"list_match={'OK' if list_match else 'MISMATCH'}"]
        if n_succ != n_total or not list_match:
            ok_all = False

        # spot checks on a random subsample
        idx = rng.choice(n_total, size=min(N_SAMPLE_FILES, n_total), replace=False)
        idx.sort()
        A = fo["design_matrix"][:].astype(np.float64)          # (N_PIXELS, n_fourier)
        sup_full = fo["support"][:].astype(bool)
        chip_mask = sup_full.any(axis=1) if sup_full.ndim == 2 else sup_full
        # G1 on the merged product, at write time.
        nz = np.any(A != 0.0, axis=1)
        g1 = np.array_equal(nz, chip_mask)
        report += [f"[{tele}] G1 support == nonzero design_matrix rows: "
                   f"{'OK' if g1 else '** FAIL **'} "
                   f"({int(chip_mask.sum())} px, bounds "
                   f"{np.asarray(fo.attrs['support_bounds']).tolist()})"]
        if not g1:
            ok_all = False
        max_dT, tfun_chip = [], []
        n_canary = 0
        for i in idx:
            Tf = fo["T_final"][i].astype(np.float64)
            max_dT.append(np.abs(Tf - T_init)[:, chip_mask].max())
            th = fo["theta"][i].astype(np.float64)             # (N_FIBERS, n_fourier)
            tf_all = np.exp(th @ A.T)
            tf_all[:, ~chip_mask] = 0.0        # the consumer contract
            tf = tf_all[:, chip_mask]
            med = np.median(tf, axis=1, keepdims=True)
            good = np.isfinite(med[:, 0]) & (med[:, 0] > 0)
            tfun_chip.append((tf[good] / med[good]).ravel())
            # G4: after the contract, no pixel may equal 1/median exactly
            n_canary += int(np.sum(tf_all[good] / med[good] == 1.0 / med[good]))
        max_dT = np.array(max_dT)
        n_unchanged = int((max_dT < 1e-8).sum())
        q = np.percentile(np.concatenate(tfun_chip), QUANTS)
        report += [
            f"[{tele}] T-update spot check ({len(idx)} files): "
            f"max|T_final-T_init| median={np.median(max_dT):.3f} "
            f"min={max_dT.min():.4f}; files with T unchanged: {n_unchanged} "
            f"{'** BUG SIGNATURE — INVESTIGATE **' if n_unchanged else '(OK)'}",
            f"[{tele}] Tfunsample quantiles {QUANTS}:",
            f"[{tele}]   merged:    " + "/".join(f"{v:.3f}" for v in q),
            f"[{tele}]   delivered: " + "/".join(f"{v:.3f}" for v in REF[tele]),
            f"[{tele}] G4 exp(0) placeholder canary (support applied in linear "
            f"space): {n_canary} hits {'(OK)' if n_canary == 0 else '** FAIL **'}",
            f"[{tele}] chi_sq_fiber p99.9 = "
            f"{np.nanpercentile(fo['chi_sq_fiber'][:], 99.9):.2f} "
            f"(screen threshold candidate; apply at consumption)",
            "",
        ]
        if n_unchanged or n_canary:
            ok_all = False

# ── 3. report + sentinel ──────────────────────────────────────────────────────
report.append("ALL CHECKS PASSED" if ok_all else "CHECKS FAILED — see above; ALL_DONE not written")
(OUT / "MERGE_REPORT.txt").write_text("\n".join(report) + "\n")
print("\n".join(report))
if ok_all:
    (OUT / "ALL_DONE").write_text(datetime.now().isoformat() + "\n")
    print(f"\nwrote {OUT}/ALL_DONE")
else:
    sys.exit(1)
