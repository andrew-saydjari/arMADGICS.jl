"""Run the acceptance gates against a PRODUCED dome-flat fit file.

    python validate_output.py out/tellurics_refit_lco.h5 \
        [--live live_lco.npz] [--edge-buffer 2] [--n-theta 20]

G1  support == nonzero-row support of design_matrix, and support declared at all
G2  all six chip boundaries inside the measured live region by <= edge_buffer px
    (skipped without --live, the liveness npz written by telluric_support.py)
G4  no pixel of a normalised Tfun sample equals 1/median(Tfun), once the
    declared support has been applied in linear space -- and the canary is shown
    to FIRE for a consumer that ignores the declaration, so the gate has teeth

Exit status 0 iff every runnable gate passes.
"""
import argparse
import sys

import numpy as np
import h5py


def _runs(flag):
    idx = np.flatnonzero(flag)
    if idx.size == 0:
        return []
    br = np.flatnonzero(np.diff(idx) != 1)
    st = np.concatenate([[idx[0]], idx[br + 1]])
    sp = np.concatenate([idx[br], [idx[-1]]]) + 1
    return [(int(a), int(b)) for a, b in zip(st, sp)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--live", default=None,
                    help="npz from telluric_support.py --npz (exposure_frac)")
    ap.add_argument("--edge-buffer", type=int, default=None,
                    help="default: the value recorded in the file's attrs")
    ap.add_argument("--min-exposure-frac", type=float, default=None)
    ap.add_argument("--n-theta", type=int, default=20,
                    help="how many (file, fiber) theta draws to run G4 on")
    args = ap.parse_args()

    fails, notes = [], []
    with h5py.File(args.path, "r") as f:
        A = f["design_matrix"][:].astype(np.float64)
        legacy = "support" not in f
        if legacy:
            # Legacy product (pre-declaration).  Grade it anyway, against the
            # support it *implies*, so it can serve as the regression baseline.
            print("G1 FAIL: no `support` dataset -- this file predates the "
                  "support declaration, so a consumer cannot tell 'no basis "
                  "here' from a real throughput.  Grading G2/G4 against the "
                  "implied (nonzero-row) support as a regression baseline.\n")
            fails.append("G1")
            sup_full = np.any(A != 0.0, axis=1)
        else:
            sup_full = f["support"][:].astype(bool)
        attrs = dict(f.attrs)
        theta = f["theta"]
        succ = f["success"][:]
        n_ok = int(succ.sum())
        rng = np.random.default_rng(20260907)
        idx_file = rng.choice(np.flatnonzero(succ), size=min(args.n_theta, n_ok),
                              replace=False) if n_ok else np.array([], dtype=int)
        idx_fib = rng.integers(0, theta.shape[1], size=idx_file.size)
        thetas = [theta[int(i), int(j)].astype(np.float64)
                  for i, j in zip(np.sort(idx_file), idx_fib)]

    sup = sup_full.any(axis=1) if sup_full.ndim == 2 else sup_full
    fill = attrs.get("off_support_fill", "zero")
    eb = args.edge_buffer if args.edge_buffer is not None \
        else int(attrs.get("support_edge_buffer", 2))

    print(f"file      : {args.path}")
    print(f"telescope : {attrs.get('support_telescope', '?')}   "
          f"per_fiber={bool(attrs.get('support_per_fiber', False))}   fill={fill}")
    print(f"canonical : {np.asarray(attrs['canonical_chip_intervals']).tolist() if 'canonical_chip_intervals' in attrs else '(not declared)'}")
    print(f"support   : {np.asarray(attrs['support_bounds']).tolist() if 'support_bounds' in attrs else _runs(sup)}  "
          f"({int(sup.sum())} px){'  [IMPLIED, not declared]' if legacy else ''}")
    print(f"provenance: thresholds {attrs.get('support_min_live_frac')}/"
          f"{attrs.get('support_min_exposure_frac')}  edge_buffer={eb}  "
          f"list_sha256={str(attrs.get('support_input_list_sha256'))[:16]}...  "
          f"n_scanned={attrs.get('support_n_exposures_scanned')}")
    print(f"fits      : {n_ok}/{len(succ)} successful\n")

    # ── G1 ────────────────────────────────────────────────────────────────────
    nz = ~np.all(np.isnan(A), axis=1) if fill == "nan" else np.any(A != 0.0, axis=1)
    if legacy:
        pass                                   # already reported above
    elif np.array_equal(nz, sup):
        print(f"G1 PASS  support == nonzero-row support of design_matrix "
              f"({int(sup.sum())} px, runs {_runs(sup)})")
    else:
        fails.append("G1")
        print(f"G1 FAIL  {int(sup.sum())} declared vs {int(nz.sum())} nonzero rows")

    # ── G2 ────────────────────────────────────────────────────────────────────
    if args.live is None:
        notes.append("G2 SKIP (no --live)")
        print("G2 SKIP  pass --live <npz from telluric_support.py --npz>")
    else:
        z = np.load(args.live)
        E_raw = z["exposure_frac"]
        thr = args.min_exposure_frac if args.min_exposure_frac is not None \
            else float(attrs.get("support_min_exposure_frac", 0.5))
        per_fiber = sup_full.ndim == 2

        def grade(support_1d, curve, label):
            live = [(a, b) for a, b in _runs(curve >= thr) if b - a >= 100]
            bounds = _runs(support_1d)
            bad = []
            if len(live) != len(bounds):
                bad.append(("chip count", len(bounds), len(live)))
            else:
                for c, ((s, e), (a, b)) in enumerate(zip(bounds, live)):
                    if not (a <= s <= a + eb):
                        bad.append((c, "start", s, (a, b)))
                    if not (b - eb <= e <= b):
                        bad.append((c, "stop", e, (a, b)))
            return bad, live, bounds

        if per_fiber and E_raw.ndim == 2:
            # Correct reference for a per-fiber support: each fiber's own curve.
            n_bad_fib, worst = 0, None
            for f in range(sup_full.shape[1]):
                bad, live, _ = grade(sup_full[:, f], E_raw[:, f], f"fiber {f}")
                if bad:
                    n_bad_fib += 1
                    worst = worst or (f, bad)
            if n_bad_fib:
                fails.append("G2")
                print(f"G2 FAIL  {n_bad_fib}/{sup_full.shape[1]} fibers have a "
                      f"boundary outside their own live region; first: {worst}")
            else:
                print(f"G2 PASS  all {sup_full.shape[1]} per-fiber supports "
                      f"inside their own live region by <= {eb} px")
            live = []
        else:
            if per_fiber:
                print("G2 NOTE  the declared support is PER FIBER but --live "
                      "carries only the global liveness curve.  Grading the "
                      "per-fiber UNION against the global curve, which is the "
                      "wrong reference: a union legitimately reaches 1-2 px "
                      "past the 50%-of-fibers boundary.  Re-derive with "
                      "--per-fiber --npz for a correct grade.")
            E = E_raw
            if E.ndim == 2:
                E = (E >= thr).mean(axis=1)
            bad, live, bounds = grade(sup, E, "global")
            n_edges = 2 * len(live)
            if bad:
                fails.append("G2")
                print(f"G2 FAIL  {len(bad)}/{n_edges} boundaries: {bad}")
            else:
                print(f"G2 PASS  {n_edges}/{n_edges} boundaries inside the live "
                      f"region by <= {eb} px   live={live} support={bounds}")
        off = int((sup & ~np.isin(np.arange(len(sup)),
                                  np.concatenate([np.arange(a, b) for a, b in live])
                                  if live else np.array([], int))).sum())
        if live:
            print(f"         support pixels with no measured data: {off}")

    # ── G4 ────────────────────────────────────────────────────────────────────
    if not thetas:
        notes.append("G4 SKIP (no successful fits)")
        print("G4 SKIP  no successful fits in this file")
    else:
        n_ignoring, n_contract = 0, 0
        for th in thetas:
            tf = np.exp(A @ th)
            vals = tf[np.isfinite(tf) & (tf != 0)]
            med = np.median(vals)
            n_ignoring += int(np.sum(np.isclose(tf / med, 1.0 / med, rtol=0, atol=0)))
            tf2 = tf.copy()
            tf2[~sup] = 0.0                       # the consumer contract
            vals2 = tf2[np.isfinite(tf2) & (tf2 != 0)]
            med2 = np.median(vals2)
            n_contract += int(np.sum(np.isclose(tf2 / med2, 1.0 / med2, rtol=0, atol=0)))
        if n_contract:
            fails.append("G4")
            print(f"G4 FAIL  {n_contract} placeholder pixels survive the support "
                  f"contract over {len(thetas)} draws")
        else:
            print(f"G4 PASS  0 placeholder pixels over {len(thetas)} Tfun draws "
                  f"after applying `support` in linear space")
        exp_hits = 0 if fill == "nan" else int((~sup).sum()) * len(thetas)
        verdict = "canary has teeth" if n_ignoring == exp_hits else \
            f"UNEXPECTED (got {n_ignoring}, expected {exp_hits})"
        print(f"         consumer IGNORING `support`: {n_ignoring} placeholder "
              f"pixels ({exp_hits} expected for fill={fill}) -- {verdict}")

    print()
    for n in notes:
        print(n)
    if fails:
        print(f"FAILED: {', '.join(fails)}")
        return 1
    print("ALL RUNNABLE GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
