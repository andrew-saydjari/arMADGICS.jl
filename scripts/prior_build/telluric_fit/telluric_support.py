"""Data-derived pixel support + FIXED canonical Fourier basis for the dome-flat fit.

Why this module exists
----------------------
The previous fitter hardcoded

    CHIP_INDICES = [(340-20, 3430+20), (3675-20, 6240+20), (6430-20, 8490+20)]

which was telescope-blind and epoch-blind, and used it for BOTH jobs at once:

  (a) it selected which pixels got fit, and
  (b) it set the Fourier period / node placement of the continuum basis
      (``fourier_design_matrix_1d(ei - si, ...)`` -- the basis was defined ON the
      window, so the window's length and origin *were* the basis).

Two consequences.  First, the window was simply wrong for LCO: the LCO dome-flat
blue chip starts ~53 px bluer than APO's, so 44 px of 78-95%-live LCO data fell
outside the window, where the design matrix is identically zero and therefore
``Tfun = exp(A @ theta) = exp(0) = 1.0`` exactly -- a placeholder that reads
downstream as a real (tiny, after median normalisation) throughput.  Second,
because of (b), you cannot fix the window per telescope without silently
redefining what mode ``k`` means: a per-telescope window gives a per-telescope
period, and ``theta[.., k]`` stops being comparable across telescopes/fibers.

The fix is to separate the two jobs:

  * BASIS -- built once on a FIXED canonical interval per chip
    (``CANONICAL_CHIP_INTERVALS`` below).  Never data-derived, never per
    telescope, never per fiber.  Mode k's period is ``(2 N_canon + 1) / k`` px
    and its nodes sit at fixed absolute pixels, for every fiber, telescope and
    epoch.  ``theta`` coefficients are therefore directly comparable.
  * SUPPORT -- derived from the dome-flat stack, per telescope (optionally per
    fiber).  It selects which ROWS of the design matrix participate; it never
    touches the COLUMNS.  Widening or narrowing a support cannot move a node.

The support is then declared explicitly in the output file (``support`` dataset
+ provenance attrs) so a consumer can distinguish "no basis here" from "the
throughput really is this small" -- the missing declaration is what let the
exp(0) placeholder reach production.

Reference: 2026_09_07 FINDINGS.md (root-cause analysis of the LCO 67-px prior
hole) and 2026_07_19/2026_08_31 refactor plans.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

N_PIXELS = 8700

# Log-lambda grid shared by ar1D products and this fit (kept here so support
# bounds can be reported in Angstroms without importing the fitter).
RESAMPLED_WL = 10 ** (np.arange(N_PIXELS) * 6e-6 + 4.17825)

# ── THE FIXED CANONICAL BASIS INTERVALS ───────────────────────────────────────
# 0-based, half-open [start, stop) pixel intervals on the 8700-px grid, one per
# chip.  These define the Fourier basis and NOTHING else.
#
# FROZEN.  Changing any of these numbers redefines every mode and invalidates
# comparability of `theta` with every previously fit epoch, so it must never be
# done to "fit the data better" -- that is what the (row-selecting) support is
# for.  They were chosen (2026-09-07) as round numbers that
#   * contain every measured per-telescope, per-epoch, per-fiber dome-flat
#     footprint at both APO and LCO over MJD 57643-61160 with >= 20 px of slack
#     on each side, and
#   * stay disjoint, with the chip gaps comfortably inside the dead zones.
# `derive_support` ASSERTS containment, so an epoch that ever drifts outside
# these bounds fails loudly instead of being silently clipped.
#
# MEASURED extremes they had to contain (60 dome flats/telescope, MJD
# 57643-61160, over telescope x epoch-tertile x fiber; see SUPPORT_REPORT.md):
#   chip 0  bluest start 274 (LCO per-fiber)   reddest stop 3443 (APO per-fiber)
#   chip 1  bluest start 3635 (LCO)            reddest stop 6251 (APO per-fiber)
#   chip 2  bluest start 6399 (APO MJD>59200)  reddest stop 8520 (LCO)
CANONICAL_CHIP_INTERVALS = (
    (240, 3500),
    (3600, 6300),
    (6350, 8560),
)

SUPPORT_DERIVATION_VERSION = 1


# ── Canonical Fourier basis ───────────────────────────────────────────────────
def fourier_columns(pix, start, stop, n_modes):
    """Fourier columns evaluated at absolute pixel indices ``pix``.

    The basis is defined on the canonical interval ``[start, stop)`` of length
    ``N = stop - start``, with the same convention as
    ``domeflats.fourier_design_matrix_1d`` ("half the period"):

        t(p) = (p - start) / (2 N + 1)
        cols = [1, cos(2 pi k t), sin(2 pi k t)  for k = 1 .. n_modes]

    Evaluating at ``pix = arange(start, stop)`` reproduces
    ``fourier_design_matrix_1d(N, n_modes)`` exactly -- this is the same basis,
    just expressed as a function of ABSOLUTE pixel so that it can be evaluated
    on (or restricted to) any subset without changing period or phase.
    """
    N = stop - start
    t = (np.asarray(pix, dtype=np.float64) - start) / (2.0 * N + 1.0)
    cols = [np.ones_like(t)]
    for k in range(1, n_modes + 1):
        cols.append(np.cos(2.0 * np.pi * k * t))
        cols.append(np.sin(2.0 * np.pi * k * t))
    return np.column_stack(cols)


def n_fourier_per_chip(n_modes):
    return 2 * n_modes + 1


def n_fourier_total(n_modes, intervals=CANONICAL_CHIP_INTERVALS):
    return len(intervals) * n_fourier_per_chip(n_modes)


def canonical_design_matrix(n_modes, intervals=CANONICAL_CHIP_INTERVALS,
                            n_pixels=N_PIXELS):
    """(n_pixels, n_fourier) basis on the canonical intervals, before support.

    Rows outside every canonical interval are zero (those pixels are in a chip
    gap / off the detector and carry no continuum model by construction).
    """
    A = np.zeros((n_pixels, n_fourier_total(n_modes, intervals)), dtype=np.float64)
    npc = n_fourier_per_chip(n_modes)
    for i, (si, ei) in enumerate(intervals):
        A[si:ei, i * npc:(i + 1) * npc] = fourier_columns(np.arange(si, ei), si, ei, n_modes)
    return A


def canonical_prior_precision(n_modes, s_pixels, amplitude,
                              intervals=CANONICAL_CHIP_INTERVALS):
    """Matern prior precision vector over the CANONICAL modes.

    The length-scale is converted to mode units with the canonical chip length,
    not the (per-telescope) support length -- the prior lives on the canonical
    coefficients, so it must be the same everywhere for `theta` to be
    comparable.  Imported lazily to keep this module numpy-only for tests.
    """
    from domeflats import matern_prior_variance_1d

    lam = np.zeros(n_fourier_total(n_modes, intervals))
    npc = n_fourier_per_chip(n_modes)
    for i, (si, ei) in enumerate(intervals):
        lam[i * npc:(i + 1) * npc] = amplitude / np.array(
            matern_prior_variance_1d(n_modes, float(s_pixels) / (ei - si))
        )
    return lam


def mode_node_pixels(k, n_modes=None, intervals=CANONICAL_CHIP_INTERVALS):
    """Absolute pixel positions of the zero crossings of the cos-k column.

    Used by the mode-anchoring test: these depend only on the canonical
    interval, so they are identical for every telescope, fiber and support.
    """
    out = []
    for (si, ei) in intervals:
        N = ei - si
        # cos(2 pi k t) = 0  <=>  t = (2 m + 1) / (4 k)
        nodes = []
        m = 0
        while True:
            p = si + (2 * m + 1) / (4.0 * k) * (2.0 * N + 1.0)
            if p >= ei:
                break
            nodes.append(p)
            m += 1
        out.append(np.array(nodes))
    return out


# ── Data-derived support ──────────────────────────────────────────────────────
@dataclass
class SupportSpec:
    """Everything needed to rebuild the support, plus how it was derived."""
    telescope: str
    bounds: list                  # [[si, ei], ...] 0-based half-open, per chip
    canonical_intervals: list
    min_live_frac: float
    min_exposure_frac: float
    edge_buffer: int
    n_exposures_scanned: int
    n_exposures_in_list: int
    input_list: str
    input_list_sha256: str
    version: int = SUPPORT_DERIVATION_VERSION
    per_fiber: bool = False
    n_fiber: int = 0
    # per-fiber bounds, only when per_fiber: [[ [si,ei] x nchip ] x nfiber]
    fiber_bounds: list | None = None

    def support_mask(self, n_pixels=N_PIXELS):
        """bool (n_pixels,) or (n_pixels, n_fiber)."""
        if not self.per_fiber:
            m = np.zeros(n_pixels, dtype=bool)
            for si, ei in self.bounds:
                m[si:ei] = True
            return m
        m = np.zeros((n_pixels, self.n_fiber), dtype=bool)
        for f, fb in enumerate(self.fiber_bounds):
            for si, ei in fb:
                m[si:ei, f] = True
        return m

    def to_json(self, path):
        Path(path).write_text(json.dumps(asdict(self), indent=2) + "\n")

    @staticmethod
    def from_json(path):
        d = json.loads(Path(path).read_text())
        d["bounds"] = [tuple(b) for b in d["bounds"]]
        d["canonical_intervals"] = [tuple(b) for b in d["canonical_intervals"]]
        return SupportSpec(**d)


def sha256_of_lines(lines):
    h = hashlib.sha256()
    for line in lines:
        h.update(line.strip().encode())
        h.update(b"\n")
    return h.hexdigest()


def _runs(flag):
    """Contiguous True runs of a bool array as [start, stop) pairs."""
    idx = np.flatnonzero(flag)
    if idx.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(idx) != 1)
    starts = np.concatenate([[idx[0]], idx[breaks + 1]])
    stops = np.concatenate([idx[breaks], [idx[-1]]]) + 1
    return list(zip(starts.tolist(), stops.tolist()))


def accumulate_live(paths, min_live_frac=0.5, per_fiber=False, n_pixels=N_PIXELS,
                    progress=False):
    """Scan dome flats and accumulate liveness statistics.

    A pixel is "live" in one exposure/fiber when ``mask_1d != 0 & ivar_1d > 0``
    -- exactly the condition the fitter itself uses to build its weights, so the
    support cannot claim pixels the fit would have thrown away.

    Returns ``(exposure_frac, mean_fiber_frac, n_scanned, n_fiber)`` where

      exposure_frac    (n_pixels,)  fraction of exposures in which at least
                                    ``min_live_frac`` of fibers were live
                       (n_pixels, n_fiber) if per_fiber: fraction of exposures
                                    in which THAT fiber was live
      mean_fiber_frac  (n_pixels,)  mean over exposures of the live-fiber
                                    fraction (diagnostic only)
    """
    import h5py

    acc = None
    mean_ff = np.zeros(n_pixels)
    n = 0
    n_fiber = 0
    for j, p in enumerate(paths):
        if not Path(p).exists():
            continue
        try:
            with h5py.File(p, "r") as f:
                ivar = f["ivar_1d"][:]
                mask = f["mask_1d"][:]
        except Exception:
            continue
        if ivar.shape[0] == n_pixels and ivar.ndim == 2 and ivar.shape[1] != n_pixels:
            ivar, mask = ivar.T, mask.T
        good = (mask != 0) & (ivar > 0)          # (n_fiber, n_pixels)
        if n_fiber == 0:
            n_fiber = good.shape[0]
            acc = (np.zeros((n_pixels, n_fiber)) if per_fiber
                   else np.zeros(n_pixels))
        elif good.shape[0] != n_fiber:
            continue
        ff = good.mean(axis=0)
        mean_ff += ff
        if per_fiber:
            acc += good.T
        else:
            acc += (ff >= min_live_frac)
        n += 1
        if progress and (n % 10 == 0):
            print(f"  scanned {n} ({j + 1}/{len(paths)})", flush=True)
    if n == 0:
        raise RuntimeError("no readable dome flats in the given list")
    return acc / n, mean_ff / n, n, n_fiber


# Default inward trim at each chip edge.  The OLD code padded OUTWARD by 20 px,
# which is what pushed the fit window past the real data at every APO edge.  The
# sign is the point; the magnitude should be as small as the data allows, since
# AKS's directive is to "preserve as much wavelength range as possible".
# The (min_live_frac,
# min_exposure_frac) thresholds are what actually decide whether a pixel has
# usable data; this buffer only absorbs the +/-1 px quantisation of the
# threshold crossing.  MEASURED cost of larger values, on the LCO red edge that
# is one of the three runs of dead prior pixels this fix exists to recover:
# the live-fiber fraction there is 0.90 / 0.88 / 0.66 / 0.33 at px
# 8509 / 8511 / 8513 / 8515, so edge_buffer=5 would throw away 8510-8512 at
# ~88% live and recover NONE of that run; edge_buffer=20 would cost 12 x 20 =
# 240 px of grid, 3.5x the 67 px the whole fix recovers.  See SUPPORT_REPORT.md
# for the full sensitivity table -- one flag changes it.
DEFAULT_EDGE_BUFFER = 2


def bounds_from_live(exposure_frac, min_exposure_frac=0.5,
                     edge_buffer=DEFAULT_EDGE_BUFFER,
                     intervals=CANONICAL_CHIP_INTERVALS, min_run=100,
                     clip_to_canonical=False):
    """Turn a per-pixel liveness curve into one [start, stop) per chip.

    * threshold at ``min_exposure_frac``;
    * keep contiguous runs longer than ``min_run`` px (drops isolated hot
      pixels in the gaps and single dead columns are bridged by the run, not
      by an interior hole -- support is contiguous per chip by construction);
    * ERODE each run by ``edge_buffer`` px on both sides.  The old code DILATED
      by 20 px, which is what pushed the fit window past the real data at every
      APO edge; the sign of this margin is the point.
    * assert each surviving run lands inside exactly one canonical interval.
    """
    alive = exposure_frac >= min_exposure_frac
    raw = [(a, b) for a, b in _runs(alive) if (b - a) >= min_run]
    if len(raw) != len(intervals):
        raise RuntimeError(
            f"expected {len(intervals)} chip runs, found {len(raw)}: {raw}"
        )
    out = []
    for (a, b), (ci, ce) in zip(raw, intervals):
        si, ei = a + edge_buffer, b - edge_buffer
        if ei - si < min_run:
            raise RuntimeError(f"chip run {(a, b)} vanishes under erosion {edge_buffer}")
        if not (ci <= si and ei <= ce):
            if clip_to_canonical:
                si, ei = max(si, ci), min(ei, ce)
            else:
                raise RuntimeError(
                    f"derived support [{si},{ei}) escapes canonical interval "
                    f"[{ci},{ce}) -- the canonical basis is FROZEN; investigate "
                    f"the data before touching CANONICAL_CHIP_INTERVALS"
                )
        out.append((int(si), int(ei)))
    return out, raw


def derive_support(list_path, telescope, min_live_frac=0.5, min_exposure_frac=0.5,
                   edge_buffer=DEFAULT_EDGE_BUFFER, n_sample=None, per_fiber=False,
                   intervals=CANONICAL_CHIP_INTERVALS, progress=False):
    """Derive a :class:`SupportSpec` from a dome-flat input list.

    ``n_sample`` evenly subsamples the (mjd, exposure)-sorted list; ``None``
    scans everything.  60 exposures spanning the full MJD range already pin
    every boundary to +/-1 px; the parameter is exposed so the full list can be
    scanned for the record.

    MEASURED caveat (SUPPORT_REPORT.md section 4): the APO footprint shifts
    ~30 px blueward at MJD ~59200 and stays there, so a single global support
    is the *median exposure's* footprint, not every exposure's.  At
    ``min_exposure_frac=0.5`` that is the older epoch; 0.35 takes the union
    instead.  A per-epoch support cannot be expressed by the current
    one-design_matrix-per-product layout.
    """
    lines = [l.strip() for l in Path(list_path).read_text().splitlines() if l.strip()]
    n_all = len(lines)
    if n_sample is not None and n_sample < n_all:
        sel = np.unique(np.round(np.linspace(0, n_all - 1, n_sample)).astype(int))
        paths = [lines[i] for i in sel]
    else:
        paths = lines

    exp_frac, mean_ff, n_scanned, n_fiber = accumulate_live(
        paths, min_live_frac=min_live_frac, per_fiber=per_fiber, progress=progress
    )

    if per_fiber:
        # global curve = fraction of FIBERS whose own exposure-liveness clears
        # min_exposure_frac; thresholded at min_live_frac.  Same two thresholds
        # as the global path, applied in the other order.
        global_live = (exp_frac >= min_exposure_frac).mean(axis=1)
        global_bounds, _ = bounds_from_live(
            global_live, min_exposure_frac=min_live_frac, edge_buffer=edge_buffer,
            intervals=intervals,
        )
        fiber_bounds, n_fallback = [], 0
        for f in range(n_fiber):
            try:
                fb, _ = bounds_from_live(exp_frac[:, f], min_exposure_frac,
                                         edge_buffer, intervals)
            except RuntimeError:
                # Dead or degenerate fiber (extra runs from an interior gap, a
                # run escaping canonical, ...).  Fall back to the global bounds;
                # the per-exposure ivar mask still excludes its dead pixels.
                fb = list(global_bounds)
                n_fallback += 1
            fiber_bounds.append([list(x) for x in fb])
        if n_fallback:
            print(f"  per-fiber: {n_fallback}/{n_fiber} fibers fell back to the "
                  f"global bounds (no clean 3-chip run of their own)")
        bounds = global_bounds
    else:
        bounds, _ = bounds_from_live(exp_frac, min_exposure_frac, edge_buffer, intervals)
        fiber_bounds = None

    return SupportSpec(
        telescope=telescope,
        bounds=[list(b) for b in bounds],
        canonical_intervals=[list(b) for b in intervals],
        min_live_frac=float(min_live_frac),
        min_exposure_frac=float(min_exposure_frac),
        edge_buffer=int(edge_buffer),
        n_exposures_scanned=int(n_scanned),
        n_exposures_in_list=int(n_all),
        input_list=str(Path(list_path).resolve()),
        input_list_sha256=sha256_of_lines(lines),
        per_fiber=bool(per_fiber),
        n_fiber=int(n_fiber),
        fiber_bounds=fiber_bounds,
    ), exp_frac, mean_ff


def design_matrix_for_support(n_modes, support, intervals=CANONICAL_CHIP_INTERVALS,
                              n_pixels=N_PIXELS):
    """Canonical basis with rows outside the (global) support zeroed.

    Rows are SELECTED, columns are untouched -- this is the whole point.  For a
    per-fiber support the emitted design matrix uses the union over fibers and
    the per-fiber mask is carried separately in ``support``; per-fiber row
    selection is applied through the weights (see fit_domeflats.py).
    """
    A = canonical_design_matrix(n_modes, intervals, n_pixels)
    sup = np.asarray(support)
    if sup.ndim == 2:
        sup = sup.any(axis=1)
    A[~sup, :] = 0.0
    return A


def check_support_matches_design(support, A):
    """G1: the declared support IS the nonzero-row support of the matrix."""
    sup = np.asarray(support)
    if sup.ndim == 2:
        sup = sup.any(axis=1)
    nz = np.any(A != 0.0, axis=1)
    return bool(np.array_equal(sup, nz)), sup, nz


# ── CLI ───────────────────────────────────────────────────────────────────────
def _main():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--list", required=True, help="dome-flat input list (one path/line)")
    ap.add_argument("--telescope", required=True, choices=["apo", "lco"])
    ap.add_argument("--out", required=True, help="output support JSON")
    ap.add_argument("--n-sample", type=int, default=60,
                    help="evenly subsample the list (0 = scan all)")
    ap.add_argument("--min-live-frac", type=float, default=0.5)
    ap.add_argument("--min-exposure-frac", type=float, default=0.5)
    ap.add_argument("--edge-buffer", type=int, default=DEFAULT_EDGE_BUFFER)
    ap.add_argument("--per-fiber", action="store_true")
    ap.add_argument("--npz", default=None, help="also dump the liveness curves")
    args = ap.parse_args()

    spec, exp_frac, mean_ff = derive_support(
        args.list, args.telescope,
        min_live_frac=args.min_live_frac,
        min_exposure_frac=args.min_exposure_frac,
        edge_buffer=args.edge_buffer,
        n_sample=(None if args.n_sample == 0 else args.n_sample),
        per_fiber=args.per_fiber,
        progress=True,
    )
    spec.to_json(args.out)
    if args.npz:
        np.savez_compressed(args.npz, exposure_frac=exp_frac, mean_fiber_frac=mean_ff)

    print(f"\n{args.telescope}: scanned {spec.n_exposures_scanned} of "
          f"{spec.n_exposures_in_list} exposures")
    for i, (si, ei) in enumerate(spec.bounds):
        ci, ce = spec.canonical_intervals[i]
        print(f"  chip {i}: support [{si}, {ei})  ({ei - si} px, "
              f"{RESAMPLED_WL[si]:.1f}-{RESAMPLED_WL[ei - 1]:.1f} A)"
              f"   canonical [{ci}, {ce})  slack {si - ci}/{ce - ei} px")
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    _main()
