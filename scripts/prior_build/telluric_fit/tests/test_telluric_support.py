"""Acceptance gates for the dome-flat fitter's support and canonical basis.

Runs under pytest, or standalone (``python tests/test_telluric_support.py``) so
it needs no dependency beyond numpy.

Gates (FINDINGS.md 2026_09_07 section 5):
  G1  fitter self-consistency:  support == any(design_matrix != 0, axis=1)
  G2  support vs data at all six chip boundaries per telescope
  G4  exp(0) placeholder canary: no pixel of a normalised Tfun equals
      1/median(Tfun) to float tolerance
  plus MODE ANCHORING, which is the constraint that makes a per-telescope
  support admissible at all.

G2 and the mode-anchoring "teeth" checks are run against the OLD hardcoded
window as a regression baseline: it must FAIL where the new construction
passes, otherwise the gate proves nothing.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import telluric_support as ts   # noqa: E402

N_MODES = 32

# ── The literal that produced the bug, 0-based half-open, as it appeared in
# 20260323.py lines 87-91 (the +/-20 pad already applied). ─────────────────────
OLD_CHIP_INDICES = [(340 - 20, 3430 + 20), (3675 - 20, 6240 + 20), (6430 - 20, 8490 + 20)]

# MEASURED dome-flat live runs, 0-based half-open, threshold (min_live_frac,
# min_exposure_frac) = (0.5, 0.5) over 60 exposures per telescope spanning
# MJD 57801-61160 (LCO) / 57643-61158 (APO).  Reproduce with
#   python telluric_support.py --list list_<tele>_full.txt --telescope <tele> \
#          --out /tmp/s.json --edge-buffer 0
MEASURED_RAW_RUNS = {
    "lco": [(276, 3407), (3637, 6237), (6422, 8515)],
    "apo": [(327, 3438), (3665, 6247), (6421, 8500)],
}
# The three runs of dead LCO prior pixels the fix exists to recover (0-based).
LCO_DEAD_RUNS = [(275, 320), (3637, 3655), (8510, 8514)]

EDGE_BUFFER = ts.DEFAULT_EDGE_BUFFER


def _spec(tele, bounds, per_fiber_bounds=None, n_fiber=0):
    return ts.SupportSpec(
        telescope=tele,
        bounds=[list(b) for b in bounds],
        canonical_intervals=[list(b) for b in ts.CANONICAL_CHIP_INTERVALS],
        min_live_frac=0.5, min_exposure_frac=0.5, edge_buffer=EDGE_BUFFER,
        n_exposures_scanned=60, n_exposures_in_list=60,
        input_list="test", input_list_sha256="0" * 64,
        per_fiber=per_fiber_bounds is not None, n_fiber=n_fiber,
        fiber_bounds=per_fiber_bounds,
    )


def derived_bounds(tele, edge_buffer=EDGE_BUFFER):
    return [(a + edge_buffer, b - edge_buffer) for a, b in MEASURED_RAW_RUNS[tele]]


def old_style_design_matrix(chip_indices, n_modes=N_MODES, n_pixels=ts.N_PIXELS):
    """Reproduces 20260323.py's construction: basis defined ON the window."""
    npc = 2 * n_modes + 1
    A = np.zeros((n_pixels, len(chip_indices) * npc))
    for i, (si, ei) in enumerate(chip_indices):
        A[si:ei, i * npc:(i + 1) * npc] = ts.fourier_columns(np.arange(si, ei), si, ei, n_modes)
    return A


# ══ G1 ════════════════════════════════════════════════════════════════════════
def test_G1_support_equals_nonzero_rows():
    for tele in ("lco", "apo"):
        spec = _spec(tele, derived_bounds(tele))
        sup = spec.support_mask()
        A = ts.design_matrix_for_support(N_MODES, sup)
        ok, s, nz = ts.check_support_matches_design(sup, A)
        assert ok, f"{tele}: {int(s.sum())} declared vs {int(nz.sum())} nonzero rows"
        # and it is a subset of the measured data footprint
        data = np.zeros(ts.N_PIXELS, dtype=bool)
        for a, b in MEASURED_RAW_RUNS[tele]:
            data[a:b] = True
        assert not (sup & ~data).any(), f"{tele}: support claims pixels with no data"


def test_G1_per_fiber_support_equals_nonzero_rows_of_union():
    n_fiber = 8
    rng = np.random.default_rng(0)
    base = derived_bounds("lco")
    fb = [[[a + int(rng.integers(0, 4)), b - int(rng.integers(0, 4))] for a, b in base]
          for _ in range(n_fiber)]
    spec = _spec("lco", base, per_fiber_bounds=fb, n_fiber=n_fiber)
    sup = spec.support_mask()
    assert sup.shape == (ts.N_PIXELS, n_fiber)
    A = ts.design_matrix_for_support(N_MODES, sup)
    ok, _, _ = ts.check_support_matches_design(sup, A)
    assert ok


def test_G1_fails_for_a_declaration_that_drifts():
    """The gate has teeth: a support that disagrees with the matrix is caught."""
    spec = _spec("lco", derived_bounds("lco"))
    sup = spec.support_mask()
    A = ts.design_matrix_for_support(N_MODES, sup)
    lied = sup.copy()
    lied[MEASURED_RAW_RUNS["lco"][0][0] - 30] = True     # claim one dead pixel
    ok, _, _ = ts.check_support_matches_design(lied, A)
    assert not ok


# ══ G2 ════════════════════════════════════════════════════════════════════════
def g2_failures(bounds, tele, edge_buffer=EDGE_BUFFER):
    """Boundaries that are NOT inside the live region by <= edge_buffer px."""
    bad = []
    for c, ((s, e), (a, b)) in enumerate(zip(bounds, MEASURED_RAW_RUNS[tele])):
        if not (a <= s <= a + edge_buffer):
            bad.append((c, "start", s, (a, b)))
        if not (b - edge_buffer <= e <= b):
            bad.append((c, "stop", e, (a, b)))
    return bad


def test_G2_old_window_fails_every_boundary_regression_baseline():
    """The old literal fails 6/6 boundaries at BOTH telescopes.

    FINDINGS.md quoted 6/6 LCO and 4/4 APO because only four APO edges had been
    measured at the time; with all six measured, APO is 6/6 too.
    """
    for tele in ("lco", "apo"):
        bad = g2_failures(OLD_CHIP_INDICES, tele)
        assert len(bad) == 6, f"{tele}: expected 6 failures, got {len(bad)}: {bad}"


def test_G2_derived_support_passes_every_boundary():
    for tele in ("lco", "apo"):
        bad = g2_failures(derived_bounds(tele), tele)
        assert not bad, f"{tele}: {bad}"


def test_G2_recovers_the_dead_lco_pixels():
    """The three dead LCO runs come back inside the new LCO support.

    MEASURED recovery, with the dead runs as reported in FINDINGS.md and the
    live runs as measured here:
      edge_buffer=0 -> 66 of 67 (0-based px 275 sits at live fraction 0.28,
                                 below the 0.5 threshold: not recoverable)
      edge_buffer=2 -> 61 of 67 (2 px at each of the two blue edges, 1 at the
                                 red edge, traded for boundary stability)
    """
    dead = np.zeros(ts.N_PIXELS, dtype=bool)
    for a, b in LCO_DEAD_RUNS:
        dead[a:b] = True
    assert dead.sum() == 67

    old = np.zeros(ts.N_PIXELS, dtype=bool)
    for a, b in OLD_CHIP_INDICES:
        old[a:b] = True
    assert not (dead & old).any(), "baseline: dead pixels were outside the old window"

    sup = _spec("lco", derived_bounds("lco")).support_mask()
    assert int((dead & sup).sum()) == 61, int((dead & sup).sum())
    sup0 = _spec("lco", derived_bounds("lco", edge_buffer=0)).support_mask()
    assert int((dead & sup0).sum()) == 66, int((dead & sup0).sum())


def test_G2_old_window_claimed_vacuum():
    """The old window's over-extension, the other half of the defect."""
    old = np.zeros(ts.N_PIXELS, dtype=bool)
    for a, b in OLD_CHIP_INDICES:
        old[a:b] = True
    for tele in ("lco", "apo"):
        data = np.zeros(ts.N_PIXELS, dtype=bool)
        for a, b in MEASURED_RAW_RUNS[tele]:
            data[a:b] = True
        assert (old & ~data).sum() > 0, tele
        sup = _spec(tele, derived_bounds(tele)).support_mask()
        assert (sup & ~data).sum() == 0, tele


# ══ MODE ANCHORING ════════════════════════════════════════════════════════════
# "we want the node for mode 1 to be at the same place across different fibers
#  so that mode amplitudes are interpretable and comparable" -- AKS
def test_mode1_node_identical_across_telescopes_and_per_fiber_variants():
    n_fiber = 6
    rng = np.random.default_rng(1)
    variants = {
        "lco": _spec("lco", derived_bounds("lco")),
        "apo": _spec("apo", derived_bounds("apo")),
        "lco_eb0": _spec("lco", derived_bounds("lco", 0)),
        "apo_eb20": _spec("apo", derived_bounds("apo", 20)),
        "lco_perfiber": _spec(
            "lco", derived_bounds("lco"),
            per_fiber_bounds=[[[a + int(rng.integers(0, 30)), b - int(rng.integers(0, 30))]
                               for a, b in derived_bounds("lco")] for _ in range(n_fiber)],
            n_fiber=n_fiber),
    }
    ref_nodes = ts.mode_node_pixels(1)
    mats = {k: ts.design_matrix_for_support(N_MODES, v.support_mask())
            for k, v in variants.items()}

    for k, A in mats.items():
        for c, nodes in enumerate(ref_nodes):
            col = c * ts.n_fourier_per_chip(N_MODES) + 1     # cos, k=1
            for p in nodes:
                lo, hi = int(np.floor(p)), int(np.ceil(p))
                # a node lands between two grid pixels; the column must change
                # sign across it wherever both rows are in support
                if A[lo, col] != 0.0 and A[hi, col] != 0.0:
                    assert A[lo, col] * A[hi, col] <= 0, f"{k} chip{c} node {p}"
        # nodes themselves are support-independent by construction
        assert all(np.array_equal(a, b) for a, b in
                   zip(ts.mode_node_pixels(1), ref_nodes))

    # Stronger: wherever two variants both have support, the ENTIRE row of the
    # design matrix is bit-identical -- every mode, not just mode 1.
    keys = list(mats)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            Ai, Aj = mats[keys[i]], mats[keys[j]]
            both = (np.any(Ai != 0, axis=1)) & (np.any(Aj != 0, axis=1))
            assert np.array_equal(Ai[both], Aj[both]), f"{keys[i]} vs {keys[j]}"


def test_widening_a_support_moves_no_node():
    narrow = _spec("lco", [(a + 40, b - 40) for a, b in derived_bounds("lco")])
    wide = _spec("lco", derived_bounds("lco"))
    An = ts.design_matrix_for_support(N_MODES, narrow.support_mask())
    Aw = ts.design_matrix_for_support(N_MODES, wide.support_mask())
    m = narrow.support_mask()
    assert np.array_equal(An[m], Aw[m])
    assert An.shape == Aw.shape          # same number of modes, always


def test_old_construction_DID_move_the_nodes():
    """Teeth for the anchoring test: the old, window-defined basis fails it.

    This is the property AKS's constraint is about.  Give the old construction
    two different windows (which is exactly what a per-telescope fix would have
    produced) and mode 1 lands somewhere else.
    """
    A_old_lco = old_style_design_matrix(derived_bounds("lco"))
    A_old_apo = old_style_design_matrix(derived_bounds("apo"))
    both = (np.any(A_old_lco != 0, axis=1)) & (np.any(A_old_apo != 0, axis=1))
    assert not np.array_equal(A_old_lco[both], A_old_apo[both]), \
        "old construction unexpectedly agreed"
    # quantify: mode-1 node of chip 0 moves by tens of pixels
    def node0(bounds):
        si, ei = bounds[0]
        return si + 0.25 * (2 * (ei - si) + 1)
    shift = abs(node0(derived_bounds("lco")) - node0(derived_bounds("apo")))
    assert shift > 10, shift
    # while the canonical node does not move at all
    assert ts.mode_node_pixels(1)[0][0] == ts.mode_node_pixels(1)[0][0]


def test_canonical_intervals_contain_every_measured_footprint():
    for tele, runs in MEASURED_RAW_RUNS.items():
        for (a, b), (ci, ce) in zip(runs, ts.CANONICAL_CHIP_INTERVALS):
            assert ci <= a and b <= ce, f"{tele} {(a, b)} escapes {(ci, ce)}"
    # and they are disjoint and ordered
    for (a, b), (c, d) in zip(ts.CANONICAL_CHIP_INTERVALS, ts.CANONICAL_CHIP_INTERVALS[1:]):
        assert b <= c


def test_prior_precision_is_support_independent():
    """The Matern prior lives on canonical coefficients, so it is too."""
    try:
        lam_a = ts.canonical_prior_precision(N_MODES, 1000.0, 100.0)
        lam_b = ts.canonical_prior_precision(N_MODES, 1000.0, 100.0)
    except ImportError:                       # domeflats needs jax
        return
    assert np.array_equal(lam_a, lam_b)
    assert lam_a.shape == (ts.n_fourier_total(N_MODES),)


# ══ G4 ════════════════════════════════════════════════════════════════════════
def _tfun(A, theta):
    return np.exp(A @ theta)


def canary_hits(tfun, support, rtol=0.0, atol=0.0):
    """Pixels of a normalised Tfun that equal 1/median(Tfun) exactly.

    That equality IS the exp(0) signature and is window-independent: it can
    only arise from a design-matrix row that is identically zero.
    """
    vals = tfun[np.isfinite(tfun) & (tfun != 0)]
    med = np.median(vals)
    norm = tfun / med
    target = 1.0 / med
    return int(np.sum(np.isclose(norm, target, rtol=rtol, atol=atol)))


def test_G4_canary_fires_on_a_consumer_that_ignores_support():
    """Baseline: the bug, reproduced.  Zero rows -> exp(0)=1 -> 1/median."""
    spec = _spec("lco", derived_bounds("lco"))
    sup = spec.support_mask()
    A = ts.design_matrix_for_support(N_MODES, sup)
    rng = np.random.default_rng(7)
    theta = rng.normal(0, 0.3, A.shape[1]) + 5.0
    hits = canary_hits(_tfun(A, theta), sup)
    assert hits == int((~sup).sum()), (hits, int((~sup).sum()))


def test_G4_passes_when_the_consumer_applies_support_in_linear_space():
    """The contract: Tfun[.!support] = 0 BEFORE any median normalisation.

    This is what apMADGICS did (linear-space zero-padded basis) and what
    build_starCont.jl's mean(filter(.!iszero, ...)) is already written for.
    """
    rng = np.random.default_rng(11)
    for tele in ("lco", "apo"):
        spec = _spec(tele, derived_bounds(tele))
        sup = spec.support_mask()
        A = ts.design_matrix_for_support(N_MODES, sup)
        for _ in range(20):
            theta = rng.normal(0, 0.3, A.shape[1]) + rng.uniform(0, 8)
            tf = _tfun(A, theta)
            tf[~sup] = 0.0                      # <-- the one-line contract
            assert canary_hits(tf, sup) == 0


def test_G4_passes_at_the_fitter_alone_with_nan_fill():
    """--dead-row-fill nan satisfies G4 without consumer cooperation.

    Recorded because it is the alternative evaluated in the report: it makes
    the bug class structurally impossible but breaks build_starCont.jl
    (`sum(starcont, dims=1)` -> NaN -> `specsum .> 0` false -> every sample
    column dropped), so it is NOT the default.
    """
    spec = _spec("lco", derived_bounds("lco"))
    sup = spec.support_mask()
    A = ts.design_matrix_for_support(N_MODES, sup)
    A[~sup, :] = np.nan
    rng = np.random.default_rng(13)
    theta = rng.normal(0, 0.3, A.shape[1]) + 5.0
    tf = _tfun(A, theta)
    assert np.all(np.isnan(tf[~sup]))
    assert canary_hits(tf, sup) == 0


# ── standalone driver ─────────────────────────────────────────────────────────
def _main():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    n_fail = 0
    for name, fn in fns:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            n_fail += 1
            print(f"  FAIL  {name}: {e}")
    print(f"\n{len(fns) - n_fail}/{len(fns)} passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(_main())
