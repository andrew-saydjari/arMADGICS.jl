"""Batch dome-flat continuum + sky fit  (successor to 20260323.py).

Model (fibre f, pixel p):

    log_flux[f, p] = A_cont[p] @ theta[f]
                   - gamma[f, :n_abs] @ T[:n_abs, p]   (absorption)
                   + gamma[f, n_abs:] @ T[n_abs:, p]   (emission)

WHAT CHANGED vs 20260323.py
---------------------------
20260323.py built the continuum design matrix from a single hardcoded,
telescope-blind, epoch-blind literal::

    CHIP_INDICES = [(340-20, 3430+20), (3675-20, 6240+20), (6430-20, 8490+20)]
    _A_cont_np = np.zeros((N_PIXELS, n_fourier))
    for i, (si, ei) in enumerate(CHIP_INDICES):
        _A_cont_np[si:ei, sj:ej] = fourier_design_matrix_1d(ei - si, N_MODES)

That literal did two jobs at once: it chose which pixels were fit AND it set the
Fourier period/phase (the basis was built ON the window).  Outside it the design
matrix is identically zero, so ``Tfun = exp(A @ theta) = exp(0) = 1.0`` exactly
-- a placeholder that downstream divides by ``median(Tfun)`` and consumes as a
throughput of ~0.003.  It was 44 px too red at LCO's blue chip edge, which is
where the 67 dead LCO prior pixels came from.

Here the two jobs are separated (see telluric_support.py):

  * the BASIS is built once on FIXED canonical intervals, identical for every
    telescope, fiber and epoch, so mode k always means the same thing and
    ``theta`` stays comparable;
  * the SUPPORT is derived from the dome-flat stack per telescope (optionally
    per fiber) and selects which ROWS take part -- never which columns exist.

and the support is DECLARED in the output (``support`` dataset + provenance
attrs, asserted at write time to equal the nonzero-row support of
``design_matrix``), so a consumer can tell "no basis here" from "the throughput
really is that small".

Consumer note
-------------
``design_matrix`` rows outside the support are exact zeros (unchanged
convention, ``--dead-row-fill nan`` available for audit runs -- see README).
exp() cannot represent "no support", so a consumer MUST apply ``support`` in
linear space::

    Tfun = exp.(Atell * theta);  Tfun[.!support] .= 0.0

which restores the apMADGICS semantics that ``nanzeromedian`` and
``build_starCont.jl``'s ``mean(filter(.!iszero, ...))`` are already written for.

Output layout:
    design_matrix (N_PIXELS, n_fourier)           f4
    support       (N_PIXELS,) or (N_PIXELS, N_FIBERS)  bool   <-- NEW
    theta         (N_FILES, N_FIBERS, n_fourier)  f4
    gamma         (N_FILES, N_FIBERS, N_TELLURIC) f4
    chi_sq        (N_FILES, N_PIXELS)             f4
    median_resid  (N_FILES, N_PIXELS)             f4
    chi_sq_fiber  (N_FILES, N_FIBERS)             f4
    stage         (N_FILES,)                      i1
    T_final       (N_FILES, N_TELLURIC, N_PIXELS) f4   only with --update-t
    paths         (N_FILES,)                      bytes
    success       (N_FILES,)                      bool
    attrs: n_fourier, n_abs, n_em, N_TELLURIC, N_PIXELS, n_modes,
           canonical_chip_intervals, support_bounds, support_telescope,
           support_min_live_frac, support_min_exposure_frac,
           support_edge_buffer, support_input_list, support_input_list_sha256,
           support_n_exposures_scanned, support_derivation_version,
           support_per_fiber, off_support_fill

Usage:
    python telluric_support.py --list list_lco_full.txt --telescope lco \
        --out support_lco.json                       # once per telescope
    python fit_domeflats.py --support support_lco.json --input list.txt \
        --output out.h5 --no-figures --s2 --s2-iters 50 --lambda-t 1e4 \
        --update-t --s-pixels 1000
"""

import argparse
import time
from pathlib import Path

import jax
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp
import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import pickle
from functools import partial
import telluric_support as ts

# ── Hyperparameters ────────────────────────────────────────────────────────────
N_PIXELS    = ts.N_PIXELS
N_MODES     = 32
S_PIXELS    = 1000       # Matern length-scale in pixels (handles gap in middle chip)
AMPLITUDE   = 100       # Matern amplitude
SMALL       = 1e-12
HUBER_K     = 5.0
FISTA_ITERS = 100       # FISTA iterations per NNLS sub-problem
S1_IRLS     = 8         # stage-1 outer IRLS iterations
S2_ITERS    = 20        # stage-2 outer alternating iterations (default)

# ── I/O defaults (overridable via CLI) ────────────────────────────────────────
_DEFAULT_INPUT   = "/mnt/home/asaydjari/ceph/scratch/2026_02_16/ar1Dfname_list.txt"
_DEFAULT_OUTPUT  = "domeflat_fits.h5"
_DEFAULT_FIGDIR  = "figures"

# ── Design matrices (support applied by configure()) ──────────────────────────
# NOTE (2026-09-07): 20260323.py also built, at import time,
#     _wl, _ = get_telluric_models();  K_lsf = instrument_lsf_sparse_matrix(...)
# Neither name is referenced anywhere else in that file -- dead code that read
# ~6 MB of FITS and built an 8700-row sparse LSF matrix on every import.
# Removed; nothing numerical depends on it.
resampled_wl = ts.RESAMPLED_WL

n_fourier_per_chip = ts.n_fourier_per_chip(N_MODES)
n_fourier          = ts.n_fourier_total(N_MODES)

# Populated by configure().  Declared here so the module imports cleanly (the
# tests import it without a support file).
SUPPORT_SPEC        = None
support_mask        = None     # (N_PIXELS,) global support (union over fibers)
support_mask_full   = None     # (N_PIXELS,) or (N_PIXELS, N_FIBERS) as declared
support_bounds      = None
_A_cont_np          = None
A_cont              = None
Lambda              = None
Lambda_diag         = None
_chip_pixel_indices = None
_fiber_support_j    = None     # (N_FIBERS, N_PIXELS) float weights, or None
OFF_SUPPORT_FILL    = "zero"


def configure(spec, s_pixels=S_PIXELS, amplitude=AMPLITUDE,
              dead_row_fill="zero", template_dir=".", load_t=True):
    """Install a support spec: build the design matrix, prior and index sets.

    Must be called before any jitted update function is traced.
    """
    global SUPPORT_SPEC, support_mask, support_mask_full, support_bounds
    global _A_cont_np, A_cont, _chip_pixel_indices, _fiber_support_j
    global OFF_SUPPORT_FILL

    if load_t:
        load_templates(template_dir)
    SUPPORT_SPEC      = spec
    support_bounds    = [tuple(b) for b in spec.bounds]
    support_mask_full = spec.support_mask(N_PIXELS)
    support_mask      = (support_mask_full.any(axis=1) if support_mask_full.ndim == 2
                         else support_mask_full)
    OFF_SUPPORT_FILL  = dead_row_fill

    # ROWS are selected; COLUMNS (the canonical basis) are never touched.
    _A_cont_np = ts.design_matrix_for_support(
        N_MODES, support_mask, intervals=tuple(map(tuple, spec.canonical_intervals)),
        n_pixels=N_PIXELS,
    )
    # G1 is a property of the zero-filled matrix; the nan variant is a strictly
    # louder relabelling of the same rows, checked against the same mask.
    ok, _, _ = ts.check_support_matches_design(support_mask, _A_cont_np)
    if not ok:
        raise RuntimeError("G1 failed at configure(): support != nonzero rows")

    A_cont = jnp.array(_A_cont_np, dtype=jnp.float64)
    _chip_pixel_indices = np.where(support_mask)[0]

    if support_mask_full.ndim == 2:
        # Per-fiber support is applied through the WEIGHTS, not through a
        # per-fiber design matrix: zeroing cinv drops the row from that fiber's
        # normal equations exactly, at no memory cost, and leaves the columns
        # (hence the mode definitions) untouched.
        _fiber_support_j = jnp.array(support_mask_full.T.astype(np.float64))
    else:
        _fiber_support_j = None

    _build_lambda(s_pixels, amplitude)


def _build_lambda(s_pixels, amplitude):
    """Matern prior precision over the CANONICAL modes (see telluric_support)."""
    global Lambda, Lambda_diag
    lam = ts.canonical_prior_precision(
        N_MODES, s_pixels, amplitude,
        intervals=tuple(map(tuple, SUPPORT_SPEC.canonical_intervals)),
    )
    Lambda      = jnp.array(lam, dtype=jnp.float64)
    Lambda_diag = jnp.diag(Lambda)


def design_matrix_for_output():
    """design_matrix as written to h5 (applies --dead-row-fill)."""
    A = np.array(_A_cont_np, dtype=np.float32)
    if OFF_SUPPORT_FILL == "nan":
        A[~support_mask, :] = np.nan
    return A


# ── Sky templates ──────────────────────────────────────────────────────────────
# Loaded from CWD by load_templates(), which configure() calls.  Deferred (they
# were module-level in 20260323.py) so the module can be imported for tests from
# a directory that does not hold the frozen .pkl inputs.
n_abs = n_em = N_TELLURIC = None
T_init = _A_sky_abs_np = _A_sky_em_np = _T_signed = _TT_mat = None


def load_templates(template_dir="."):
    global n_abs, n_em, N_TELLURIC, T_init
    global _A_sky_abs_np, _A_sky_em_np, _T_signed, _TT_mat
    d = Path(template_dir)
    with open(d / "tellurics_init.pkl", "rb") as fp:
        _A_sky_abs_np = -pickle.load(fp)["H"].T   # (N_PIXELS, n_abs)
    with open(d / "nmf_sky_lines.pkl", "rb") as fp:
        _A_sky_em_np = pickle.load(fp).T          # (N_PIXELS, n_em)
    n_abs      = _A_sky_abs_np.shape[1]
    n_em       = _A_sky_em_np.shape[1]
    N_TELLURIC = n_abs + n_em
    with open(d / "T_init.pkl", "rb") as fp:
        T_init = pickle.load(fp)   # (N_TELLURIC, N_PIXELS), already normalised
    # Precompute T-dependent gram matrices (valid when T = T_init).
    _T_signed = jnp.concatenate([-T_init[:n_abs].T, T_init[n_abs:].T], axis=1)
    _TT_mat   = (_T_signed[:, :, None] * _T_signed[:, None, :]).reshape(N_PIXELS, -1)


# ── FISTA (non-negative least squares) ────────────────────────────────────────
def _nnls_fista(H, g, x_init):
    """FISTA for  min_{x >= 0}  0.5 x^T H x - g^T x.  Warm-started from x_init."""
    alpha = 1.0 / (jnp.linalg.norm(H, ord="fro") + 1e-12)
    x0    = jnp.maximum(x_init, 0.0)

    def step(carry, _):
        x, y, t = carry
        x_new = jnp.maximum(y - alpha * (H @ y - g), 0.0)
        t_new = 0.5 * (1.0 + jnp.sqrt(1.0 + 4.0 * t * t))
        y_new = x_new + ((t - 1.0) / t_new) * (x_new - x)
        return (x_new, y_new, t_new), None

    (x, _, _), _ = jax.lax.scan(step, (x0, x0, 1.0), None, length=FISTA_ITERS)
    return x


# ── Fast update functions (T fixed = T_init, precomputed gram matrix) ─────────
@jax.jit
def _update_theta(log_f, cinv, gamma):
    """Batched weighted ridge solve for theta (T fixed = T_init)."""
    correction = gamma[:, :n_abs] @ T_init[:n_abs] - gamma[:, n_abs:] @ T_init[n_abs:]
    rhs   = log_f + correction
    g_all = (cinv * rhs) @ A_cont
    H_all = jnp.einsum("fp,pm,pn->fmn", cinv, A_cont, A_cont) + Lambda_diag
    return jnp.linalg.solve(H_all, g_all[..., None])[..., 0]   # (N_F, n_fourier)


@jax.jit
def _update_gamma(D, cinv, gamma_prev):
    """NNLS per fibre for all sky/telluric amplitudes (T fixed = T_init)."""
    g_all = (cinv * D) @ _T_signed                               # (N_F, K)
    H_all = (cinv @ _TT_mat).reshape(-1, N_TELLURIC, N_TELLURIC) # (N_F, K, K)
    return jax.vmap(_nnls_fista)(H_all, g_all, gamma_prev)       # (N_F, K)


@jax.jit
def _huber_cinv(log_f, cinv, theta, gamma):
    """Huber-reweighted cinv (T fixed = T_init)."""
    correction = gamma[:, :n_abs] @ T_init[:n_abs] - gamma[:, n_abs:] @ T_init[n_abs:]
    r = (log_f - theta @ A_cont.T + correction) * jnp.sqrt(cinv)
    w = jnp.where(jnp.abs(r) <= HUBER_K, 1.0, HUBER_K / (jnp.abs(r) + 1e-12))
    return cinv * w


# ── General update functions (T passed explicitly, used when updating T) ───────
@jax.jit
def _update_theta_T(log_f, cinv, gamma, T):
    correction = gamma[:, :n_abs] @ T[:n_abs] - gamma[:, n_abs:] @ T[n_abs:]
    rhs   = log_f + correction
    g_all = (cinv * rhs) @ A_cont
    H_all = jnp.einsum("fp,pm,pn->fmn", cinv, A_cont, A_cont) + Lambda_diag
    return jnp.linalg.solve(H_all, g_all[..., None])[..., 0]


@jax.jit
def _update_gamma_T(D, cinv, gamma_prev, T):
    K        = N_TELLURIC
    T_signed = jnp.concatenate([-T[:n_abs].T, T[n_abs:].T], axis=1)
    g_all    = (cinv * D) @ T_signed
    TT_mat   = (T_signed[:, :, None] * T_signed[:, None, :]).reshape(N_PIXELS, -1)
    H_all    = (cinv @ TT_mat).reshape(-1, K, K)
    return jax.vmap(_nnls_fista)(H_all, g_all, gamma_prev)


@jax.jit
def _huber_cinv_T(log_f, cinv, theta, gamma, T):
    correction = gamma[:, :n_abs] @ T[:n_abs] - gamma[:, n_abs:] @ T[n_abs:]
    r = (log_f - theta @ A_cont.T + correction) * jnp.sqrt(cinv)
    w = jnp.where(jnp.abs(r) <= HUBER_K, 1.0, HUBER_K / (jnp.abs(r) + 1e-12))
    return cinv * w


@jax.jit
def _update_T(D, cinv, gamma, T_ref, T_prev, lambda_t):
    """Direct solve per support pixel with L2 prior toward T_ref.

    Off-support pixels stay at T_ref.  (Was "per chip pixel" against the old
    hardcoded window; it now follows the declared support.)
    """
    K            = N_TELLURIC
    Gamma_signed = jnp.concatenate([-gamma[:, :n_abs], gamma[:, n_abs:]], axis=1)
    D_chip       = D[:, _chip_pixel_indices]
    cinv_chip    = cinv[:, _chip_pixel_indices]
    T_ref_chip   = T_ref[:, _chip_pixel_indices]
    g_all  = (Gamma_signed.T @ (cinv_chip * D_chip)).T + lambda_t * T_ref_chip.T
    GG     = Gamma_signed[:, :, None] * Gamma_signed[:, None, :]
    H_all  = (cinv_chip.T @ GG.reshape(-1, K * K)).reshape(-1, K, K) + lambda_t * jnp.eye(K)
    T_chip = jax.vmap(jnp.linalg.solve)(H_all, g_all)
    T_chip = jnp.maximum(T_chip, 0.0)
    return T_ref.at[:, _chip_pixel_indices].set(T_chip.T)


# ── Stage 1: IRLS, all K components, T fixed ─────────────────────────────────
def stage1_irls(log_f, cinv):
    """Cold-start IRLS fit: theta + all gamma components, T fixed = T_init."""
    n_f   = log_f.shape[0]
    gamma = jnp.full((n_f, N_TELLURIC), SMALL)
    theta = _update_theta(log_f, cinv, gamma)   # cold start from plain WLS
    for _ in range(S1_IRLS):
        cinv_rob = _huber_cinv(log_f, cinv, theta, gamma)
        theta    = _update_theta(log_f, cinv_rob, gamma)
        D        = log_f - theta @ A_cont.T
        gamma    = _update_gamma(D, cinv_rob, gamma)
    return theta, gamma


# ── Stage 2: optional, more iterations, optional T update ─────────────────────
@partial(jax.jit, static_argnames=["n_iter", "update_t"])
def joint_fit(log_f, cinv, theta_init, gamma_init, n_iter=S2_ITERS,
              update_t=False, lambda_t=1e4):
    """Continue IRLS from a warm start; optionally update T."""
    theta = theta_init
    gamma = gamma_init
    T     = jnp.array(T_init)
    T_ref = jnp.array(T_init)   # fixed reference for the L2 prior

    for _ in range(n_iter):
        if update_t:
            cinv_rob = _huber_cinv_T(log_f, cinv, theta, gamma, T)
            theta    = _update_theta_T(log_f, cinv_rob, gamma, T)
            D        = log_f - theta @ A_cont.T
            gamma    = _update_gamma_T(D, cinv_rob, gamma, T)
            T        = _update_T(D, cinv_rob, gamma, T_ref, T, lambda_t)
        else:
            cinv_rob = _huber_cinv(log_f, cinv, theta, gamma)
            theta    = _update_theta(log_f, cinv_rob, gamma)
            D        = log_f - theta @ A_cont.T
            gamma    = _update_gamma(D, cinv_rob, gamma)

    return theta, gamma, T


# ── Per-file preprocessing and fit ────────────────────────────────────────────
def preprocess(path):
    with h5.File(path, "r") as f:
        flux = f["flux_1d"][:]
        ivar = f["ivar_1d"][:]
        mask = f["mask_1d"][:].astype(bool) & (ivar > 0)
    # Off-support pixels never enter the fit: their design-matrix rows are zero,
    # so they contribute nothing to the normal equations anyway, but masking
    # here also keeps chi_sq/median_resid honest about what was fit.
    mask &= support_mask[None, :]
    ivar[~mask] = 0.0
    flux[~mask] = 0.0
    flux   = np.clip(flux, SMALL, None)
    snr_sq = (flux ** 2 * ivar).astype(np.float64)
    log_f  = np.log(flux).astype(np.float64)
    c      = 0.5 / np.maximum(snr_sq, SMALL)
    log_f += c
    ivar_log = np.where(mask, snr_sq / (1.0 + c), 0.0)
    cinv = jnp.array(ivar_log)
    if _fiber_support_j is not None:
        cinv = cinv * _fiber_support_j          # per-fiber row selection
    return flux, ivar, mask, jnp.array(log_f), cinv


def fit_file(path, run_s2=False, s2_iters=S2_ITERS, update_t=False, lambda_t=1e4):
    """Fit one dome-flat file.  Returns (theta, gamma, T_out, flux, ivar, mask)."""
    flux, ivar, mask, log_f, cinv = preprocess(path)
    n_f = log_f.shape[0]

    theta, gamma = stage1_irls(log_f, cinv)
    T_out = T_init   # numpy array; unchanged unless stage 2 updates it

    if run_s2:
        theta, gamma, T_j = joint_fit(
            log_f, cinv, theta, gamma,
            n_iter=s2_iters, update_t=update_t, lambda_t=lambda_t,
        )
        # E3/P2 FIX (2026-09-01): restore the T_out reassignment that was
        # commented out in 20260323.py (the "T_final == T_init on rerun" bug).
        T_out = np.array(T_j, dtype=np.float64)

    return theta, gamma, T_out, flux, ivar, mask


# ── QA figure ─────────────────────────────────────────────────────────────────
def make_figure(path, theta_np, gamma_np, flux, ivar, mask, T_use=None, settings=None):
    """Four-row QA figure for one dome flat (panels follow the derived support)."""
    if T_use is None:
        T_use = T_init

    n_f     = flux.shape[0]
    theta_j = jnp.array(theta_np, dtype=jnp.float64)
    gamma_j = jnp.array(gamma_np, dtype=jnp.float64)
    T_j     = jnp.array(T_use,    dtype=jnp.float64)

    correction = np.array(
        gamma_j[:, :n_abs] @ T_j[:n_abs] - gamma_j[:, n_abs:] @ T_j[n_abs:]
    )
    model     = np.exp(np.array(theta_j @ A_cont.T) - correction)
    continuum = np.exp(np.array(theta_j @ A_cont.T))
    sky       = (
        np.exp(-np.array(gamma_j[:, :n_abs] @ T_j[:n_abs]))
        * np.exp( np.array(gamma_j[:, n_abs:] @ T_j[n_abs:]))
    )

    resid = (flux - model) * np.sqrt(np.where(mask, ivar, 0.0))
    resid[~mask] = np.nan

    fig = plt.figure(figsize=(24, 13))
    gs  = gridspec.GridSpec(
        4, 5, figure=fig, height_ratios=[2.5, 1, 1, 1],
        hspace=0.45, wspace=0.3, top=0.94,
    )
    ax_data      = fig.add_subplot(gs[0, 0])
    ax_continuum = fig.add_subplot(gs[0, 1])
    ax_sky       = fig.add_subplot(gs[0, 2])
    ax_residuals = fig.add_subplot(gs[0, 3])
    ax_z_hist    = fig.add_subplot(gs[0, 4])

    gs_spec = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, :], wspace=0.08)
    ax_spec_chips = [fig.add_subplot(gs_spec[0, c]) for c in range(3)]

    gs_rand = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2, :], wspace=0.08)
    gs_261  = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[3, :], wspace=0.08)
    ax_rand_chips = [fig.add_subplot(gs_rand[0, c]) for c in range(3)]
    ax_261_chips  = [fig.add_subplot(gs_261[0,  c]) for c in range(3)]

    for ax, data, title in [
        (ax_data,      flux,      "Flux"),
        (ax_continuum, continuum, "Transfer function model"),
        (ax_sky,       sky,       "Sky model"),
    ]:
        valid  = data[mask]
        vmin, vmax = (np.nanpercentile(valid, [1, 99]) if valid.size else (0, 1))
        ax.imshow(data, aspect="auto", origin="lower", vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_ylabel("Fiber")
        ax.set_xticks([])

    im = ax_residuals.imshow(
        resid, aspect="auto", origin="lower", vmin=-3, vmax=3,
        cmap="RdBu_r", interpolation="nearest",
    )
    plt.colorbar(im, ax=ax_residuals, fraction=0.046, location="top")
    ax_residuals.set_title("Residuals (σ)")
    ax_residuals.set_xticks([])

    z  = resid[np.isfinite(resid)]
    xg = np.linspace(-5, 5, 100)
    ax_z_hist.hist(z,        bins=xg, histtype="step", density=True, color="#666666", label="nominal σ")
    ax_z_hist.hist(z / 2.3,  bins=xg, histtype="step", density=True, color="k",       label="2.3× σ", ls="--")
    ax_z_hist.plot(xg, np.exp(-0.5 * xg ** 2) / np.sqrt(2 * np.pi),
                   c="tab:red", ls="--", lw=0.8, label="N(0,1)")
    ax_z_hist.set_xlim(xg[[0, -1]])
    ax_z_hist.set_xlabel("Residual (σ)")
    ax_z_hist.set_title("Residual distribution")
    ax_z_hist.set_yticks([])
    ax_z_hist.legend(fontsize=8)

    with np.errstate(all="ignore"):
        med_resid = np.nanmedian(resid, axis=0)
    for chip_idx, (ax, (si, ei)) in enumerate(zip(ax_spec_chips, support_bounds)):
        wl_c = resampled_wl[si:ei]
        ax.plot(wl_c, med_resid[si:ei], lw=0.5, c="k")
        ax.axhline(0, c="tab:blue", ls="--", lw=0.8)
        ax.set_xlim(wl_c[[0, -1]])
        ax.set_ylim(-3, +3)
        ax.set_xlabel("Wavelength (Å)")
        if chip_idx == 0:
            ax.set_ylabel("Median residual (σ)")
        else:
            ax.set_yticklabels([])
        if chip_idx == 1:
            ax.set_title("Median residual spectrum (all fibres)")

    total_cinv = ivar.sum(axis=1)
    valid_f    = np.where(total_cinv > 0)[0]
    rng        = np.random.default_rng(42)
    rand_f     = int(rng.choice(valid_f)) if valid_f.size else 0
    f261       = min(261, n_f - 1)

    for axes_chips, fi, label in [
        (ax_rand_chips, rand_f, f"Random fibre (index {rand_f})"),
        (ax_261_chips,  f261,   f"Fibre {f261}" + ("  [gap in middle chip]" if f261 == 261 else "")),
    ]:
        d = flux[fi].copy();     d[~mask[fi]] = np.nan
        m = model[fi].copy()
        c = continuum[fi].copy()

        for chip_idx, (ax, (si, ei)) in enumerate(zip(axes_chips, support_bounds)):
            wl_c = resampled_wl[si:ei]
            ax.plot(wl_c, d[si:ei], lw=0.4, c="k", alpha=0.8,
                    label="data" if chip_idx == 2 else None)
            flux_ylim = ax.get_ylim()
            ax.plot(wl_c, m[si:ei], lw=0.9, c="tab:red",
                    label="model" if chip_idx == 2 else None)
            ax.plot(wl_c, c[si:ei], lw=0.7, c="tab:orange", ls="--",
                    label="continuum" if chip_idx == 2 else None)
            ax.set_ylim(flux_ylim)
            ax.set_xlim(wl_c[[0, -1]])
            ax.set_xlabel("Wavelength (Å)")
            if chip_idx == 0:
                ax.set_ylabel("Flux")
            else:
                ax.set_yticklabels([])
            if chip_idx == 1:
                ax.set_title(label)
            if chip_idx == 2:
                ax.legend(fontsize=8, loc="upper right")

    s = settings or {}
    settings_str = (
        f"  {SUPPORT_SPEC.telescope}  support={support_bounds}"
        f"  s_pixels={s.get('s_pixels', S_PIXELS)}"
        f"  amplitude={s.get('amplitude', AMPLITUDE)}"
        f"  n_modes={N_MODES}"
        f"  s2={s.get('s2', False)}"
        f"  s2_iters={s.get('s2_iters', S2_ITERS)}"
        f"  update_t={s.get('update_t', False)}"
        f"  lambda_t={s.get('lambda_t', 1e4)}"
    )
    fig.suptitle(f"{Path(path).name}:  {settings_str}", fontsize=8)
    fig.tight_layout()
    return fig


# ── Output-file provenance ────────────────────────────────────────────────────
def write_support_and_provenance(fout, A_out):
    """Create `design_matrix` + `support` and stamp the provenance attrs.

    Asserts G1 (support == nonzero-row support of design_matrix) AT WRITE TIME:
    the whole point of the declaration is that it cannot drift from the matrix.
    """
    spec = SUPPORT_SPEC
    fout.create_dataset("design_matrix", data=A_out)
    fout.create_dataset("support", data=support_mask_full)

    # G1, on exactly the bytes that go to disk.
    if OFF_SUPPORT_FILL == "nan":
        nz = ~np.all(np.isnan(A_out), axis=1)
    else:
        nz = np.any(A_out.astype(np.float64) != 0.0, axis=1)
    if not np.array_equal(nz, support_mask):
        raise RuntimeError(
            "G1 FAILED at write time: declared support != design_matrix support "
            f"({int(nz.sum())} nonzero rows vs {int(support_mask.sum())} declared)"
        )

    a = fout.attrs
    a["n_modes"]                    = N_MODES
    a["canonical_chip_intervals"]   = np.array(spec.canonical_intervals, dtype=np.int32)
    a["support_bounds"]             = np.array(spec.bounds, dtype=np.int32)
    a["support_telescope"]          = spec.telescope
    a["support_min_live_frac"]      = spec.min_live_frac
    a["support_min_exposure_frac"]  = spec.min_exposure_frac
    a["support_edge_buffer"]        = spec.edge_buffer
    a["support_input_list"]         = spec.input_list
    a["support_input_list_sha256"]  = spec.input_list_sha256
    a["support_n_exposures_scanned"] = spec.n_exposures_scanned
    a["support_n_exposures_in_list"] = spec.n_exposures_in_list
    a["support_derivation_version"] = spec.version
    a["support_per_fiber"]          = bool(spec.per_fiber)
    a["off_support_fill"]           = OFF_SUPPORT_FILL
    a["off_support_note"] = (
        "exp() cannot represent 'no support': apply `support` in LINEAR space, "
        "Tfun[.!support] = 0, before any median normalisation."
    )


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Batch dome-flat fit.")
    parser.add_argument("--support",    required=True,
                        help="support JSON from telluric_support.py (per telescope)")
    parser.add_argument("--input",      default=_DEFAULT_INPUT)
    parser.add_argument("--output",     default=_DEFAULT_OUTPUT)
    parser.add_argument("--figures",    default=_DEFAULT_FIGDIR)
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--start",      type=int, default=0)
    parser.add_argument("--end",        type=int, default=None)
    parser.add_argument("--overwrite",  action="store_true")
    parser.add_argument("--s2",         action="store_true")
    parser.add_argument("--s2-iters",   type=int,   default=S2_ITERS)
    parser.add_argument("--update-t",   action="store_true")
    parser.add_argument("--lambda-t",   type=float, default=1e4)
    parser.add_argument("--s-pixels",   type=float, default=S_PIXELS)
    parser.add_argument("--amplitude",  type=float, default=AMPLITUDE)
    parser.add_argument("--dead-row-fill", choices=["zero", "nan"], default="zero",
                        help="value written into off-support design_matrix rows. "
                             "'zero' keeps the delivered convention (consumers "
                             "MUST apply `support` in linear space); 'nan' makes "
                             "an ignoring consumer fail loudly but BREAKS "
                             "build_starCont.jl until it does. Audit runs only.")
    args = parser.parse_args()

    spec = ts.SupportSpec.from_json(args.support)
    configure(spec, s_pixels=args.s_pixels, amplitude=args.amplitude,
              dead_row_fill=args.dead_row_fill)

    figure_dir = Path(args.figures)
    if not args.no_figures:
        figure_dir.mkdir(exist_ok=True)

    paths   = np.loadtxt(args.input, dtype=str)
    N_FILES = len(paths)

    with h5.File(paths[0], "r") as _f:
        N_FIBERS = _f["flux_1d"].shape[0]

    if spec.per_fiber and spec.n_fiber != N_FIBERS:
        raise SystemExit(f"support has n_fiber={spec.n_fiber}, data has {N_FIBERS}")

    print(f"{N_FILES} files  |  N_FIBERS={N_FIBERS}  |  output → {args.output}")
    print(f"n_fourier={n_fourier}  n_abs={n_abs}  n_em={n_em}  N_TELLURIC={N_TELLURIC}")
    print(f"support[{spec.telescope}] = {support_bounds}  "
          f"({int(support_mask.sum())} px)  per_fiber={spec.per_fiber}  "
          f"edge_buffer={spec.edge_buffer}")
    print(f"canonical = {[tuple(x) for x in spec.canonical_intervals]}  (basis, FROZEN)")
    print(f"stage2={args.s2}  update_t={args.update_t}  "
          f"s2_iters={args.s2_iters}  lambda_t={args.lambda_t}")

    A_out = design_matrix_for_output()

    if not Path(args.output).exists():
        with h5.File(args.output, "w") as fout:
            write_support_and_provenance(fout, A_out)
            fout.create_dataset("theta",
                (N_FILES, N_FIBERS, n_fourier), dtype="f4",
                fillvalue=np.nan, chunks=(1, N_FIBERS, n_fourier))
            fout.create_dataset("gamma",
                (N_FILES, N_FIBERS, N_TELLURIC), dtype="f4",
                fillvalue=np.nan, chunks=(1, N_FIBERS, N_TELLURIC))
            fout.create_dataset("chi_sq",
                (N_FILES, N_PIXELS), dtype="f4",
                fillvalue=np.nan, chunks=(1, N_PIXELS))
            fout.create_dataset("median_resid",
                (N_FILES, N_PIXELS), dtype="f4",
                fillvalue=np.nan, chunks=(1, N_PIXELS))
            fout.create_dataset("chi_sq_fiber",
                (N_FILES, N_FIBERS), dtype="f4",
                fillvalue=np.nan, chunks=(1, N_FIBERS))
            fout.create_dataset("stage", data=np.zeros(N_FILES, dtype=np.int8))
            fout.create_dataset("paths",
                data=np.array([str(p).encode() for p in paths]))
            fout.create_dataset("success", data=np.zeros(N_FILES, dtype=bool))
            fout.attrs["n_fourier"]  = n_fourier
            fout.attrs["n_abs"]      = n_abs
            fout.attrs["n_em"]       = n_em
            fout.attrs["N_TELLURIC"] = N_TELLURIC
            fout.attrs["N_PIXELS"]   = N_PIXELS
        print("Created output file.")
    else:
        # Resuming: the support baked into the file must match the one we were
        # handed, or the rows already fit mean something different.
        with h5.File(args.output, "r") as fchk:
            if "support" not in fchk:
                raise SystemExit(
                    f"{args.output} predates the support declaration; refusing to "
                    "append incompatible rows. Fit to a fresh output file."
                )
            if not np.array_equal(fchk["support"][:], support_mask_full):
                raise SystemExit(
                    f"{args.output} was fit with a DIFFERENT support; refusing to "
                    "mix. Use the same --support JSON or a fresh output file."
                )
            if not np.array_equal(
                np.asarray(fchk.attrs["canonical_chip_intervals"]),
                np.array(spec.canonical_intervals, dtype=np.int32),
            ):
                raise SystemExit(f"{args.output} used different canonical intervals")
        print("Resuming from existing output file (support matches).")

    if args.update_t:
        with h5.File(args.output, "r+") as fout:
            if "T_final" not in fout:
                fout.create_dataset("T_final",
                    (N_FILES, N_TELLURIC, N_PIXELS), dtype="f4",
                    fillvalue=np.nan, chunks=(1, N_TELLURIC, N_PIXELS))
                print("Added T_final dataset.")

    with h5.File(args.output, "r") as fchk:
        success_all = np.array(fchk["success"])

    start = args.start
    end   = args.end if args.end is not None else N_FILES

    todo_all = (np.arange(start, end) if args.overwrite
                else np.where(~success_all[start:end])[0] + start)
    if todo_all.size == 0:
        print("All requested files already done (use --overwrite to re-run).")
        return

    warmup_idx = int(todo_all[0])
    print(f"JIT warmup (file {warmup_idx}) …", flush=True)
    t0 = time.time()
    _warmup_result = fit_file(
        paths[warmup_idx], run_s2=args.s2, s2_iters=args.s2_iters,
        update_t=args.update_t, lambda_t=args.lambda_t,
    )
    print(f"  done in {time.time() - t0:.1f}s")

    n_done   = int(success_all[:start].sum())
    n_failed = 0
    t_loop   = time.time()
    stage_val = 2 if args.s2 else 1

    with h5.File(args.output, "r+") as fout:
        success = np.array(fout["success"])

        for i in range(start, end):
            if success[i] and not args.overwrite:
                n_done += 1
                continue

            t0 = time.time()
            try:
                if i == warmup_idx:
                    theta, gamma, T_out, flux, ivar, mask = _warmup_result
                else:
                    theta, gamma, T_out, flux, ivar, mask = fit_file(
                        paths[i], run_s2=args.s2, s2_iters=args.s2_iters,
                        update_t=args.update_t, lambda_t=args.lambda_t,
                    )

                if flux.shape[0] != N_FIBERS:
                    raise ValueError(
                        f"N_FIBERS mismatch: got {flux.shape[0]}, expected {N_FIBERS}"
                    )

                T_use = T_out if args.update_t else T_init
                model = np.exp(
                    theta.astype(np.float64) @ np.array(_A_cont_np).T
                    - (gamma[:, :n_abs].astype(np.float64) @ T_use[:n_abs]
                       - gamma[:, n_abs:].astype(np.float64) @ T_use[n_abs:])
                )

                fout["theta"][i]        = theta
                fout["gamma"][i]        = gamma
                fout["chi_sq"][i]       = np.mean((flux - model) ** 2 * ivar, axis=0)
                fout["median_resid"][i] = np.median((flux - model) * np.sqrt(ivar), axis=0)
                _n_good = np.maximum((ivar > 0).sum(axis=1), 1)
                fout["chi_sq_fiber"][i] = (
                    np.sum((flux - model) ** 2 * ivar, axis=1) / _n_good
                ).astype(np.float32)
                fout["stage"][i]        = stage_val
                fout["success"][i]      = True
                if args.update_t:
                    fout["T_final"][i]  = T_out

                n_done += 1

                if not args.no_figures:
                    fig = make_figure(
                        paths[i], theta, gamma, flux, ivar, mask,
                        T_use=T_use if args.update_t else None,
                        settings=dict(
                            s_pixels=args.s_pixels, amplitude=args.amplitude,
                            s2=args.s2, s2_iters=args.s2_iters,
                            update_t=args.update_t, lambda_t=args.lambda_t,
                        ),
                    )
                    fig_path = figure_dir / (Path(paths[i]).stem + ".png")
                    fig.savefig(fig_path, dpi=120, bbox_inches="tight")
                    plt.close(fig)

                elapsed = time.time() - t0
                rate    = n_done / max(time.time() - t_loop, 1e-6)
                eta     = (todo_all.size - n_done) / rate if rate > 0 else float("inf")
                print(
                    f"  [{i+1:5d}/{N_FILES}]  {Path(paths[i]).name}"
                    f"  {elapsed:.2f}s"
                    f"  ({rate:.2f} files/s  ETA {eta/60:.0f}min)",
                    flush=True,
                )

            except Exception as exc:
                n_failed += 1
                print(
                    f"  [{i+1:5d}/{N_FILES}]  FAILED  {Path(paths[i]).name}: {exc}",
                    flush=True,
                )

    print(f"\nDone. {n_done} succeeded, {n_failed} failed.")


if __name__ == "__main__":
    main()
