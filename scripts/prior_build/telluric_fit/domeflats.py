
import numpy as np
import jax.numpy as jnp
from astropy.io import fits

def get_telluric_models():

    # Load all molecular spectra and stack into one matrix
    files = ["CH4.fits", "CO2.fits", "H2O.fits"]
    all_tellurics = []
    labels = []
    for f in files:
        data = fits.getdata(f).astype(np.float64)
        n, m, p = data.shape
        all_tellurics.append(np.hstack([np.ones((n*m, 252)), data.reshape(n * m, p)]))
        labels.append(f"{f}: ({n}, {m}, {p})")

    P = 19500
    wl = 15100 + 0.1 * np.arange(P)

    wl_pad = np.arange(1, 253)[::-1] * -0.1
    wl = np.concatenate([wl_pad + wl[0], wl])
    return (wl, all_tellurics)




def fourier_design_matrix_1d(N, P):
    """(N, 2P+1) Fourier design matrix over half the period."""
    t = jnp.arange(N) / (2 * N + 1)
    cols = [jnp.ones(N)]
    for k in range(1, P + 1):
        cols.append(jnp.cos(2 * jnp.pi * k * t))
        cols.append(jnp.sin(2 * jnp.pi * k * t))
    return jnp.column_stack(cols)


def matern_prior_variance_1d(P, s):
    """Matérn-1 power-spectral-density prior over Fourier modes."""
    freqs = jnp.arange(P + 1)
    psd   = 1.0 / (1.0 + s**2 * (2 * jnp.pi * freqs)**2)
    return jnp.concatenate([jnp.array([psd[0]]), jnp.repeat(psd[1:], 2)])
