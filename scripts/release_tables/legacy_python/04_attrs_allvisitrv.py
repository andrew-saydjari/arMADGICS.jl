# 04_attrs_allvisitrv.py -- attach the description attributes to the 21-column allVisitRV file.
#
# CHAIN POSITION: allVisitRV step 4 of 4 (after 03_write_allvisitrv.jl).
# USAGE:   python 04_attrs_allvisitrv.py ALLEXPOSURE.h5 ALLVISITRV.h5   (second file modified in place)
#          ALLEXPOSURE supplies the `categories` coding of `category`.
# SOURCE:  work/2026_09_16/allvisitrv/s11_attrs_v2.py (0.4.0). Logic and texts unchanged; the two
#          paths became argv[1], argv[2].
#
# The 14 metadata descriptions and the `categories` coding are taken from the exposure product /
# the standing text; the two scatter columns each get a description that states the relationship
# between them and which one to use when.  The five velocity columns carry no attributes.
import sys

import h5py

SRC = sys.argv[1]
NEW = sys.argv[2]
VEL = {"sdss_id", "tele", "mjd", "mjd_mid", "v_rad_mean"}

PLAIN = (
    "sample standard deviation (n-1) of the member v_rad, km/s. Defined for every row. "
    "For n_exp == 2 it is |dv|/sqrt(2). It is biased low at small n_exp: under a Gaussian "
    "E[s] = sqrt(2/pi)*sigma at n=2, about 20 percent low, and about 11 percent low at n=4. "
    "Use v_rad_scatter_jack instead wherever it is finite; use this column for the "
    "two-exposure visits, where the jackknife is undefined, and allow for the bias."
)

JACK = (
    "jackknifed (bias-corrected) scatter of the member v_rad, km/s, and the preferred scatter "
    "wherever it is finite. "
    "sigma_jack = N*s(x) - (N-1)*mean_i s(x without x_i), where "
    "s(x) = sqrt( sum_i (x_i - xbar)^2 / (N-1) ) is the ordinary sample standard deviation and "
    "the mean runs over the N leave-one-out subsamples. This is Eq. 8 of Saydjari et al. 2025, "
    "AJ 169, 167 (arXiv:2408.07126), Section 4, applied to the statistic of their Eq. 10. "
    "It removes the small-sample bias of the plain standard deviation in v_rad_scatter. "
    "NaN for n_exp == 2, where the estimator is undefined: the leave-one-out subsamples of a "
    "2-element sample are singletons, whose (N-1)-normalised standard deviation is 0/0. "
    "4,703,259 of 6,670,336 rows (70.5 percent) are NaN for that reason; use v_rad_scatter "
    "for those rows."
)

DESC = {
    "n_exp":         "number of member exposures in the visit",
    "snr_mean":      "mean of the member exposures' snr",
    "snr_visit":     "quadrature sum of the member exposures' snr (visit-coadd-equivalent depth)",
    "v_rad_scatter": PLAIN,
    "v_rad_scatter_jack": JACK,
    "e_v_rad_visit": "mean member e_v_rad / sqrt(n_exp), km/s; NAIVE propagation only - inherits the known missing high-S/N systematic floor, so it UNDERSTATES the true error at high S/N",
    "ra":            "right ascension, deg (constant within a visit, taken from the first member)",
    "dec":           "declination, deg (constant within a visit)",
    "glon":          "Galactic longitude, deg (constant within a visit)",
    "glat":          "Galactic latitude, deg (constant within a visit)",
    "fiber":         "product fiber index = adjfiberindx - 300*(tele-1), 1-300 (constant within a visit)",
    "era":           "1 = plate era, 2 = FPS era",
    "field_id":      "almanac field_id (-1 sentinel where absent)",
    "plate_id":      "almanac plate_id (-1 sentinel in the FPS era)",
    "config_id":     "almanac config_id (-1 sentinel in the plate era)",
    "category":      "target class code; see the 'categories' attribute",
}

with h5py.File(SRC, "r") as s:
    cats = s["category"].attrs["categories"]

with h5py.File(NEW, "r+") as f:
    assert len(f.keys()) == 21, len(f.keys())
    assert "v_rad_scatter" in f and "v_rad_scatter_jack" in f
    for k, v in DESC.items():
        assert k in f and k not in VEL, k
        f[k].attrs["description"] = v
    f["category"].attrs["categories"] = cats
    for k in VEL:
        assert len(f[k].attrs) == 0, "velocity dataset %s gained an attribute!" % k
    print("21 datasets; %d carry a description" % sum(1 for k in f if "description" in f[k].attrs))
    print("category coding =", list(f["category"].attrs["categories"]))
    for k in ("v_rad_scatter", "v_rad_scatter_jack"):
        print("\n%s:\n   %s" % (k, f[k].attrs["description"]))
