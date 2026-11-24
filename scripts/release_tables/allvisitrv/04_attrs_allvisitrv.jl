# 04_attrs_allvisitrv.jl -- attach the description attributes to the 21-column allVisitRV file.
#
# CHAIN POSITION: allVisitRV step 4 of 4 (after 03_write_allvisitrv.jl).
# USAGE:   julia --project=scripts/release_tables 04_attrs_allvisitrv.jl ALLEXPOSURE.h5 ALLVISITRV.h5
#          (second file modified in place). ALLEXPOSURE supplies the `categories` coding of `category`.
# INPUTS:  ARGS[1] allExposure h5 (reads the `categories` attribute of `category`)
#          ARGS[2] allVisitRV h5 with 21 datasets and no dataset attributes yet
# OUTPUTS: ARGS[2], in place: a `description` attribute on 16 datasets and `categories` on `category`
# SOURCE:  work/2026_09_16/allvisitrv/s11_attrs_v2.py (0.4.0), via legacy_python/04_attrs_allvisitrv.py.
#
# Julia port of legacy_python/04_attrs_allvisitrv.py (2026-11-24); reproduces it exactly (attribute values and
# HDF5 string flavour, checked with compare_h5.py and compare_h5.jl). The descriptions are written as h5py
# writes a Python str: scalar, variable-length, UTF-8, NULLTERM. `categories` is copied with the
# HDF5 type it has in allExposure (there: a 1-D vlen UTF-8 array, which is also what h5py writes
# back). Only difference in principle: if allExposure ever stored `categories` as vlen ASCII,
# h5py would re-write it as UTF-8 while this copy keeps ASCII.
#
# The 14 metadata descriptions and the `categories` coding are taken from the exposure product /
# the standing text; the two scatter columns each get a description that states the relationship
# between them and which one to use when.  The five velocity columns carry no attributes.
using HDF5
include(joinpath(@__DIR__, "..", "lib", "h5attrs.jl"))

const VEL = ("sdss_id", "tele", "mjd", "mjd_mid", "v_rad_mean")

const PLAIN = "sample standard deviation (n-1) of the member v_rad, km/s. Defined for every row. " *
    "For n_exp == 2 it is |dv|/sqrt(2). It is biased low at small n_exp: under a Gaussian " *
    "E[s] = sqrt(2/pi)*sigma at n=2, about 20 percent low, and about 11 percent low at n=4. " *
    "Use v_rad_scatter_jack instead wherever it is finite; use this column for the " *
    "two-exposure visits, where the jackknife is undefined, and allow for the bias."

const JACK = "jackknifed (bias-corrected) scatter of the member v_rad, km/s, and the preferred scatter " *
    "wherever it is finite. " *
    "sigma_jack = N*s(x) - (N-1)*mean_i s(x without x_i), where " *
    "s(x) = sqrt( sum_i (x_i - xbar)^2 / (N-1) ) is the ordinary sample standard deviation and " *
    "the mean runs over the N leave-one-out subsamples. This is Eq. 8 of Saydjari et al. 2025, " *
    "AJ 169, 167 (arXiv:2408.07126), Section 4, applied to the statistic of their Eq. 10. " *
    "It removes the small-sample bias of the plain standard deviation in v_rad_scatter. " *
    "NaN for n_exp == 2, where the estimator is undefined: the leave-one-out subsamples of a " *
    "2-element sample are singletons, whose (N-1)-normalised standard deviation is 0/0. " *
    "4,703,259 of 6,670,336 rows (70.5 percent) are NaN for that reason; use v_rad_scatter " *
    "for those rows."

# insertion order as in the Python dict
const DESC = [
    "n_exp"              => "number of member exposures in the visit",
    "snr_mean"           => "mean of the member exposures' snr",
    "snr_visit"          => "quadrature sum of the member exposures' snr (visit-coadd-equivalent depth)",
    "v_rad_scatter"      => PLAIN,
    "v_rad_scatter_jack" => JACK,
    "e_v_rad_visit"      => "mean member e_v_rad / sqrt(n_exp), km/s; NAIVE propagation only - inherits the known missing high-S/N systematic floor, so it UNDERSTATES the true error at high S/N",
    "ra"                 => "right ascension, deg (constant within a visit, taken from the first member)",
    "dec"                => "declination, deg (constant within a visit)",
    "glon"               => "Galactic longitude, deg (constant within a visit)",
    "glat"               => "Galactic latitude, deg (constant within a visit)",
    "fiber"              => "product fiber index = adjfiberindx - 300*(tele-1), 1-300 (constant within a visit)",
    "era"                => "1 = plate era, 2 = FPS era",
    "field_id"           => "almanac field_id (-1 sentinel where absent)",
    "plate_id"           => "almanac plate_id (-1 sentinel in the FPS era)",
    "config_id"          => "almanac config_id (-1 sentinel in the plate era)",
    "category"           => "target class code; see the 'categories' attribute",
]

function add_attrs(src::AbstractString, new::AbstractString)
    h5open(src, "r") do s
        h5open(new, "r+") do f
            ks = keys(f)
            @assert length(ks) == 21 length(ks)
            @assert haskey(f, "v_rad_scatter") && haskey(f, "v_rad_scatter_jack")
            for (k, v) in DESC
                @assert haskey(f, k) && !(k in VEL) k
                put_str_attr!(f[k], "description", v)
            end
            copy_attr_raw!(s["category"], "categories", f["category"])
            for k in VEL
                @assert length(attributes(f[k])) == 0 "velocity dataset $k gained an attribute!"
            end
            ndesc = count(k -> haskey(attributes(f[k]), "description"), keys(f))
            println("21 datasets; $ndesc carry a description")
            cats = read_attribute(f["category"], "categories")
            println("category coding = [", join(("'" * c * "'" for c in cats), ", "), "]")   # Python list repr
            for k in ("v_rad_scatter", "v_rad_scatter_jack")
                println("\n$k:\n   ", read_attribute(f[k], "description"))
            end
        end
    end
    return nothing
end

function main(args)
    length(args) == 2 || error("usage: 04_attrs_allvisitrv.jl ALLEXPOSURE.h5 ALLVISITRV.h5")
    t = @elapsed add_attrs(args[1], args[2])
    println("work time: ", round(t; digits = 2), " s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
