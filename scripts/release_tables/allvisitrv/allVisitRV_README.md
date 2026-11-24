# allVisitRV — pass-1 DR21

`allVisitRV_pass1_dr21.h5` is a visit-level radial-velocity table: 6,670,336 rows, 24 columns,
one row per star per visit. It is built from `allExposure_pass1_dr21.h5` in this directory by
averaging each star's exposures within a visit. It contains only visits whose exposures split
evenly between the two dither positions, so the dither-dependent velocity term cancels in the
mean. Velocities are barycentric, in km/s. Three columns carry velocities with the per-fiber
zero-point offsets removed; see **Corrected velocities** below.

Rows are sorted by `(tele, mjd, sdss_id, mjd_mid)`. Datasets are chunked at 65,536 rows with the
shuffle filter and gzip level 9. `h5py` reads `tele` as `|S3` (`b"apo"`). The file is 512 MB.
Every column except the five velocity columns carries a `description` attribute.

---

## Columns

| column | type | what it is |
|---|---|---|
| `sdss_id` | Int64 | star identifier (`sdss_id` from the exposure product, not `sdss_id_workup`) |
| `tele` | `S3` | `"apo"` or `"lco"` |
| `mjd` | Int32 | SDSS observing-night MJD |
| `mjd_mid` | Float64 | mean of the member exposures' `mjd_mid_exposure` — the epoch the velocity refers to |
| `v_rad_mean` | Float64 | unweighted mean of the member `v_rad`, km/s, barycentric |
| `n_exp` | Int16 | number of member exposures. Always even, ≥ 2, max 44 |
| `snr_mean` | Float64 | mean of the member `snr` |
| `snr_visit` | Float64 | quadrature sum of the member `snr` — the depth an equivalent visit coadd would have |
| `v_rad_scatter` | Float64 | sample standard deviation (n−1) of the member `v_rad`, km/s. Defined for every row |
| `v_rad_scatter_jack` | Float64 | jackknifed (debiased) version of the same scatter, km/s. NaN when `n_exp == 2` |
| `e_v_rad_visit` | Float64 | mean member `e_v_rad` / √`n_exp`, km/s. Naive propagation — read the caveat below |
| `v_rad_corr_mean` | Float64 | `v_rad_mean` with the per-fiber zero-point offset removed, km/s. NaN on four fibers |
| `v_rad_corr_scatter` | Float64 | sample standard deviation (n−1) of the member `v_rad_corr`, km/s |
| `v_rad_corr_scatter_jack` | Float64 | jackknifed version of the same, km/s. NaN when `n_exp == 2` |
| `ra`, `dec` | Float64 | ICRS position, deg |
| `glon`, `glat` | Float64 | Galactic coordinates, deg |
| `fiber` | Int16 | fiber index, `adjfiberindx - 300*(tele-1)`, 1–300 |
| `era` | Int8 | 1 = plate era, 2 = FPS era |
| `field_id` | Int32 | almanac field id (−1 where absent) |
| `plate_id` | Int32 | almanac plate id (−1 in the FPS era) |
| `config_id` | Int32 | almanac configuration id (−1 in the plate era) |
| `category` | UInt8 | target class: `0 ""`, `1 "science"`, `2 "standard_apogee"`, `3 "sky_apogee"`, `4 "unplugged"`. The coding is also the `categories` attribute on the column. Only 1 and 2 occur |

`ra`, `dec`, `glon`, `glat`, `fiber`, `era`, `field_id`, `plate_id`, `config_id` and `category`
are the same for every exposure in a visit, so they are carried across rather than averaged. This
was checked on all 6,670,336 visits, not sampled.

`(sdss_id, tele, mjd, mjd_mid)` is unique. `(sdss_id, tele, mjd)` is not: 418,610 rows are second
or later visits to the same star on the same night.

---

## Which scatter column to use

Use `v_rad_scatter_jack` where it is finite. It is the plain standard deviation with its
small-sample bias removed, following Saydjari et al. 2025, AJ 169, 167 (arXiv:2408.07126),
Section 4. With `x` the member velocities and `N` = `n_exp`:

```
s(x)       = sqrt( (1/(N-1)) * sum_i ( x_i - mean(x) )^2 )       plain standard deviation
sigma_jack = N * s(x) - (N-1) * (1/N) * sum_i s(x without x_i)   debiased (their Eq. 8)
```

The correction is largest for short visits: it raises the median scatter by 12 percent at
`n_exp = 4`, 7 percent at 6, 5 percent at 8, and less beyond that.

`v_rad_scatter_jack` is **undefined when `n_exp == 2`**, because dropping one of two exposures
leaves a single point with no standard deviation. Two-exposure visits are 70.5 percent of the
table (4,703,259 rows), so most of the column is NaN. For those rows use `v_rad_scatter`, which
is `|Δv|/√2`, and allow for the fact that it **underestimates the true scatter by about
20 percent** — for Gaussian noise `E[s] = √(2/π)·σ` at `n = 2`. If you need an approximate
correction, multiply by `√(π/2)` ≈ 1.25.

`v_rad_scatter` is defined for every row, so it is the column to use if you want one scatter
definition across the whole table. It is biased low at every small `n_exp`: about 11 percent at
`n = 4`, less as `n` grows.

Neither column depends on the pipeline error model, which is what makes them a useful empirical
alternative to `e_v_rad_visit`. Both are dominated by real velocity variation for binaries and by
measurement noise for everything else.

| | `v_rad_scatter` | `v_rad_scatter_jack` |
|---|---|---|
| rows with a value | 6,670,336 | 1,967,077 (`n_exp ≥ 4`) |
| median, km/s | 0.0696 | 0.137 |
| range, km/s | 3.1e-8 – 436 | 3.0e-4 – 435 |

The two medians differ mostly because they cover different rows. On the 1,967,077 rows where both
are defined the medians are 0.128 and 0.137.

---

## Corrected velocities

Each fiber has its own radial-velocity zero point. `v_rad_corr_mean` is `v_rad_mean` with that
zero point subtracted:

```
v_rad_corr_mean = v_rad_mean - offset(tele, fiber, era)
```

The offsets are in `fiber_offsets.csv` in this directory, one per `(telescope, fiber, era)`. They
are an **empirical calibration**. They were measured from stars observed on more than one fiber,
telescope or era: those stars tie the fibers together, and a joint least-squares fit for the star
velocities and the fiber zero points at the same time gives the offsets. Stars seen on only one
fiber carry no information about the offsets and contribute nothing to the fit.

One consequence matters. Because a single set of offsets covers both telescopes and both eras,
the **APO/LCO difference and the plate→FPS step are absorbed into the offsets too**. You do not
need to apply a separate telescope or era correction to `v_rad_corr_mean`. The reason the
telescopes differ by about 250 m/s is not understood; the correction removes the difference
without explaining it.

`v_rad_corr_scatter` and `v_rad_corr_scatter_jack` follow the same two definitions as
`v_rad_scatter` and `v_rad_scatter_jack`, applied to the corrected member velocities. **They are
almost identical to the uncorrected columns** — medians 0.0691 and 0.1367 km/s against 0.0696 and
0.137. A visit lives on one fiber, so a constant per-fiber offset cancels in the scatter within a
visit. The columns are here so that `v_rad_corr_mean` has a matching spread, not because the
correction changes them. The correction shows up between visits, not inside one.

**Four fibers could not be calibrated**, all at LCO: `fiber` 121, 240, 270 and 182 (FIBERID 180,
61, 31 and 119). The first three have velocities far outside the range of every other fiber and
too few linking stars to pin down; these are the same three fibers Saydjari et al. 2025 (AJ 169,
167) could not calibrate and dropped from the DR17 analysis. The fourth, `fiber` 182, is not
obviously bad — it simply has too little data in the FPS era to measure its offset. All four
carry NaN in all three corrected columns. That is **36,872 rows**, 0.55 percent of the table. The
uncorrected columns on those rows are unaffected.

---

## What a visit is

A visit is one star's exposures inside one **almanac object sequence**. An object sequence is an
uninterrupted run of science exposures on a single pointing, taken from
`raw/<observatory>/<mjd>/sequences/objects` in the pass-1 almanac. It breaks whenever a
calibration exposure interrupts the run, an exposure number is missing, or the pointing changes.
The survey has 48,572 such sequences over 3,757 nights, and every row of the exposure product
falls inside one.

Visits are grouped on `(object sequence, sdss_id)`. Rows with `sdss_id <= 0` are sky or
unassigned fibers and are dropped.

**Member exposures are masked first.** The release default mask is

```
qa_bitmask == 0  &&  ingestBit == 0  &&  RV_flag == 0  &&  bary_flag == 0  &&  isfinite(v_rad)
```

which keeps 25,516,821 of 26,715,540 exposure rows. See `RECOMMENDED_CUTS.md`.

**Only balanced visits are kept.** After masking, the surviving exposures must sit at exactly two
distinct dither positions in equal numbers. AB, ABBA, ABAB and so on all qualify. If the mask
removed an exposure and broke the balance, the whole visit is dropped — the surviving partner is
not discarded to rebalance it.

The dither position is the commanded `dithered_pixels` value. Only eight values occur across the
survey and at most two on any night, so the raw value is already the position label. A value of 0
means the commanded position was not recorded; a visit containing one cannot be shown to be
balanced and is excluded.

### Counts

| | visits |
|---|---|
| grouped by `(object sequence, sdss_id)` | 7,244,864 |
| no exposure survived the mask | −238,039 |
| all surviving exposures at one dither position | −260,945 |
| uneven split between the two positions | −74,541 |
| dither position not recorded | −1,003 |
| **in this file** | **6,670,336** |

3,908,568 distinct stars; 3,272,198 rows from APO and 3,398,138 from LCO.

Exposures per visit: `n_exp` = 2 is 70.5 percent of rows, 4 is 13.0, 8 is 7.7, 6 is 4.3, 12 is
1.6, 10 is 1.5, and everything above 12 is 1.4 percent together.

1,442,791 rows are plate era and 5,227,545 are FPS era. All plate-era rows have
`config_id == -1`; 5,160,637 FPS-era rows have `plate_id == -1`. The other 66,908 FPS-era rows
are plate-plugged configurations observed after the switch and carry both ids. No row has both
sentinels.

---

## Times and the mean

`mjd_mid` is the mean of the member exposures' `mjd_mid_exposure`, the frame-averaged UTC
mid-time on chip c. That is the epoch at which the barycentric correction was evaluated, so it is
the time `v_rad_mean` refers to. `mjd` is the observing-night MJD, so `mjd_mid - mjd` runs from
about −0.19 to +0.69 d; the two are not meant to agree more closely than that.

`v_rad_mean` is a plain unweighted mean, and so is `mjd_mid`. This is deliberate. The
dither-dependent velocity term has opposite sign at the two dither positions, so equal weights
cancel it exactly in a balanced visit. Inverse-variance weighting gives the higher-S/N position
more weight and leaves a residual proportional to the imbalance, which puts back the systematic
the ABBA cadence exists to remove. Inverse-variance weighting is the better choice for visits of
arbitrary composition; this table has only balanced visits, so plain means are correct here.

`v_rad_corr_mean` is a plain unweighted mean of the corrected member velocities, over the same
members, for the same reason.

---

## On error bars

`e_v_rad_visit` is `mean(e_v_rad)/√n_exp` and nothing more. It inherits the pipeline's missing
systematic floor at high S/N, so it understates the true uncertainty for bright stars. Use it as
a relative weight or a rough cut, not as a calibrated error bar. The scatter columns are the
empirical alternative.

Ranges, for reference: `snr_mean` 0.013 – 629 (median 60.2), `snr_visit` 0.020 – 1655
(median 96.4), `e_v_rad_visit` 5.3e-4 – 8.2e4 km/s (median 0.0252). The low-S/N and large-error
tails are real rows that pass the default mask. They are what `snr_mean`, `snr_visit` and
`e_v_rad_visit` are here to let you cut on.
