#!/bin/bash
# fetch_iers.sh -- re-derive the vendored IERS tables and verify them (provenance; not run by the build).
#
# The barycentric correction depends on Earth orientation (UT1-UTC, polar motion), interpolated in an
# IERS-A table. 0.4.0 used astropy's IERS_Auto at run time (2026-09-13) and did not record the table.
# IERS_Auto combines two files: the IERS-A table (finals2000A.all: rapid service + 1 year of predictions)
# and the IERS-B table (eopc04.1962-now: EOP 20 C04 final values), whose values replace the A table's
# Bulletin B columns. Tested on 2026-11-24 against all 26,715,540 rows of 0.4.0, the IERS-A table of
#   the PyPI wheel astropy-iers-data 0.2026.9.7.0.56.14 (uploaded 2026-09-07)
# reproduces every v_barycentric_correction bit (with the IERS-B table of that wheel or of 0.2026.8.31;
# they agree on the 0.4.0 dates). Other IERS-A tables do not: 0.2026.8.31 (the one installed) differs in
# 1,387 rows of MJD 61224-61230, 0.2026.8.24 in 4,318, 0.2026.9.14 and a table downloaded from
# datacenter.iers.org on 2026-11-24 in 90,224 rows of MJD 61222-61230.
# Both files of the 0.2026.9.7 wheel are vendored here and pinned for 0.4.1:
#   finals2000A.all  3.76 MB  values for MJD 41684-61652, measured (non-predicted) before MJD 61280
#   eopc04.1962-now  5.17 MB  final values to MJD 61252 (used up to the A table's last Bulletin B row, 61222)
# The builder never downloads anything; it reads this file (or the one given with --iers), records its
# path and sha256 in the product and the build report, and refuses a table whose sha256 does not match
# SHA256SUMS when one sits beside it.
#
# Usage: fetch_iers.sh OUTDIR     (writes OUTDIR/finals2000A.all and OUTDIR/eopc04.1962-now and
#                                         checks both against SHA256SUMS)
set -euo pipefail
OUT=${1:?usage: fetch_iers.sh OUTDIR}
HERE=$(cd "$(dirname "$0")" && pwd)
WHEEL_URL=https://files.pythonhosted.org/packages/fb/5f/7872fc61dc9cbf27abb5303967dcdbc016fa80bc5cd483064fe74b59ab02/astropy_iers_data-0.2026.9.7.0.56.14-py3-none-any.whl
WHEEL_SHA256=b24396f527feb4ec4b4ef0dcba820246ac7ecb1eed91ebb9e767d64774349659
mkdir -p "$OUT"
TMP=$(mktemp -d "$OUT/.fetch.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
curl -sS --fail --limit-rate 50M -o "$TMP/w.whl" "$WHEEL_URL"
echo "$WHEEL_SHA256  $TMP/w.whl" | sha256sum -c -
unzip -q -j "$TMP/w.whl" astropy_iers_data/data/finals2000A.all astropy_iers_data/data/eopc04.1962-now -d "$TMP"
mv "$TMP/finals2000A.all" "$TMP/eopc04.1962-now" "$OUT/"
(cd "$OUT" && sha256sum -c "$HERE/SHA256SUMS")
