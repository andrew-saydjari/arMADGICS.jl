#!/bin/bash
# build_metadata.sh -- decorate a COPY of the almanac file with target metadata (/meta,
# /meta_not_found, /meta_carton_pks) using `almanac add metadata`, restricted to the release's
# allExposure sdss_id set. Resumable: re-run the same command after an interruption.
#
# CHAIN POSITION: metadata step 2 of 2 (after make_restrict_ids.py; before alltargets/).
# SOURCE: work/2026_11_18/targets_040/scripts/build_040.sh + RUN_INFO.txt (0.4.0). For 0.4.0 the
#   code was the almanac fork's `metadata-resumable` branch (PR #14) at commit
#   4e3c218de6998355b452a39f9ff7363f0ff61375, Python 3.12.11, sdss-semaphore 0.2.8, 8 query
#   workers, batch size 2000, run on the head node against the SDSS database through an SSH tunnel.
#   0.4.0 decorated the copy IN PLACE (no --output); the copy is made here if it does not exist.
#
# Usage: build_metadata.sh <almanac.h5> <restrict_ids.npy> <workdir> [query_workers=8]
# Environment:
#   ALMANAC_SRC     path to an almanac checkout at the wanted commit (PYTHONPATH=$ALMANAC_SRC/src)
#   ALMANAC_PYTHON  Python interpreter with almanac's dependencies
#   PGPASSFILE      database password file (never printed); default ~/.pgpass
# Output: <workdir>/<basename of almanac>_meta.h5 (decorated copy), <workdir>/checkpoint/, build.log
set -euo pipefail
ALM=${1:?almanac file}; IDS=${2:?restrict ids .npy}; WD=${3:?work dir}; W=${4:-8}
: "${ALMANAC_SRC:?set ALMANAC_SRC to an almanac checkout}"; : "${ALMANAC_PYTHON:?set ALMANAC_PYTHON}"
export PYTHONPATH=$ALMANAC_SRC/src PGPASSFILE=${PGPASSFILE:-$HOME/.pgpass} PYTHONUNBUFFERED=1
mkdir -p "$WD/checkpoint"
OUT="$WD/$(basename "${ALM%.h5}")_meta.h5"
[ -f "$OUT" ] || cp "$ALM" "$OUT"
echo "[$(date -u +%FT%TZ)] build_metadata.sh start pid=$$ workers=$W commit=$(git -C "$ALMANAC_SRC" rev-parse HEAD)"
cat "$OUT" > /dev/null          # warm the page cache (0.4.0 did this; it made the run much faster)
exec "$ALMANAC_PYTHON" -m almanac.cli add metadata "$OUT" \
    --restrict-to "$IDS" \
    --checkpoint-dir "$WD/checkpoint" --keep-checkpoint \
    --query-workers "$W" --batch-size 2000
