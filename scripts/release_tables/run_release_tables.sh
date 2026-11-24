#!/bin/bash
# run_release_tables.sh -- build the release summary tables end to end (see README.md).
#
# Every input, output and release name is an environment variable; nothing is hard-coded.
# Products go to $OUTDIR (never into git); intermediates to $WORKDIR.
#
# Required:
#   RELEASE        release name, e.g. 0.4.1 (recorded in allTargets)
#   NAME           table name suffix, e.g. pass1_dr21 -> allExposure_pass1_dr21.h5
#   OUTDIR         where the release products are written
#   WORKDIR        intermediates (identity, exposure times, kernel outputs, fiber-offset solve, ...)
#   FULL_LIST      arMADGICS raw run full_list_info.h5
#   ALMANAC        almanac allobs_<start>_<end>.h5 the workup used
#   APRED          apred 1D output root
#   WORKUP         arMADGICS workup directory
#   META_FILE      almanac file decorated by `almanac add metadata` (metadata/build_metadata.sh);
#                  only needed for the alltargets step
# Exclusion lists, one of:
#   EXCLUSION_LISTS  existing directory with bad_nights.csv, bad_fibers.csv, rv_railed.csv, or
#   WAVECAL_DIR + TRACE_CSV  inputs to exclusions/build_exclusion_lists.jl (written to $OUTDIR)
# Optional:
#   FIBER_OFFSETS  existing fiber_offsets.csv: skip the solve and use it
#   IERS, IERS_B   pinned Earth-orientation tables: IERS-A (finals2000A format) and IERS-B (eopc04),
#                  combined as astropy IERS_Auto does (defaults: the vendored allexposure/iers/finals2000A.all
#                  and eopc04.1962-now, sha256 in allexposure/iers/SHA256SUMS). Never downloaded.
#   NUMERICS       julia (default, 0.4.1: pure Julia, no Python anywhere in the chain) |
#                  legacy040 (numpy 10**x through Python: ONLY for the 0.4.0 bit-regression, AVX-512 CPUs)
#   BARY           julia (default) | astropy (Python astropy kernel, cross-check only)
#   PYTHON         Python with requirements.txt installed; only needed for NUMERICS=legacy040 or BARY=astropy
#   THREADS        Julia threads (default: SLURM_CPUS_PER_TASK, else 8)
#   WORKERS        Distributed workers for the allExposure reads (default: SLURM_CPUS_ON_NODE, else 8)
#   STEPS          space-separated subset of: identity exclusions allexposure_a allvisitrv
#                  fiber_offsets allexposure offtarget restrict_ids alltargets cuts
#                  (default: all except restrict_ids, which build_metadata.sh needs before META_FILE exists)
#
# Slurm: this script never submits anything. For a batch run, wrap it in an sbatch script with
#   #SBATCH --constraint="[genoa|icelake|rome]"  and --cpus-per-task / --mem set (README.md).
#   NOTE: NUMERICS=legacy040 reproduces 0.4.0 bits only on AVX-512 CPUs (numpy uses Intel SVML there);
#   rome nodes have no AVX-512. The default (julia) gives the same bits on every Rusty CPU type.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
: "${RELEASE:?}" "${NAME:?}" "${OUTDIR:?}" "${WORKDIR:?}" "${FULL_LIST:?}" "${ALMANAC:?}" "${APRED:?}" "${WORKUP:?}"
PYTHON=${PYTHON:-python3}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
WORKERS=${WORKERS:-${SLURM_CPUS_ON_NODE:-8}}
NUMERICS=${NUMERICS:-julia}
BARY=${BARY:-julia}
IERS=${IERS:-$HERE/allexposure/iers/finals2000A.all}
IERS_B=${IERS_B-$HERE/allexposure/iers/eopc04.1962-now}
STEPS=${STEPS:-"identity exclusions allexposure_a allvisitrv fiber_offsets allexposure offtarget alltargets cuts"}
# Ceph + HDF5: same settings as the workup's run_workup.sh (lazy I/O, no POSIX file locking)
export CEPHTWEAKS_LAZYIO=${CEPHTWEAKS_LAZYIO:-1}
export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}
export RELEASE_TABLES_PYTHON=$PYTHON

JP="julia --project=$HERE"
JA="julia --project=$HERE/allexposure"
mkdir -p "$OUTDIR" "$WORKDIR"/{cache,visits,fiber_offsets,offtarget}
want() { [[ " $STEPS " == *" $1 "* ]]; }
log() { echo "[$(date +%FT%T)] $*"; }

AE_OUT="$OUTDIR/allExposure_${NAME}.h5"
AV_OUT="$OUTDIR/allVisitRV_${NAME}.h5"
AT_OUT="$OUTDIR/allTargets_${NAME}.h5"
LISTS=${EXCLUSION_LISTS:-$OUTDIR}

if want identity; then
    log "identity"
    [ -f "$WORKDIR/cache/identity.h5" ] || $JP -t "$THREADS" "$HERE/allexposure/legacy_python/01_build_identity.jl" \
        --full-list "$FULL_LIST" --out "$WORKDIR/cache/identity.h5"
fi

if want exclusions && [ -z "${EXCLUSION_LISTS:-}" ]; then
    log "exclusion lists"
    : "${WAVECAL_DIR:?set WAVECAL_DIR and TRACE_CSV, or EXCLUSION_LISTS}" "${TRACE_CSV:?}"
    $JP "$HERE/exclusions/build_exclusion_lists.jl" --wavecal-dir "$WAVECAL_DIR" --trace-csv "$TRACE_CSV" \
        --identity "$WORKDIR/cache/identity.h5" --workup "$WORKUP" --out "$OUTDIR"
fi

build_ae() {   # $1 = output, rest = extra args
    local out=$1; shift
    $JA -t "$THREADS" "$HERE/allexposure/bin/build_allexposure.jl" \
        --full-list "$FULL_LIST" --almanac "$ALMANAC" --apred "$APRED" --workup "$WORKUP" \
        --exclusion-lists "$LISTS" --cache-dir "$WORKDIR/cache" --workers "$WORKERS" \
        --python "$PYTHON" --numerics "$NUMERICS" --bary "$BARY" --iers "$IERS" --iers-b "$IERS_B" \
        --out "$out" --report "${out%.h5}.report.json" "$@"
}

# allExposure pass A: everything except the fiber zero points (they are fitted on allVisitRV,
# which is built from this table). Pass B reuses the cache and only re-joins and re-writes.
if want allexposure_a; then
    log "allExposure pass A (no fiber offsets)"
    build_ae "$WORKDIR/allExposure_A.h5"
fi

if want allvisitrv; then
    log "allVisitRV (21 columns)"
    $JP "$HERE/allvisitrv/01_seqmap.jl" "$ALMANAC" "$WORKDIR/visits/seqmap.jls"
    $JP -t "$THREADS" "$HERE/allvisitrv/02_build_visits.jl" --allexposure "$WORKDIR/allExposure_A.h5" \
        --seqmap "$WORKDIR/visits/seqmap.jls" --out "$WORKDIR/visits/v2_cols.jls"
    $JP "$HERE/allvisitrv/03_write_allvisitrv.jl" "$WORKDIR/visits/v2_cols.jls" "$WORKDIR/visits/allVisitRV_21.h5"
    $JP "$HERE/allvisitrv/04_attrs_allvisitrv.jl" "$WORKDIR/allExposure_A.h5" "$WORKDIR/visits/allVisitRV_21.h5"
fi

if want fiber_offsets; then
    F=$WORKDIR/fiber_offsets
    $JP -t "$THREADS" "$HERE/fiber_offsets/01_build_working_tables.jl" "$F" "$WORKDIR/allExposure_A.h5" \
        "$WORKDIR/visits/allVisitRV_21.h5" "$WORKDIR/visits/seqmap.jls"
    if [ -z "${FIBER_OFFSETS:-}" ]; then
        log "fiber offsets: solve"
        $JP -t "$THREADS" "$HERE/fiber_offsets/02_solve_offsets.jl" "$F"
        $JP "$HERE/fiber_offsets/03_offset_table.jl" "$F" "$LISTS/bad_fibers.csv"
        $JP "$HERE/fiber_offsets/04_release_header.jl" "$F"
        cp "$F/fiber_offsets_bare.csv" "$OUTDIR/fiber_offsets.csv"      # 0.4.0 shipped the bare file
        cp "$F/FIBER_OFFSETS_TABLE.md" "$OUTDIR/FIBER_OFFSETS_TABLE.md"
    else
        log "fiber offsets: using $FIBER_OFFSETS"
        cp "$FIBER_OFFSETS" "$F/fiber_offsets.csv"
        [ "$FIBER_OFFSETS" -ef "$OUTDIR/fiber_offsets.csv" ] || cp "$FIBER_OFFSETS" "$OUTDIR/fiber_offsets.csv"
    fi
    $JP -t "$THREADS" "$HERE/fiber_offsets/05_apply_allvisitrv.jl" "$F" "$WORKDIR/visits/allVisitRV_21.h5" \
        "$WORKDIR/visits/allVisitRV_24.h5" "$AV_OUT"
fi

if want allexposure; then
    log "allExposure (final, with fiber offsets)"
    build_ae "$AE_OUT" --fiber-offsets "$OUTDIR/fiber_offsets.csv"
fi

if want offtarget; then
    log "offtarget: sidecar CSV and allVisitRV column"
    O=$WORKDIR/offtarget
    $JP "$HERE/offtarget/01_almanac_scan.jl" "$ALMANAC" "$O/almanac_fps_flags.h5"
    $JP "$HERE/offtarget/02_sidecar_and_visits.jl" "$AE_OUT" "$WORKDIR/visits/allVisitRV_24.h5" \
        "$O/almanac_fps_flags.h5" "$OUTDIR/offtarget_exposures.csv" "$O/offtarget_visit_flag.h5"
    $JP "$HERE/offtarget/03_stage_allvisitrv.jl" --src "$WORKDIR/visits/allVisitRV_24.h5" \
        --flag "$O/offtarget_visit_flag.h5" --dst "$AV_OUT"
fi

if want restrict_ids; then
    $JP "$HERE/metadata/make_restrict_ids.jl" "$AE_OUT" "$WORKDIR/sdss_id_restrict.npy"
    log "next: metadata/build_metadata.sh $ALMANAC $WORKDIR/sdss_id_restrict.npy <meta workdir> (needs the SDSS database)"
fi

if want alltargets; then
    : "${META_FILE:?set META_FILE (output of metadata/build_metadata.sh)}"
    log "allTargets"
    julia --project="$HERE/alltargets" -t 1 "$HERE/alltargets/build_allTargets.jl" --release "$RELEASE" \
        --allexposure "$AE_OUT" --meta "$META_FILE" \
        --carton-table "$HERE/alltargets/reference/targetdb_carton.csv" \
        --semaphore-mapping "$HERE/alltargets/reference/sdss5_target_4_with_groups.csv" \
        --semaphore-version 0.2.9 --out "$AT_OUT"
fi

if want cuts; then
    $JP "$HERE/cuts/recommended_cuts.jl" "$AE_OUT"
fi
log "done; products in $OUTDIR"
