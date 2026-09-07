# shellcheck shell=bash
# ------------------------------------------------------------------------------
# scripts/lib_julia_warm.sh -- the ONE implementation of precompile-cache warming.
#
# SOURCE this file. Do NOT copy the function out of it into a launcher.
# The SLURM_NTASKS 2x bug (fixed 2026-09-07) was caused by copying a launcher body
# between scripts and changing the header without re-checking the body against it. Every
# copy of a shared block is a future divergence; this file exists so there is nothing to
# copy. Thin wrappers over a shared body, shared bodies over a sourced library.
#
# WHAT IT FIXES (MEASURED, Slurm job 6995969, 2026-09-07, 7 icelake nodes x 32 workers):
#   224 workers started cold and contended for the same precompile pidfiles on the shared
#   GPFS depot. 219 of them precompiled ApogeeReduction independently; Julia reported 92 s
#   for the first and 154-171 s for the later ones, i.e. the phase got WORSE as workers were
#   added. Job log: "Worker loading took 15 minutes, 3 seconds, 532 milliseconds".
#   One serial pass on the head node first makes every worker hit a warm cache.
#
# COST WHEN ALREADY WARM (MEASURED, ccalin051, julia 1.11.0, 2026-09-07):
#   targeted load of the sky-prior package set   8.7 s   (274 s the first, cold, time)
#   Pkg.precompile() over the whole manifest    30.2 s   (219 s cold; AR_WARM_FULL=1 only)
#
# ASSUMPTION, stated so it can be checked rather than believed: the Julia depot is on
# SHARED storage, so warming on the head node serves every node in the allocation. That is
# true at Flatiron (the pidfiles in the job log above are under /mnt/home/.../.julia) and at
# Utah. If a site ever moves the depot to node-local disk, head-only warming buys nothing
# for the other nodes -- warm_precompile.jl prints DEPOT_PATH[1] into the job log so the
# assumption is visible there.
#
# Usage (from a launcher, after base_dir and julia_version are known):
#     source "$base_dir/scripts/lib_julia_warm.sh"
#     warm_julia_precompile "$julia_version" "$base_dir" \
#         "$base_dir/scripts/prior_build/e5_sky_run.jl" "$base_dir/src"
# Arguments:
#     $1  julia version channel, e.g. 1.11.0, or "" to use whatever `julia` resolves to.
#         MUST match the version the job then runs: precompile caches are keyed on the exact
#         Julia build, so warming with 1.11.6 does nothing for a job that runs
#         `julia +1.11.0` (MEASURED: juliaup default here is 1.11.6, every launcher pins
#         1.11.0 -- get this wrong and the warm step is pure cost with no benefit).
#     $2  project directory (passed as --project)
#     $3+ files to scan. Pass the DRIVER; the warmer follows its include() graph, so there
#         is no package list to keep in sync with the driver. Directories are also accepted
#         and are walked for .jl files.
# Env:
#     AR_WARM_PRECOMPILE=0  skip warming entirely
#     AR_WARM_FULL=1        also Pkg.precompile() the whole manifest (+~30 s warm)
#
# NEVER fails the caller. Launchers run under `set -e`; a warming failure must not take a
# job down, so every path here returns 0.
# ------------------------------------------------------------------------------

warm_julia_precompile() {
    local jver="$1"
    local proj="$2"
    shift 2

    if [ "${AR_WARM_PRECOMPILE:-1}" = "0" ]; then
        echo "warm: AR_WARM_PRECOMPILE=0, skipping precompile warming"
        return 0
    fi

    local warmer="$proj/scripts/warm_precompile.jl"
    if [ ! -f "$warmer" ]; then
        echo "warm: WARNING $warmer not found, skipping (job continues unwarmed)"
        return 0
    fi

    local t_start=$SECONDS
    echo "warm: serial precompile warm on $(hostname) before any workers spawn"
    # `|| true` twice over: once here and once as exit(0) inside the Julia script. Warming
    # is best-effort by design.
    julia +"$jver" --project="$proj" "$warmer" "$@" || \
        echo "warm: WARNING warmer exited non-zero; continuing unwarmed"
    echo "warm: done in $((SECONDS - t_start)) s"
    return 0
}
