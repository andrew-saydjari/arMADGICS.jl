# shellcheck shell=bash
# ------------------------------------------------------------------------------
# scripts/lib_lmod.sh -- run the site's lmod setup and `module` from strict-mode scripts.
#
# SOURCE this file (same convention as lib_julia_warm.sh: one implementation, no copies).
#
# WHY: our launchers run under `set -euo pipefail`. The Flatiron lmod setup
# (/etc/profile.d/modules.sh -> /mnt/sw/lmod/setup.sh, Ansible-managed, changed with the
# Rocky 9 rebuild) and the `module` shell function read variables that may be unset
# (ENABLE_LMOD, FPATH, __Init_Default_Modules, ...). Sourcing them under `nounset` aborts:
#     /mnt/sw/lmod/setup.sh: line 12: ENABLE_LMOD: unbound variable
# That is what stopped the LCO 61360 daily (2026-11-17) at "Running Spectra Workup".
# It never showed up while the dailies ran from a shell that already had `module` defined;
# from the systemd service the environment is minimal, so the profile scripts get sourced.
# Exporting the variables is not a fix (the next unset one trips it); `nounset` -- and,
# defensively, `errexit` -- must be OFF while site code runs, and restored afterwards.
#
# API:
#   lmod_relax_opts / lmod_restore_opts
#       Bracket TOP-LEVEL code that sources site scripts (sourcing must stay at top level,
#       not inside a function, so any `declare` in site code stays global). Saves whether
#       -e / -u were on, turns both off, and restores exactly that state. Not nestable.
#   lmod_ensure_module
#       NOT provided as a function on purpose (see above); the caller writes:
#           lmod_relax_opts
#           if ! type module &>/dev/null; then
#               for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
#                   [ -r "$f" ] && source "$f" && break
#               done
#           fi
#           lmod_restore_opts
#   lmod_module <args...>
#       Runs `module <args...>` with -e/-u off, restores the caller's options, and returns
#       module's exit status, so `if ! lmod_module load X; then FATAL...; fi` still works.
# ------------------------------------------------------------------------------

_lmod_saved_e=0
_lmod_saved_u=0

lmod_relax_opts() {
    case $- in *e*) _lmod_saved_e=1 ;; *) _lmod_saved_e=0 ;; esac
    case $- in *u*) _lmod_saved_u=1 ;; *) _lmod_saved_u=0 ;; esac
    set +eu
}

lmod_restore_opts() {
    if [ "$_lmod_saved_e" = 1 ]; then set -e; fi
    if [ "$_lmod_saved_u" = 1 ]; then set -u; fi
    return 0
}

lmod_module() {
    local had_e=0 had_u=0 rc=0
    case $- in *e*) had_e=1 ;; esac
    case $- in *u*) had_u=1 ;; esac
    set +eu
    module "$@" || rc=$?
    if [ "$had_e" = 1 ]; then set -e; fi
    if [ "$had_u" = 1 ]; then set -u; fi
    return "$rc"
}
