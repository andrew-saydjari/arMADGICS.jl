#!/bin/bash
# ------------------------------------------------------------------------------
# workup/test/test_module_setup.sh -- regression test for the lmod/nounset failure.
#
# 2026-11-17: the LCO 61360 daily died at "Running Spectra Workup" with
#     /mnt/sw/lmod/setup.sh: line 12: ENABLE_LMOD: unbound variable
# because run_workup.sh (set -euo pipefail) sourced the site lmod scripts in the minimal
# environment of the systemd-launched dailies, where `module` is not yet defined.
#
# This test extracts the module-setup section of workup/run_workup.sh VERBATIM (between
# the ">>> module-setup" / "<<< module-setup" markers) and runs it in a clean environment
# (env -i, set -euo pipefail), then asserts that:
#   1. it exits 0 and mpiexec is on PATH afterwards;
#   2. the strict options (errexit, nounset, pipefail) are back on afterwards;
#   3. lmod_module returns non-zero for a nonexistent module (so the FATAL check in
#      run_workup.sh still fires) and still restores the options.
#
# Skips (exit 0, prints SKIP) when the host has no lmod profile script, e.g. GitHub
# Actions runners. It therefore only proves anything on a Flatiron node.
#
# Usage: bash workup/test/test_module_setup.sh
# ------------------------------------------------------------------------------
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORKUP_DIR=$(cd "$HERE/.." && pwd)
SCRIPT=$WORKUP_DIR/run_workup.sh

have_lmod=0
for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
    [ -r "$f" ] && have_lmod=1
done
if [ "$have_lmod" = 0 ]; then
    echo "SKIP: no /etc/profile.d/modules.sh or lmod.sh on this host (no lmod); nothing to test"
    exit 0
fi

SECTION=$(sed -n '/^# >>> module-setup/,/^# <<< module-setup/p' "$SCRIPT")
if [ -z "$SECTION" ] || ! grep -q '^# <<< module-setup' <<<"$SECTION"; then
    echo "FAIL: could not extract the module-setup section from $SCRIPT (markers missing?)" >&2
    exit 1
fi

fails=0
check() {  # check <name> <expected-rc> <bash -c body>
    local name=$1 want=$2 body=$3 out rc=0
    out=$(env -i HOME="$HOME" USER="$USER" LOGNAME="$USER" PATH=/usr/bin:/bin \
            WORKUP_DIR="$WORKUP_DIR" SECTION="$SECTION" \
            bash -c "set -euo pipefail; $body" 2>&1) || rc=$?
    if [ "$rc" = "$want" ]; then
        echo "PASS: $name"
    else
        echo "FAIL: $name (exit $rc, expected $want)"
        fails=$((fails + 1))
    fi
    printf '%s\n' "$out" | sed 's/^/    | /'
}

# 1 + 2: the real section, clean env, strict mode; mpiexec found; options restored.
check "run_workup.sh module-setup in clean env (env -i, set -euo pipefail)" 0 '
    [ -z "$(type -t module)" ] || { echo "module unexpectedly predefined"; exit 9; }
    eval "$SECTION"
    echo "mpiexec: $(command -v mpiexec)"
    case $- in *e*) ;; *) echo "errexit not restored"; exit 8 ;; esac
    case $- in *u*) ;; *) echo "nounset not restored"; exit 8 ;; esac
    shopt -qo pipefail || { echo "pipefail not restored"; exit 8; }
    echo "options after section: $- (pipefail on)"
'

# 3: a failing module load is reported (non-zero) and options are still restored.
check "lmod_module on a nonexistent module returns non-zero, options restored" 0 '
    source "$(dirname "$WORKUP_DIR")/scripts/lib_lmod.sh"
    lmod_relax_opts
    for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
        [ -r "$f" ] && source "$f" && break
    done
    lmod_restore_opts
    if lmod_module load this-module-does-not-exist/0.0 2>/dev/null; then
        echo "lmod_module unexpectedly succeeded"; exit 7
    fi
    case $- in *e*) ;; *) echo "errexit not restored"; exit 8 ;; esac
    case $- in *u*) ;; *) echo "nounset not restored"; exit 8 ;; esac
    echo "failing load detected; options: $-"
'

if [ "$fails" -gt 0 ]; then
    echo "$fails test(s) FAILED"
    exit 1
fi
echo "all module-setup tests passed"
