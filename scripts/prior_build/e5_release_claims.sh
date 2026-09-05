#!/bin/bash
# E5: release STALE per-fiber build claims left by a dead run (claims whose
# owner host has no live build and whose three outputs do not all exist).
# Usage: ./e5_release_claims.sh <built_dir> <dead-hostname>
# Only run this when you are certain no build is running on <dead-hostname>.
set -u
BUILT=${1:?built_dir required}
DEADHOST=${2:?dead hostname required}
n=0
for c in "$BUILT"/.claims/[0-9][0-9][0-9]; do
  [ -d "$c" ] || continue
  fib=$(basename "$c")
  grep -q "^$DEADHOST " "$c/owner" 2>/dev/null || continue
  if [ -f "$BUILT/APOGEE_skycont_svd_30_f$fib.h5" ] && \
     [ -f "$BUILT/APOGEE_skyline_faint_svd_120_f$fib.h5" ] && \
     [ -f "$BUILT/APOGEE_skyline_faint_GSPICE_svd_120_f$fib.h5" ]; then
    continue  # fiber complete; claim harmless
  fi
  rm -rf "$c"
  n=$((n+1))
  echo "released stale claim: fiber $fib (owner $DEADHOST, outputs incomplete)"
done
echo "released $n stale claims in $BUILT/.claims"
