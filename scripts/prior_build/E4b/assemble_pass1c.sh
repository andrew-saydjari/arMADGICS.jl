#!/bin/bash
# pass-1c corpus + built-prior view under the FINAL AKS-approved (2026-09-05)
# tfunlist cut policy. Supersedes starCont_pass1b for pass-1 (pass-1b retained;
# nothing deleted). APO: all 300 fibers resampled + rebuilt under the final
# policy. LCO: pass-1b reused wholesale (final lco list is byte-identical to
# the intelligent list) — links resolve pass-1b's targets so pass-1c stands
# alone. MANIFEST.md is written next to this view by the pass-1c run.
set -euo pipefail
NEW_S=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_final/tell_prior_disk
NEW_B=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_final/built
P1B=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_pass1b
ROOT=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_pass1c
mkdir -p "$ROOT/tell_prior_disk" "$ROOT/built_apo" "$ROOT/built_lco"

for i in $(seq 1 300); do
    i3=$(printf "%03d" "$i")
    ln -sfn "$NEW_S/starCont_$i3.jdat" "$ROOT/tell_prior_disk/starCont_$i3.jdat"
    ln -sfn "$NEW_B/APOGEE_starcont_svd_60_f$i3.h5" "$ROOT/built_apo/APOGEE_starcont_svd_60_f$i3.h5"
done
for i in $(seq 301 600); do
    i3=$(printf "%03d" "$i")
    ln -sfn "$(readlink -f "$P1B/tell_prior_disk/starCont_$i3.jdat")" "$ROOT/tell_prior_disk/starCont_$i3.jdat"
    ln -sfn "$(readlink -f "$P1B/built_lco/APOGEE_starcont_svd_60_f$i3.h5")" "$ROOT/built_lco/APOGEE_starcont_svd_60_f$i3.h5"
done
for d in tell_prior_disk built_apo built_lco; do
    n=$(ls "$ROOT/$d" | wc -l); nb=$(find "$ROOT/$d" -maxdepth 1 -xtype l | wc -l)
    echo "$d: $n links, $nb broken"
done
echo ASSEMBLE-PASS1C-DONE
