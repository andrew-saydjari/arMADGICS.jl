#!/bin/bash
# pass-1b corpus + built-prior view under the ADOPTED intelligent bright-cut
# policy (gap-exposure screen + APO 10k ceiling). Supersedes starCont_pass1
# (C2_LCO=3000) for pass-1. Everything referenced is retained; nothing deleted.
set -euo pipefail
APO_S=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/prior_outputs/starCont_20260903/tell_prior_disk
LCO_S=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_20260904_lcoC2/tell_prior_disk
NEW_S=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_intelligent/tell_prior_disk
OLD_BA=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_apo
OLD_BL=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_lco
NEW_B=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_intelligent/built
ROOT=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_pass1b
CHANGED="388 448 459 519"
mkdir -p "$ROOT/tell_prior_disk" "$ROOT/built_apo" "$ROOT/built_lco"

is_changed() { for c in $CHANGED; do [ "$1" -eq "$c" ] && return 0; done; return 1; }

for i in $(seq 1 600); do
    i3=$(printf "%03d" "$i")
    if [ "$i" -le 300 ]; then src=$APO_S
    elif is_changed "$i"; then src=$NEW_S
    else src=$LCO_S; fi
    ln -sfn "$src/starCont_$i3.jdat" "$ROOT/tell_prior_disk/starCont_$i3.jdat"
done
for i in $(seq 1 300); do
    i3=$(printf "%03d" "$i")
    ln -sfn "$OLD_BA/APOGEE_starcont_svd_60_f$i3.h5" "$ROOT/built_apo/APOGEE_starcont_svd_60_f$i3.h5"
done
for i in $(seq 301 600); do
    i3=$(printf "%03d" "$i")
    if is_changed "$i"; then src=$NEW_B; else src=$OLD_BL; fi
    ln -sfn "$src/APOGEE_starcont_svd_60_f$i3.h5" "$ROOT/built_lco/APOGEE_starcont_svd_60_f$i3.h5"
done
for d in tell_prior_disk built_apo built_lco; do
    n=$(ls "$ROOT/$d" | wc -l); nb=$(find "$ROOT/$d" -maxdepth 1 -xtype l | wc -l)
    echo "$d: $n links, $nb broken"
done
echo ASSEMBLE-PASS1B-DONE
