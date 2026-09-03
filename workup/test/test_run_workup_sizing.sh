#!/bin/bash
# Unit tests for run_workup.sh's WORKUP_RANKS=auto arithmetic
# (sizing_report_and_choose), with mocked Slurm / node environments.
# Run:  bash workup/test/test_run_workup_sizing.sh
set -u
EP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/run_workup.sh"
fails=0
npass=0

# choose <expected_total> <expected_per_node> <case name> [ENV=VAL ...] -- <batch_bytes> <nbatch>
choose() {
    local exp_total=$1 exp_rpn=$2 name=$3; shift 3
    local envs=()
    while [ "$1" != "--" ]; do envs+=("$1"); shift; done
    shift
    local bb=$1 nb=$2
    local out
    out=$(env -i PATH="$PATH" HOME="$HOME" "${envs[@]}" WORKUP_SIZING_LIB=1 bash -c '
        source "'"$EP"'"
        sizing_report_and_choose "'"$bb"'" "'"$nb"'" >/dev/null
        echo "$AUTO_RANKS $AUTO_RANKS_PER_NODE"')
    if [ "$out" = "$exp_total $exp_rpn" ]; then
        npass=$((npass + 1))
        echo "PASS: $name -> $out"
    else
        fails=$((fails + 1))
        echo "FAIL: $name -> got '$out', expected '$exp_total $exp_rpn'"
    fi
}

# per-batch payload for all cases: the real production batch, 95,265,504 B
# → 91 MB (ceil) → per_rank = 1500 + 3×91 = 1773 MB.
BB=95265504

# 1. single-node Slurm, mem-per-node: floor(0.9*64000/1773)=32, cpus 16 → 16
choose 16 16 "slurm 1-node mem-per-node, cpu-capped" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_NODE=64000 SLURM_JOB_CPUS_PER_NODE=16 -- $BB 10000

# 2. single-node Slurm, mem-per-cpu: 2000*8=16000; floor(0.9*16000/1773)=8, cpus 8 → 8
choose 8 8 "slurm 1-node mem-per-cpu" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_CPU=2000 SLURM_JOB_CPUS_PER_NODE=8 -- $BB 10000

# 3. mem-capped below cpus: floor(0.9*4000/1773)=2 → 2
choose 2 2 "slurm mem-capped" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_NODE=4000 SLURM_JOB_CPUS_PER_NODE=16 -- $BB 10000

# 4. floor 1: floor(0.9*1000/1773)=0 → 1
choose 1 1 "slurm floor-1" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_NODE=1000 SLURM_JOB_CPUS_PER_NODE=16 -- $BB 10000

# 5. work-unit cap: plenty of resources, only 3 batches → 3
choose 3 3 "work-unit cap" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_NODE=256000 SLURM_JOB_CPUS_PER_NODE=64 -- $BB 3

# 6. MULTI-NODE: 2 nodes × (floor(0.9*64000/1773)=32, cpus 16 → 16) = 32 total, 16/node
#    (memory must NOT pool: 2×64000 MB is 128 GB but ranks/node still 16)
choose 32 16 "slurm 2-node: per-node then multiply" \
    SLURM_JOB_ID=1 SLURM_NNODES=2 SLURM_MEM_PER_NODE=64000 SLURM_JOB_CPUS_PER_NODE='16(x2)' -- $BB 10000

# 7. multi-node mem-capped: 4 nodes, floor(0.9*8000/1773)=4/node → 16 total
choose 16 4 "slurm 4-node mem-capped" \
    SLURM_JOB_ID=1 SLURM_NNODES=4 SLURM_MEM_PER_NODE=8000 SLURM_JOB_CPUS_PER_NODE='32(x4)' -- $BB 10000

# 8. multi-node + work-unit cap: 2×16=32 raw, 20 batches → 20 total, ceil(20/2)=10/node
choose 20 10 "2-node work-unit cap shrinks per-node" \
    SLURM_JOB_ID=1 SLURM_NNODES=2 SLURM_MEM_PER_NODE=64000 SLURM_JOB_CPUS_PER_NODE='16(x2)' -- $BB 20

# 9. heterogeneous cpus: "16(x2),8" → min 8/node; 3 nodes → 24
choose 24 8 "heterogeneous cpus uses min" \
    SLURM_JOB_ID=1 SLURM_NNODES=3 SLURM_MEM_PER_NODE=64000 SLURM_JOB_CPUS_PER_NODE='16(x2),8' -- $BB 10000

# 10. non-Slurm node (test hooks): 32000 MB avail, 4 cpus → min(floor(28800/1773)=16, 4)=4
choose 4 4 "local node, cpu-capped" \
    WORKUP_TEST_MEMAVAIL_MB=32000 WORKUP_TEST_NPROC=4 -- $BB 10000

# 11. WORKUP_MEM_FRACTION override: 0.5×64000/1773 → 18, cpus 64 → 18
choose 18 18 "mem-fraction override" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_NODE=64000 SLURM_JOB_CPUS_PER_NODE=64 WORKUP_MEM_FRACTION=0.5 -- $BB 10000

# 12. fat-key future: 1 GB batches → per_rank = 1500+3×1024 = 4572 MB;
#     floor(0.9*64000/4572)=12, cpus 64 → 12 (RAM-capping protects)
choose 12 12 "fat batches RAM-capped" \
    SLURM_JOB_ID=1 SLURM_NNODES=1 SLURM_MEM_PER_NODE=64000 SLURM_JOB_CPUS_PER_NODE=64 -- 1073741824 10000

echo
if [ "$fails" -eq 0 ]; then
    echo "All $npass sizing tests passed."
    exit 0
else
    echo "$fails sizing test(s) FAILED ($npass passed)."
    exit 1
fi
