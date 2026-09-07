"""Split the full input lists into single-telescope contiguous shards.

Each shard is a contiguous slice of list_<tele>_full.txt (which is sorted by
(mjd, exposure)), so shard outputs concatenate back into full-list order and
every shard output is single-telescope (matching the delivered per-telescope
layout consumed downstream by sample_starCont.jl).

Defaults: 5 APO + 3 LCO shards = 8 GPU tasks, sized so the largest shard is
~2360 files (~14.3 h at the validated 21.8 s/file A6000 rate; A100/H100 f64
throughput is higher, so this is an upper bound).

Writes shards/list_shard_<k>.txt and shards/shard_map.tsv
(columns: shard_id  tele  n_files  list_path  output_path).
"""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--n-apo", type=int, default=5)
parser.add_argument("--n-lco", type=int, default=3)
args = parser.parse_args()

base = Path(__file__).resolve().parent
shdir = base / "shards"
shdir.mkdir(exist_ok=True)

plan = []  # (tele, paths_slice)
for tele, n_shards in (("apo", args.n_apo), ("lco", args.n_lco)):
    paths = [l for l in (base / f"list_{tele}_full.txt").read_text().splitlines() if l.strip()]
    n = len(paths)
    for k in range(n_shards):
        lo = k * n // n_shards
        hi = (k + 1) * n // n_shards
        plan.append((tele, paths[lo:hi]))

rows = []
for sid, (tele, chunk) in enumerate(plan):
    lp = shdir / f"list_shard_{sid}.txt"
    lp.write_text("\n".join(chunk) + "\n")
    op = base / "out" / f"shard_{sid}_{tele}.h5"
    rows.append(f"{sid}\t{tele}\t{len(chunk)}\t{lp}\t{op}")
    print(rows[-1])

(shdir / "shard_map.tsv").write_text("\n".join(rows) + "\n")
print(f"wrote {len(plan)} shards; max shard size = {max(len(c) for _, c in plan)} files")
