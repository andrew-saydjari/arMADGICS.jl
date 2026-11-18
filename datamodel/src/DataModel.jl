"""
    DataModel

Data-model generator for the arMADGICS / APOGEE DR21 release products.

Each run:

1. introspects the files named in a release config (`releases/<version>.toml`),
2. writes a schema snapshot (`snapshots/<version>.toml`: product -> column -> dtype/shape),
3. diffs that snapshot against the previous release's snapshot and against the curated
   column registry (`registry/<product>.toml`), printing a loud banner for every removed,
   new or changed column and for every column with no description,
4. renders static HTML (index + one page per product).

See `generate.jl` for the command-line entry point.
"""
module DataModel

using CSV
using Dates
using HDF5
using Printf
using TOML

export ColumnSchema, ProductSchema, Change
export introspect_product, load_release, load_registry, snapshot_dict, read_snapshot,
       write_snapshot, diff_snapshots, registry_check, resolve_column, render_site,
       find_previous_snapshot, report_changes, run_release

include("schema.jl")
include("introspect.jl")
include("registry.jl")
include("harvest.jl")
include("diff.jl")
include("render.jl")
include("run.jl")

end # module
