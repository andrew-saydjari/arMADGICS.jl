#!/usr/bin/env julia
#
# Data-model generator for an arMADGICS / APOGEE release.
#
#   julia --project=. generate.jl releases/0.4.0.toml --outdir <html dir>
#
# Options
#   --outdir DIR            where to write index.html and the product pages (required)
#   --previous FILE         snapshot to diff against (default: newest snapshots/<v>.toml with v < release)
#   --snapshot-dir DIR      where snapshots live / are written (default: datamodel/snapshots)
#   --no-snapshot           do not write snapshots/<release>.toml (dry run)
#   --accept-changes        acknowledge schema changes vs the previous release (exit 0 instead of 2)
#   --allow-undocumented    publish with UNDOCUMENTED badges (exit 0 instead of 3)
#
# Every removed, new or retyped column is printed as a banner, e.g.
#   !!!  COLUMN REMOVED: allExposure_pass1_dr21.v_helio — are we sure we don't want this column?
#   !!!  NEW COLUMN: allExposure_pass1_dr21.foo (float64) — needs a description in datamodel/registry/allExposure_pass1_dr21.toml
# and the same list is shown in a red box at the top of index.html.
#
# Exit codes: 0 ok / acknowledged, 2 unacknowledged schema changes, 3 documentation gaps.

using DataModel

function parse_cli(args)
    opts = Dict{String,Any}("accept" => false, "allow_undoc" => false, "snapshot" => true)
    pos = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--outdir"
            opts["outdir"] = args[i+1]; i += 2
        elseif a == "--previous"
            opts["previous"] = args[i+1]; i += 2
        elseif a == "--snapshot-dir"
            opts["snapshot_dir"] = args[i+1]; i += 2
        elseif a == "--no-snapshot"
            opts["snapshot"] = false; i += 1
        elseif a == "--accept-changes"
            opts["accept"] = true; i += 1
        elseif a == "--allow-undocumented"
            opts["allow_undoc"] = true; i += 1
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            push!(pos, a); i += 1
        end
    end
    length(pos) == 1 || error("usage: julia --project=. generate.jl releases/<version>.toml --outdir DIR [options]")
    haskey(opts, "outdir") || error("--outdir is required")
    return pos[1], opts
end

function main(args)
    cfg, o = parse_cli(args)
    kw = Dict{Symbol,Any}(:outdir => o["outdir"], :accept_changes => o["accept"],
                          :allow_undocumented => o["allow_undoc"], :write_snapshot_file => o["snapshot"])
    haskey(o, "previous") && (kw[:previous] = o["previous"])
    haskey(o, "snapshot_dir") && (kw[:snapshot_dir] = o["snapshot_dir"])
    code, _ = run_release(cfg; kw...)
    return code
end

exit(main(ARGS))
