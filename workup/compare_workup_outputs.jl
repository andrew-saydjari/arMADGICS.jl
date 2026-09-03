## compare_workup_outputs.jl — dataset-by-dataset bitwise comparison of two
## workup output directories (h5diff mindset, but NaN-aware and streaming).
#
# For every arMADGICS_out_<key>.h5 in --a: demand the same file in --b, the
# same key-dataset shape and eltype, and elementwise `isequal` (bitwise for
# floats up to NaN payload; NaN == NaN) over row-chunks. Also compares the
# `hdr` dataset attrs (producer git provenance). Extra non-key datasets in
# either side (e.g. missing_row_mask, workup file attrs) are listed, not
# compared. Exit status 0 iff everything matches.
#
# Usage:
#   nice -n 10 julia --project=. compare_workup_outputs.jl --a DIR_A --b DIR_B [--chunk 20000]

using ArgParse
using HDF5
using Printf

function parse_commandline()
    s = ArgParseSettings(description = "Bitwise comparison of two workup output dirs")
    @add_arg_table s begin
        "--a"
        arg_type = String
        required = true
        "--b"
        arg_type = String
        required = true
        "--chunk"
        help = "rows per comparison chunk"
        arg_type = Int
        default = 20000
    end
    return parse_args(s)
end

const parg = parse_commandline()
const dira = abspath(expanduser(parg["a"]))
const dirb = abspath(expanduser(parg["b"]))
const CHUNK = parg["chunk"]

const OUT_RE = r"^arMADGICS_out_(.+)\.h5$"

function keyfiles(dir)
    out = Dict{String, String}()
    for fn in sort(readdir(dir))
        m = match(OUT_RE, fn)
        isnothing(m) && continue
        out[String(m.captures[1])] = joinpath(dir, fn)
    end
    return out
end

"Streaming elementwise isequal over the last (row) axis; returns (ok, first_bad_row)."
function datasets_equal(da, db; chunk = CHUNK)
    n = size(da)[end]
    nd = ndims(da)
    idx = ntuple(_ -> Colon(), nd - 1)
    i = 1
    while i <= n
        j = min(i + chunk - 1, n)
        A = da[idx..., i:j]
        B = db[idx..., i:j]
        if !isequal(A, B)
            # locate the first differing row for the report
            for (kk, r) in enumerate(i:j)
                av = selectdim(A, nd, kk)
                bv = selectdim(B, nd, kk)
                isequal(collect(av), collect(bv)) || return (false, r)
            end
            return (false, i)   # unreachable, defensive
        end
        i = j + 1
    end
    return (true, 0)
end

fa = keyfiles(dira)
fb = keyfiles(dirb)

only_a = sort(collect(setdiff(keys(fa), keys(fb))))
only_b = sort(collect(setdiff(keys(fb), keys(fa))))
shared = sort(collect(intersect(keys(fa), keys(fb))))

println("Comparing $(length(shared)) shared keys")
isempty(only_a) || println("  keys only in A: $(join(only_a, ", "))")
isempty(only_b) || println("  keys only in B: $(join(only_b, ", "))")

nfail = 0
hdr_mismatch = 0
for k in shared
    global nfail, hdr_mismatch
    ha = h5open(fa[k], "r")
    hb = h5open(fb[k], "r")
    try
        status = "PASS"
        detail = ""
        if !haskey(ha, k) || !haskey(hb, k)
            status = "FAIL"; detail = "key dataset missing"
        else
            da, db = ha[k], hb[k]
            if size(da) != size(db)
                status = "FAIL"; detail = "shape $(size(da)) vs $(size(db))"
            elseif eltype(da) != eltype(db)
                status = "FAIL"; detail = "eltype $(eltype(da)) vs $(eltype(db))"
            else
                ok, badrow = datasets_equal(da, db)
                ok || (status = "FAIL"; detail = "first differing row: $badrow")
            end
        end
        # hdr provenance attrs
        hstat = ""
        if haskey(ha, "hdr") && haskey(hb, "hdr")
            aat = Dict(a => read_attribute(ha["hdr"], a) for a in keys(attributes(ha["hdr"])))
            bat = Dict(a => read_attribute(hb["hdr"], a) for a in keys(attributes(hb["hdr"])))
            if !isequal(aat, bat)
                hstat = " [hdr attrs DIFFER]"
                hdr_mismatch += 1
            end
        else
            hstat = " [hdr dataset missing on one side]"
            hdr_mismatch += 1
        end
        status == "FAIL" && (nfail += 1)
        # note extra datasets (informational)
        extra_a = sort([n for n in keys(ha) if n != k && n != "hdr"])
        extra_b = sort([n for n in keys(hb) if n != k && n != "hdr"])
        extras = (isempty(extra_a) && isempty(extra_b)) ? "" :
            " (extra dsets A: [$(join(extra_a, ","))] B: [$(join(extra_b, ","))])"
        @printf("%-32s %s%s%s%s\n", k, status,
            isempty(detail) ? "" : " — " * detail, hstat, extras)
    finally
        close(ha); close(hb)
    end
end

println()
if nfail == 0 && hdr_mismatch == 0 && isempty(only_a) && isempty(only_b)
    println("VERDICT: IDENTICAL — all $(length(shared)) key datasets bitwise-equal (isequal, NaN-aware), hdr provenance attrs equal.")
    exit(0)
else
    println("VERDICT: DIFFER — $nfail key dataset failures, $hdr_mismatch hdr attr mismatches, $(length(only_a)) keys only in A, $(length(only_b)) keys only in B.")
    exit(1)
end
