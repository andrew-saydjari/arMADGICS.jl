## validate_workup.jl — W3: standalone row-matching validator for the
## arMADGICS batch corpus (and, by extension, any workup built on it).
#
# For K sampled rows plus ALL rows at fiber/batch boundaries (first and last
# row of every present batch), takes (tele, mjd, expnum, adjfiberindx) from
# full_list_info.h5 — the identity source — re-reads the corresponding
# ar1Dunical_<tele>_<mjd>_<expnum>_object.h5 fiber row, and demands EXACT
# (bitwise, NaN-aware) equality with the batch file's stored flux/fluxivar.
#
# The producer stores flux = nanify(fspec[simplemsk], simplemsk) (pipeline.jl
# ~line 246): the ar1Dunical values verbatim on simplemsk pixels, NaN
# elsewhere, with simplemsk saved in the same batch. So the exactness demand
# is: flux[simplemsk] == flux_1d[:,fiber][simplemsk] bitwise AND
# flux[!simplemsk] all-NaN — a free end-to-end identity check through the
# whole chain. Any row shift (the workup "matching issue") breaks it loudly.
#
# Also cross-checks batch_info.txt vs full_list_info.h5 row-for-row, verifies
# the contiguous-slice ordering contract, reconciles the discovered batch set
# against the identity-derived expected set, and compares nsamp vs the sum of
# actual batch sizes.
#
# Usage (example, partial 2026_05_01 corpus):
#   nice -n 10 julia --project=. validate_workup.jl \
#       --rawdir /mnt/ceph/.../outdir/arMADGICS/raw_57600_61160 \
#       --redux  /mnt/ceph/.../outdir \
#       --K 500 --allow-missing \
#       --out /path/to/W3_report.md

import Pkg
using Dates
using Distributed
using Printf
using Random

using ArgParse

function parse_commandline()
    s = ArgParseSettings(description = "W3 row-matching validator for arMADGICS batch corpora")
    @add_arg_table s begin
        "--rawdir"
        help = "arMADGICS raw batch dir (contains NNN/ fiber dirs, batch_info.txt, full_list_info.h5)"
        arg_type = String
        required = true
        "--redux"
        help = "reduxBase (the outdir containing apred/<mjd>/ar1Dunical_*.h5)"
        arg_type = String
        required = true
        "--out"
        help = "markdown report output path"
        arg_type = String
        required = true
        "--failcsv"
        help = "optional CSV path for per-row failure records"
        arg_type = String
        default = ""
        "--K"
        help = "number of random interior rows to sample (in addition to all boundary rows)"
        arg_type = Int
        default = 500
        "--seed"
        help = "RNG seed for row sampling"
        arg_type = Int
        default = 2026
        "--batchsize"
        help = "producer batch size"
        arg_type = Int
        default = 100
        "--allow-missing"
        help = "tolerate missing batch files (partial corpus); missing rows are reported, never silently skipped"
        action = :store_true
        "--nworkers"
        help = "worker processes (validation is ceph-I/O-bound; more than ~16 rarely helps)"
        arg_type = Int
        default = 3
        "--max-batches"
        help = "limit number of present batches to row-check (0 = all; for smoke tests)"
        arg_type = Int
        default = 0
        "--skip-crosscheck"
        help = "skip the 26M-row batch_info.txt crosscheck (for smoke tests only)"
        action = :store_true
        "--chunk"
        help = "batches per checkpointed work chunk"
        arg_type = Int
        default = 2000
        "--resume"
        help = "resume row checks from <out>.ckpt if present (corpus-level checks rerun; row selection is deterministic given the same seed/K)"
        action = :store_true
    end
    return parse_args(s)
end

const parg = parse_commandline()
const t0 = now()

nw = parg["nworkers"]  # no hard cap: the ≤4-proc limit during initial validation was a
                       # shared-node scheduling constraint of that session, not a code property
nw > 0 && addprocs(nw; exeflags = "--project=$(Base.active_project())")

const RCPATH = joinpath(@__DIR__, "RowContract.jl")
@everywhere using HDF5
@everywhere include($RCPATH)
@everywhere using .RowContract

@everywhere begin
    """
    One work unit: a present batch plus the rows in it to check.
    Row spec: (within, globalrow, telecode, mjd, expnum, category).
    """
    struct BatchWork
        id::BatchId
        path::String
        nrow_expect::Int
        rows::Vector{NTuple{6, Int}}
    end

    struct RowResult
        globalrow::Int
        adjfiberindx::Int
        startind::Int
        within::Int
        telecode::Int
        mjd::Int
        expnum::Int
        category::Int          # 1=boundary, 2=sampled
        npix_msk::Int
        flux_ok::Bool
        ivar_ok::Bool
        fiber_ok::Bool
        cnt_ok::Bool
        msk_subset_ok::Bool
        input_ok::Bool         # ar1Dunical readable
        shift_info::String     # non-empty only on flux mismatch
    end

    struct BatchResult
        id::BatchId
        openable::Bool
        integrity::Vector{String}
        nrow_actual::Int
        rows::Vector{RowResult}
    end

    function get_1Duni_name(reduxBase, tele, mjd, expnum)
        return joinpath(reduxBase, "apred", string(mjd),
            join(["ar1Dunical", tele, string(mjd), lpad(expnum, 4, "0"), "object.h5"], "_"))
    end

    adjfiberindx2fiberindx(a::Int) = a > 300 ? a - 300 : a

    "Bitwise, NaN-aware equality of the masked copy against the source row."
    function masked_copy_matches(copyvec, srcvec, msk)
        @inbounds for p in eachindex(copyvec)
            if msk[p]
                isequal(copyvec[p], srcvec[p]) || return false
            else
                isnan(copyvec[p]) || return false
            end
        end
        return true
    end

    """
    On a flux mismatch, scan all 300 fibers of the ar1Dunical file for
    columns whose masked pixels DO reproduce the batch copy — direct
    evidence of a row/fiber shift (vs a value-level discrepancy).
    """
    function shift_diagnostic(fluxrow, msk, flux_all, expect_fib)
        hits = Int[]
        nfib = size(flux_all, 2)
        any(msk) || return "no unmasked pixels to compare"
        for j in 1:nfib
            ok = true
            @inbounds for p in eachindex(msk)
                if msk[p] && !isequal(fluxrow[p], flux_all[p, j])
                    ok = false
                    break
                end
            end
            ok && push!(hits, j)
        end
        if isempty(hits)
            return "matches NO fiber column of the identified ar1Dunical file (value-level mismatch or wrong exposure)"
        else
            return "matches fiber column(s) $(hits) instead of expected fiber $(expect_fib) — fiber/row SHIFT"
        end
    end

    function process_batch(bw::BatchWork; ref_keyinfo = nothing, reduxBase = "")
        rows = RowResult[]
        integrity = String[]
        nrow_actual = -1
        local fb
        try
            fb = h5open(bw.path, "r")
        catch e
            return BatchResult(bw.id, false, ["unreadable: $(sprint(showerror, e))"], -1, rows)
        end
        try
            # integrity: every dataset's last axis == identity-derived length
            for k in keys(fb)
                k == "hdr" && continue
                obj = fb[k]
                obj isa HDF5.Dataset || continue
                sz = size(obj)
                nlast = isempty(sz) ? -1 : sz[end]
                if nlast != bw.nrow_expect
                    push!(integrity, "key $k: last-axis $nlast != expected $(bw.nrow_expect)")
                end
                if !isnothing(ref_keyinfo) && haskey(ref_keyinfo, k)
                    rv = ref_keyinfo[k]
                    (eltype(obj) == rv.eltype && sz[1:end-1] == rv.shape[1:end-1]) ||
                        push!(integrity, "key $k: shape/type $(sz)/$(eltype(obj)) inconsistent with reference")
                end
            end
            haskey(fb, "hdr") || push!(integrity, "missing hdr dataset")
            # reference keys absent from this file (e.g. batch truncated by a
            # mid-write kill) are as important as wrong shapes
            if !isnothing(ref_keyinfo)
                for k in keys(ref_keyinfo)
                    haskey(fb, k) || push!(integrity, "missing key $k")
                end
            end
            nrow_actual = haskey(fb, "adjfiberindx") ? length(fb["adjfiberindx"]) : -1

            d_flux = fb["flux"]
            d_ivar = fb["fluxivar"]
            d_msk = fb["simplemsk"]
            d_fib = fb["adjfiberindx"]
            d_cnt = fb["data_pix_cnt"]

            # group rows by (tele,mjd,expnum) so each ar1Dunical is opened once
            byexp = Dict{Tuple{Int, Int, Int}, Vector{NTuple{6, Int}}}()
            for r in bw.rows
                push!(get!(byexp, (r[3], r[4], r[5]), NTuple{6, Int}[]), r)
            end

            fibindx = adjfiberindx2fiberindx(bw.id.adjfiberindx)
            for ((tc, mjd, expnum), rs) in byexp
                tele = tele_name(tc)
                arpath = get_1Duni_name(reduxBase, tele, mjd, expnum)
                input_ok = true
                fspec = fivar = fmsk = nothing
                far = nothing
                try
                    # column hyperslabs only: reading the full {300,8700}
                    # arrays is ~300x more I/O and only needed for the
                    # shift diagnostic on a mismatch
                    far = h5open(arpath, "r")
                    fspec = far["flux_1d"][:, fibindx]
                    fivar = far["ivar_1d"][:, fibindx]
                    fmsk = far["mask_1d"][:, fibindx]
                catch e
                    input_ok = false
                end
                for r in rs
                    within, grow, _, _, _, cat = r[1], r[2], r[3], r[4], r[5], r[6]
                    fluxrow = d_flux[:, within]
                    ivarrow = d_ivar[:, within]
                    mskrow = d_msk[:, within] .== 1
                    fibrow = Int(d_fib[within])
                    cntrow = Int(d_cnt[within])
                    fiber_ok = (fibrow == bw.id.adjfiberindx)
                    cnt_ok = (cntrow == count(mskrow))
                    if !input_ok
                        push!(rows, RowResult(grow, bw.id.adjfiberindx, bw.id.startind,
                            within, tc, mjd, expnum, cat, count(mskrow),
                            false, false, fiber_ok, cnt_ok, false, false,
                            "ar1Dunical unreadable: $arpath"))
                        continue
                    end
                    flux_ok = masked_copy_matches(fluxrow, fspec, mskrow)
                    ivar_ok = masked_copy_matches(ivarrow, fivar, mskrow)
                    msk_subset_ok = !any(mskrow .& .!Bool.(fmsk))
                    shift_info = ""
                    if !flux_ok
                        # only now pay for the full flux array
                        flux_all = try
                            read(far["flux_1d"])
                        catch
                            nothing
                        end
                        shift_info = isnothing(flux_all) ?
                            "flux mismatch (full-file re-read failed)" :
                            shift_diagnostic(fluxrow, mskrow, flux_all, fibindx)
                    end
                    push!(rows, RowResult(grow, bw.id.adjfiberindx, bw.id.startind,
                        within, tc, mjd, expnum, cat, count(mskrow),
                        flux_ok, ivar_ok, fiber_ok, cnt_ok, msk_subset_ok, true,
                        shift_info))
                end
                isnothing(far) || try close(far) catch end
            end
        catch e
            push!(integrity, "error during row checks: $(sprint(showerror, e))")
        finally
            close(fb)
        end
        return BatchResult(bw.id, true, integrity, nrow_actual, rows)
    end
end

# ============================================================== main (master)

using .RowContract

rawdir = abspath(expanduser(parg["rawdir"]))
reduxBase = abspath(expanduser(parg["redux"]))
outpath = abspath(expanduser(parg["out"]))
mkpath(dirname(outpath))

logmsg(s) = (println("[$(now())] $s"); flush(stdout))

logmsg("W3 validator starting. rawdir=$rawdir redux=$reduxBase")

# ---- 1. identity load + ordering contract ----------------------------------
logmsg("Loading full_list_info.h5 (identity source)...")
fli = load_full_list_info(joinpath(rawdir, "full_list_info.h5"))
logmsg("  nsamp = $(fli.nsamp)")
fidx = verify_fiber_blocks(fli; batchsize = parg["batchsize"])
n_fibers_nonempty = count(>(0), fidx.counts)
logmsg("Ordering contract VERIFIED: $(n_fibers_nonempty) non-empty fibers in contiguous increasing blocks.")

# ---- 2. batch_info.txt crosscheck ------------------------------------------
cc = nothing
if !parg["skip-crosscheck"]
    logmsg("Cross-checking batch_info.txt row-for-row against full_list_info.h5...")
    cc = crosscheck_batch_info(joinpath(rawdir, "batch_info.txt"), fli, fidx)
    logmsg("  rows=$(cc.nrows) mismatches=$(cc.n_mismatch) header_total_batches=$(cc.header_total_batches)")
end

# ---- 3. discovery + reconciliation ------------------------------------------
logmsg("Discovering batch files...")
discovered = discover_batches(rawdir)
expected = expected_batches(fidx)
logmsg("  discovered=$(length(discovered)) expected(full corpus)=$(length(expected))")
rec = reconcile_batches(discovered, expected; allow_missing = parg["allow-missing"])
logmsg("  present=$(length(rec.present)) missing=$(length(rec.missing)) extra=$(length(rec.extra))")

# per-fiber completeness for the report
fiber_present = zeros(Int, 600)
fiber_expected = zeros(Int, 600)
for id in expected
    fiber_expected[id.adjfiberindx] += 1
end
for id in rec.present
    fiber_present[id.adjfiberindx] += 1
end

# rows covered by present batches
present_row_count = sum(length(batch_row_range(fidx, id)) for id in rec.present)
missing_msk = missing_rows_mask(fidx, rec.missing)
@assert present_row_count + count(missing_msk) == fli.nsamp

# ---- 4. reference key discovery + hdr provenance ---------------------------
refpath = rec.paths[rec.present[1]]
ref_keyinfo = discover_keys(refpath)
hdr = read_hdr_attrs(refpath)
logmsg("Reference batch $(batch_filename(rec.present[1])): $(length(ref_keyinfo)) keys; hdr attrs: $(join(sort(collect(keys(hdr))), ", "))")

# ---- 5. row selection ------------------------------------------------------
check_batches = rec.present
if parg["max-batches"] > 0 && length(check_batches) > parg["max-batches"]
    check_batches = check_batches[1:parg["max-batches"]]
    logmsg("SMOKE MODE: limiting row checks to first $(length(check_batches)) present batches")
end

rng = MersenneTwister(parg["seed"])
rows_by_batch = Dict{BatchId, Vector{NTuple{6, Int}}}()
function add_row!(id::BatchId, within::Int, cat::Int)
    wr = batch_within_range(fidx, id)
    grow = fidx.offsets[id.adjfiberindx] + id.startind + within - 1
    spec = (within, grow, Int(fli.tele[grow]), Int(fli.mjd[grow]),
        Int(fli.expnum[grow]), cat)
    v = get!(rows_by_batch, id, NTuple{6, Int}[])
    any(x -> x[1] == within, v) && return
    push!(v, spec)
end

# all fiber/batch boundary rows: first + last row of every present batch
for id in check_batches
    wr = batch_within_range(fidx, id)
    add_row!(id, 1, 1)
    add_row!(id, length(wr), 1)
end
# K random interior rows over rows covered by the checked batches
check_set = Set(check_batches)
global_rows_avail = Int[]
for id in check_batches
    append!(global_rows_avail, batch_row_range(fidx, id))
end
for grow in rand(rng, global_rows_avail, parg["K"])
    fib = Int(fli.adjfiberindx[grow])
    li = grow - fidx.offsets[fib]
    id = batch_of_linear_index(fidx, fib, li)
    add_row!(id, li - id.startind + 1, 2)
end

work = [BatchWork(id, rec.paths[id], length(batch_within_range(fidx, id)),
            sort(rows_by_batch[id]))
        for id in sort(collect(keys(rows_by_batch)))]
n_rows_to_check = sum(length(w.rows) for w in work)
logmsg("Row checks queued: $n_rows_to_check rows across $(length(work)) batches")

# ---- 6. distributed row checks (chunked, checkpointed) ----------------------
using Serialization
ckptpath = outpath * ".ckpt"
results = BatchResult[]
done_ids = Set{BatchId}()
if parg["resume"] && isfile(ckptpath)
    ck = open(deserialize, ckptpath)
    results = ck.results
    done_ids = Set(br.id for br in results)
    logmsg("Resumed checkpoint: $(length(done_ids)) batches already checked")
end
todo = [w for w in work if !(w.id in done_ids)]
logmsg("Running distributed checks on $(nworkers()) workers ($(length(todo)) batches to do)...")
process_one = bw -> process_batch(bw; ref_keyinfo = ref_keyinfo, reduxBase = reduxBase)

function run_chunks!(results, todo, ntotal, nchunk, ckptpath)
    i = 1
    while i <= length(todo)
        j = min(i + nchunk - 1, length(todo))
        tchunk0 = time()
        append!(results, pmap(process_one, todo[i:j]; batch_size = 25))
        # atomic checkpoint: write tmp then rename
        tmp = ckptpath * ".tmp"
        open(io -> serialize(io, (results = results,)), tmp, "w")
        mv(tmp, ckptpath; force = true)
        rate = (j - i + 1) / (time() - tchunk0)
        eta = round(Int, (length(todo) - j) / max(rate, 1e-9))
        logmsg(@sprintf("  checked %d / %d batches (%.1f batches/s, ETA %ds)",
            length(results), ntotal, rate, eta))
        i = j + 1
    end
end
run_chunks!(results, todo, length(work), max(parg["chunk"], 1), ckptpath)
logmsg("Row checks complete.")

# ---- 7. aggregate ----------------------------------------------------------
all_rows = RowResult[]
unopenable = String[]
integrity_bad = Tuple{BatchId, Vector{String}}[]
sum_actual_rows = 0
for br in results
    if !br.openable
        push!(unopenable, batch_filename(br.id))
    else
        global sum_actual_rows += max(br.nrow_actual, 0)
        isempty(br.integrity) || push!(integrity_bad, (br.id, br.integrity))
    end
    append!(all_rows, br.rows)
end
sum_expected_present = sum(length(batch_within_range(fidx, id)) for id in check_batches)

nrows = length(all_rows)
c(f) = count(f, all_rows)
n_input_bad = c(r -> !r.input_ok)
n_flux_ok = c(r -> r.flux_ok)
n_ivar_ok = c(r -> r.ivar_ok)
n_fiber_ok = c(r -> r.fiber_ok)
n_cnt_ok = c(r -> r.cnt_ok)
n_msub_ok = c(r -> r.msk_subset_ok)
n_all_ok = c(r -> r.flux_ok && r.ivar_ok && r.fiber_ok && r.cnt_ok && r.msk_subset_ok && r.input_ok)
n_vacuous = c(r -> r.npix_msk == 0)
n_boundary = c(r -> r.category == 1)
n_sampled = c(r -> r.category == 2)
failures = filter(r -> !(r.flux_ok && r.ivar_ok && r.fiber_ok && r.cnt_ok && r.msk_subset_ok && r.input_ok), all_rows)
shift_rows = filter(r -> occursin("SHIFT", r.shift_info), all_rows)

# optional failure CSV
if !isempty(parg["failcsv"]) && !isempty(failures)
    mkpath(dirname(abspath(parg["failcsv"])))
    open(parg["failcsv"], "w") do io
        println(io, "globalrow,adjfiberindx,startind,within,tele,mjd,expnum,category,npix_msk,flux_ok,ivar_ok,fiber_ok,cnt_ok,msk_subset_ok,input_ok,shift_info")
        for r in failures
            println(io, join([r.globalrow, r.adjfiberindx, r.startind, r.within,
                tele_name(r.telecode), r.mjd, r.expnum, r.category, r.npix_msk,
                r.flux_ok, r.ivar_ok, r.fiber_ok, r.cnt_ok, r.msk_subset_ok,
                r.input_ok, "\"$(r.shift_info)\""], ","))
        end
    end
end

# ---- 8. report -------------------------------------------------------------
pct(a, b) = b == 0 ? "n/a" : @sprintf("%.4f%%", 100a / b)
verdict_rows = (n_all_ok == nrows) ? "PASS" : ((n_input_bad == nrows - n_all_ok) ? "PASS (with unreadable inputs)" : "FAIL")
verdict_completeness = isempty(rec.missing) ? "COMPLETE" : "PARTIAL ($(length(rec.present))/$(length(expected)) batches present)"

open(outpath, "w") do io
    println(io, "# W3 row-matching validation report")
    println(io, "")
    println(io, "- Generated: $(now()) on $(gethostname())")
    println(io, "- rawdir: `$rawdir`")
    println(io, "- redux: `$reduxBase`")
    println(io, "- K (random rows requested): $(parg["K"]), seed $(parg["seed"]); plus first+last row of every checked batch")
    println(io, "- Batch producer provenance (from `hdr` dataset attrs of $(batch_filename(rec.present[1]))):")
    for k in sort(collect(keys(hdr)))
        println(io, "    - $k: `$(hdr[k])`")
    end
    println(io, "")

    println(io, "## Verdict summary")
    println(io, "")
    println(io, "| check | result |")
    println(io, "|---|---|")
    println(io, "| ordering contract (contiguous per-fiber blocks) | VERIFIED |")
    if !isnothing(cc)
        println(io, "| batch_info.txt vs full_list_info.h5 row-for-row | $(cc.n_mismatch == 0 ? "PASS ($(cc.nrows) rows)" : "FAIL ($(cc.n_mismatch) mismatches)") |")
        println(io, "| header total batches vs identity-derived | $(cc.header_total_batches == length(expected) ? "PASS ($(length(expected)))" : "MISMATCH (header $(cc.header_total_batches) vs derived $(length(expected)))") |")
    end
    println(io, "| batch-set completeness | $verdict_completeness |")
    println(io, "| extra (unexpected) batch files | $(isempty(rec.extra) ? "none" : "$(length(rec.extra)) — CONTRACT VIOLATION") |")
    println(io, "| nsamp vs sum of batch sizes (checked batches) | $(sum_actual_rows == sum_expected_present ? "PASS ($sum_actual_rows)" : "FAIL (actual $sum_actual_rows vs expected $sum_expected_present)") |")
    println(io, "| per-batch integrity (keys/shapes) | $(isempty(integrity_bad) && isempty(unopenable) ? "PASS ($(length(results)) batches)" : "FAIL ($(length(integrity_bad)) bad, $(length(unopenable)) unopenable)") |")
    println(io, "| row-matching (exact flux/fluxivar identity) | $verdict_rows |")
    println(io, "| row-shift evidence | $(isempty(shift_rows) ? "NONE" : "$(length(shift_rows)) rows — SEE BELOW") |")
    println(io, "")

    println(io, "## Corpus accounting")
    println(io, "")
    println(io, "- nsamp (full_list_info.h5 rows): $(fli.nsamp)")
    println(io, "- non-empty fibers (adjfiberindx): $n_fibers_nonempty")
    println(io, "- expected batches (full corpus, batchsize $(parg["batchsize"])): $(length(expected))")
    println(io, "- discovered batch files: $(length(discovered))")
    println(io, "- present (expected AND on disk): $(length(rec.present)) = $(pct(length(rec.present), length(expected))) of expected")
    println(io, "- missing batches: $(length(rec.missing)) covering $(count(missing_msk)) rows ($(pct(count(missing_msk), fli.nsamp)) of nsamp)")
    println(io, "- rows covered by present batches: $present_row_count")
    println(io, "")

    println(io, "### Per-fiber completeness (fibers with any expected batches)")
    println(io, "")
    lastfib_any = findlast(>(0), fiber_present)
    println(io, "- fibers with ALL batches present: $(count(f -> fiber_expected[f] > 0 && fiber_present[f] == fiber_expected[f], 1:600))")
    println(io, "- fibers with SOME batches present: $(count(f -> 0 < fiber_present[f] < fiber_expected[f], 1:600))")
    println(io, "- fibers with NO batches present: $(count(f -> fiber_expected[f] > 0 && fiber_present[f] == 0, 1:600))")
    println(io, "- highest fiber with any batch present: $(isnothing(lastfib_any) ? "none" : lastfib_any)")
    incomplete = [f for f in 1:600 if 0 < fiber_present[f] < fiber_expected[f]]
    if !isempty(incomplete)
        println(io, "")
        println(io, "| adjfiberindx | present | expected |")
        println(io, "|---|---|---|")
        for f in incomplete
            println(io, "| $f | $(fiber_present[f]) | $(fiber_expected[f]) |")
        end
    end
    println(io, "")

    if !isnothing(cc) && cc.n_mismatch > 0
        println(io, "### batch_info.txt mismatches (first $(length(cc.first_mismatches)))")
        println(io, "")
        for l in cc.first_mismatches
            println(io, "- $l")
        end
        println(io, "")
    end

    if !isempty(unopenable) || !isempty(integrity_bad)
        println(io, "## Batch integrity problems")
        println(io, "")
        for b in unopenable
            println(io, "- UNOPENABLE: $b")
        end
        for (id, probs) in integrity_bad
            println(io, "- $(batch_filename(id)): $(join(probs, "; "))")
        end
        println(io, "")
    end

    println(io, "## Row-matching results")
    println(io, "")
    println(io, "Checked $nrows rows ($(n_boundary) boundary, $(n_sampled) random-sampled) across $(length(work)) batches.")
    println(io, "")
    println(io, "| check | pass | fail | pass rate |")
    println(io, "|---|---|---|---|")
    println(io, "| flux exact match | $n_flux_ok | $(nrows - n_flux_ok) | $(pct(n_flux_ok, nrows)) |")
    println(io, "| fluxivar exact match | $n_ivar_ok | $(nrows - n_ivar_ok) | $(pct(n_ivar_ok, nrows)) |")
    println(io, "| stored adjfiberindx correct | $n_fiber_ok | $(nrows - n_fiber_ok) | $(pct(n_fiber_ok, nrows)) |")
    println(io, "| data_pix_cnt == count(simplemsk) | $n_cnt_ok | $(nrows - n_cnt_ok) | $(pct(n_cnt_ok, nrows)) |")
    println(io, "| simplemsk subset of ar1Dunical good mask | $n_msub_ok | $(nrows - n_msub_ok) | $(pct(n_msub_ok, nrows)) |")
    println(io, "| ar1Dunical input readable | $(nrows - n_input_bad) | $n_input_bad | $(pct(nrows - n_input_bad, nrows)) |")
    println(io, "| ALL checks | $n_all_ok | $(nrows - n_all_ok) | $(pct(n_all_ok, nrows)) |")
    println(io, "")
    println(io, "- rows with zero unmasked pixels (vacuous flux comparison): $n_vacuous")
    println(io, "")

    if !isempty(failures)
        println(io, "### Failed rows (full identity tuples)")
        println(io, "")
        println(io, "| globalrow | tele | mjd | expnum | adjfiberindx | batch | within | fails | shift diagnostic |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        for r in first(failures, 200)
            fails = String[]
            r.flux_ok || push!(fails, "flux")
            r.ivar_ok || push!(fails, "ivar")
            r.fiber_ok || push!(fails, "fiberid")
            r.cnt_ok || push!(fails, "cnt")
            r.msk_subset_ok || push!(fails, "msksub")
            r.input_ok || push!(fails, "input")
            println(io, "| $(r.globalrow) | $(tele_name(r.telecode)) | $(r.mjd) | $(r.expnum) | $(r.adjfiberindx) | $(r.startind) | $(r.within) | $(join(fails, "+")) | $(r.shift_info) |")
        end
        length(failures) > 200 && println(io, "")
        length(failures) > 200 && println(io, "(… $(length(failures) - 200) more failures; see CSV)")
        println(io, "")
    end

    println(io, "## Notes")
    println(io, "")
    println(io, "- Exactness is bitwise (isequal, NaN-aware): batch `flux`/`fluxivar` must equal the ar1Dunical fiber row on every `simplemsk` pixel and be NaN elsewhere (producer stores `nanify(fspec[simplemsk], simplemsk)`).")
    println(io, "- Row identity comes exclusively from full_list_info.h5 + the batch filename; discovery order is never used.")
    println(io, "- Missing batches are reported, never skipped silently; a workup over this corpus MUST use the missing-rows mask (RowContract.missing_rows_mask) with sentinel fill.")
    println(io, "- Wall time: $(Dates.canonicalize(Dates.CompoundPeriod(now() - t0)))")
end

logmsg("Report written to $outpath")
logmsg("VERDICTS: completeness=$verdict_completeness rows=$verdict_rows shifts=$(length(shift_rows))")
isfile(ckptpath) && rm(ckptpath)   # report written; checkpoint no longer needed
rmprocs(workers())
