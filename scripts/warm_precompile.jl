## Serial precompile-cache warmer. Run ONCE on the head node BEFORE any job spawns workers.
#
# WHY THIS EXISTS (MEASURED, Slurm job 6995969, 2026-09-07):
#   The E5 sky rebuild ran 224 Julia workers on 7 nodes. Every worker started cold, and all
#   224 hit the SAME precompile pidfiles on the shared GPFS depot at the same instant:
#       ApogeeReduction Being precompiled by another process (pid: 1432087, pidfile:
#       /mnt/home/asaydjari/.julia/compiled/v1.11/ApogeeReduction/....ji.pidfile)
#   219 of the 224 workers ended up precompiling ApogeeReduction independently (plus 36x
#   CairoMakie, 26x Optim, 33x two DifferentiationInterface/ArrayInterface extensions, ...).
#   Julia reported 92 s for the first and 154-171 s for the later ones -- the lock contention
#   makes the job SLOWER the more workers you add. Measured cost of that phase:
#       "Worker loading took 15 minutes, 3 seconds, 532 milliseconds"
#   Loading is not the expensive part; *concurrent cold precompilation* is. Doing it once,
#   serially, in a single process first collapses the whole phase.
#
# WHY A SCAN AND NOT A HARDCODED LIST:
#   A copied, hand-maintained package list is the same failure mode as a copied launcher
#   body -- it silently drifts from the thing it is supposed to mirror (that is exactly how
#   the SLURM_NTASKS 2x bug arose). So this script DERIVES the set: it scans the files it is
#   given for top-level `using`/`import` lines and keeps the names that are real
#   dependencies of the active project (or loadable stdlibs). Add a package to a driver and
#   the warmer picks it up with no edit here.
#
# WHY NOT JUST `Pkg.precompile()`:
#   It works, and it cannot miss anything, but it precompiles the ENTIRE manifest -- Makie,
#   Korg, Distributions extensions and everything else no sky-prior worker ever loads.
#   MEASURED on this project, 2026-09-07, ccalin051, julia 1.11.0:
#       Pkg.precompile()                 219 s cold -> 30.2 s warm  (298 packages)
#       using ApogeeReduction (targeted) 274 s cold ->  8.7 s warm  (13 packages)
#   Targeted is the default. Set AR_WARM_FULL=1 to ALSO run Pkg.precompile() when you want
#   belt-and-braces completeness at the extra ~30 s.
#
# Loading (not merely precompiling) the set is deliberate: it proves on the HEAD node, in
# under a minute, that every package the workers are about to need actually loads. A broken
# dependency then fails here instead of on worker 137 after the allocation is already spent.
#
# This script NEVER fails a job. Any error is reported and swallowed; it always exits 0.
# Warming is an optimisation, and an optimisation must not be able to take down a run.
#
# Usage:
#   julia --project=<proj> scripts/warm_precompile.jl <file-or-dir> [<file-or-dir> ...]
# Env:
#   AR_WARM_PRECOMPILE=0   skip entirely
#   AR_WARM_FULL=1         also run Pkg.precompile() over the whole manifest
# Author - Andrew Saydjari

using Dates

const T0 = time()

function main()
    if get(ENV, "AR_WARM_PRECOMPILE", "1") == "0"
        println("warm: AR_WARM_PRECOMPILE=0, skipping")
        return
    end

    proj = Base.active_project()
    println("warm: project   = ", proj)
    println("warm: julia     = ", VERSION)
    # If the depot is node-local rather than shared, warming on the head buys nothing for
    # the other nodes. Print it so that assumption is visible in the job log, not implicit.
    println("warm: depot[1]  = ", first(DEPOT_PATH))

    # ---- the project's declared dependencies (the filter for the scan) ----
    deps = Set{String}()
    try
        indeps = false
        for line in eachline(proj)
            s = strip(line)
            if startswith(s, "[")
                indeps = (s == "[deps]")
                continue
            end
            if indeps && occursin("=", s)
                push!(deps, strip(split(s, "=")[1]))
            end
        end
    catch e
        println("warm: WARNING could not parse [deps] from $proj ($e)")
    end
    println("warm: project declares ", length(deps), " direct dependencies")

    # ---- collect the files to scan, FOLLOWING include() chains ----
    # Following includes is what keeps this honest. A driver's real package set is not in
    # the driver: e5_sky_run.jl's `@everywhere` block never names ApogeeReduction -- it
    # arrives via `include("src/ingest.jl")`, whose line 5 is
    # `import ApogeeReduction: get_fibTargDict, ...`. That single package is what 219 of the
    # 224 workers in job 6995969 precompiled independently. A launcher therefore passes only
    # the driver, and the include graph supplies the rest; there is no list to keep in sync.
    projdir = dirname(proj)
    files = String[]
    seen = Set{String}()
    queue = String[]
    for a in ARGS
        if isdir(a)
            for (root, _, fs) in walkdir(a), f in fs
                endswith(f, ".jl") && push!(queue, joinpath(root, f))
            end
        elseif isfile(a)
            push!(queue, a)
        else
            println("warm: WARNING scan target does not exist, ignoring: ", a)
        end
    end
    # `include(src_code * "src/utils.jl")` / `include(joinpath(src_dir, "src/priors.jl"))`:
    # take every quoted fragment on an include line and try it against the project root and
    # against the including file's own directory. Unresolvable fragments are simply skipped.
    rx_inc = r"^\s*include\s*\("
    while !isempty(queue)
        f = popfirst!(queue)
        rf = try
            realpath(f)
        catch
            f
        end
        (rf in seen) && continue
        push!(seen, rf)
        push!(files, rf)
        try
            for line in eachline(rf)
                occursin(rx_inc, line) || continue
                for m in eachmatch(r"\"([^\"]+\.jl)\"", line)
                    frag = m.captures[1]
                    for cand in (joinpath(projdir, frag), joinpath(dirname(rf), frag), frag)
                        if isfile(cand)
                            push!(queue, cand)
                            break
                        end
                    end
                end
            end
        catch e
            println("warm: WARNING could not follow includes in $rf ($e)")
        end
    end
    isempty(files) && (println("warm: no files to scan, nothing to do"); return)
    println("warm: scanning ", length(files), " file(s) (include graph followed) for top-level using/import")

    # ---- scan for `using A, B` / `import A: x, y` at the start of a line ----
    # Deliberately conservative: only lines whose first token is using/import, and only the
    # part before any `:` (so `import ApogeeReduction: f, g` contributes ApogeeReduction and
    # not f/g). Anything that survives must also be a declared dep or a loadable stdlib, so a
    # false positive cannot do harm.
    found = Set{String}()
    rx = r"^\s*(using|import)\s+([^#]*)"
    for f in files
        try
            for line in eachline(f)
                m = match(rx, line)
                m === nothing && continue
                body = split(m.captures[2], ":")[1]
                for tok in split(body, ",")
                    n = strip(tok)
                    (isempty(n) || occursin(" ", n) || occursin(".", n)) && continue
                    occursin(r"^[A-Za-z][A-Za-z0-9_]*$", n) || continue
                    push!(found, n)
                end
            end
        catch e
            println("warm: WARNING could not scan $f ($e)")
        end
    end

    # keep declared deps; also keep names that resolve as stdlibs in this environment
    keep = String[]
    for n in sort(collect(found))
        if n in deps
            push!(keep, n)
        elseif Base.identify_package(n) !== nothing
            push!(keep, n)   # stdlib or otherwise resolvable in the active env
        end
    end
    println("warm: warming ", length(keep), " package(s): ", join(keep, " "))
    flush(stdout)

    if get(ENV, "AR_WARM_FULL", "0") == "1"
        t = time()
        try
            @eval Main using Pkg
            Base.invokelatest(Main.Pkg.precompile)
            println("warm: Pkg.precompile() full manifest OK in ", round(time() - t, digits = 1), " s")
        catch e
            println("warm: WARNING Pkg.precompile() failed after ", round(time() - t, digits = 1), " s: ", e)
        end
        flush(stdout)
    end

    nfail = 0
    for n in keep
        t = time()
        try
            Core.eval(Main, Meta.parse("using $n"))
            dt = time() - t
            # only report the ones that cost something; a warm load is ~0.0-0.5 s
            dt >= 1.0 && println(rpad("warm:   $n", 40), round(dt, digits = 1), " s")
        catch e
            nfail += 1
            println("warm:   $n FAILED after ", round(time() - t, digits = 1), " s: ",
                sprint(showerror, e)[1:min(end, 300)])
        end
        flush(stdout)
    end
    println("warm: ", length(keep) - nfail, "/", length(keep), " packages loaded, ",
        round(time() - T0, digits = 1), " s total")
    if nfail > 0
        println("warm: WARNING ", nfail, " package(s) failed to load on the head node. The ",
            "workers are about to try the same thing. Investigate before trusting this run.")
    end
end

try
    main()
catch e
    println("warm: WARNING warmer itself failed (", sprint(showerror, e)[1:min(end, 500)], ")")
    println("warm: continuing -- warming is an optimisation and must never fail a job")
end
flush(stdout)
exit(0)
