## Per-exposure sky-bundle cache
#
# WHY THIS EXISTS
# ---------------
# `getSky4visit` (src/ingest.jl) is called once per TARGET fiber. Its inner
# `sky_decomp` loop is genuinely per-fiber (it decomposes every surviving sky fiber
# against the TARGET fiber's per-fiber E5 priors) and therefore cannot be shared.
# Everything BEFORE that loop is per-EXPOSURE and is recomputed identically for every
# target fiber on the exposure:
#
#   1. `get_sky_fiber_indices`  - almanac lookup of the exposure's sky fibers
#   2. the ar1Dunical read      - `jldopen` loads the FULL 8700x300 flux/ivar/mask
#                                 arrays (44.4 MB) and then slices the sky columns
#   3. `select_sky_fibers`      - the M-SKY guard chain
#
# MEASURED on the 2026_05_01 corpus (see 2026_09_07/sky_cache/DESIGN.md): a median of
# 243 target fibers per exposure, i.e. that exposure-level work is repeated ~243x.
#
# The cache stores exactly the output of steps 1-3, with the bulk arrays reduced to the
# GUARD-SURVIVING sky columns only, so a hit replaces a 44.4 MB full-array read with a
# compact contiguous read.
#
# CONTRACT
# --------
# A hit must be indistinguishable from a miss. The cached quantities are a pure,
# deterministic function of (almanac entry, ar1Dunical data, skyZcut); nothing that
# depends on the PRIORS is cached, because nothing before the `sky_decomp` loop uses
# them. Repointing ARM_SKY_PRIOR_DIR (or any other prior) therefore cannot stale this
# cache -- see `sky_cache_key` for the full invalidation set.
#
# Disabled by default: with no cache root configured the runtime behaves exactly as it
# did before this file existed.

using SHA: sha256
using Random: randstring

"""
Cache format/semantics version. **BUMP THIS** whenever any of the following changes:
- the fields written by `write_sky_bundle!` / read by `read_sky_bundle`,
- `select_sky_fibers` / `validate_sky_fiber` (the guard verdict semantics),
- `get_sky_fiber_indices` (which almanac fibers count as sky),
- the meaning of any `SKY_*` / `SKYFIB_*` bit.
A bump makes every existing entry a miss (different key -> different path), so stale
entries can never be silently reused. `test/sky_cache.jl` pins this constant against a
checksum of the guard source so a silent semantic change fails CI.
"""
# v2: sky-fiber selection now hard-masks fibers AR flagged throughput-broken BEFORE
#     the z-cut runs (SKYFIB_RELTHRPT_BIT / SKY_RELTHRPT_FIBER_BIT). Every v1 entry
#     encodes the old, pre-filter-free verdict and must miss.
const SKY_CACHE_SCHEMA = 2

"""
    sky_cache_root()

Cache root directory, or `nothing` when caching is disabled (the default). Set
`ARM_SKY_CACHE_DIR` to enable. An empty value disables the cache.
"""
function sky_cache_root()
    v = get(ENV, "ARM_SKY_CACHE_DIR", "")
    isempty(v) ? nothing : v
end

# Once a write or mkpath fails we stop trying (a read-only or full cache directory must
# not turn into one warning per spectrum). Per-process, deliberately not thread-safe
# beyond a plain Bool store: a lost update only costs one extra failed attempt.
const SKY_CACHE_WRITE_DISABLED = Ref(false)
const SKY_CACHE_WARNED = Ref(false)

function sky_cache_warn(msg)
    if !SKY_CACHE_WARNED[]
        SKY_CACHE_WARNED[] = true
        println("sky cache: $msg (this warning is printed once per worker; the run continues by recomputing)")
        flush(stdout)
    end
end

"""
    file_identity(path)

Cheap content proxy for a large input file: `(size, mtime_ns)`. Hashing the 108 MB
ar1Dunical payload on every lookup would cost more than the work being cached, and the
pipeline never rewrites an ar1Dunical file in place without changing its mtime (AR.jl
writes to a new path and renames). Returns `nothing` if the file is missing, which
forces a miss.
"""
function file_identity(path)
    st = try
        stat(path)
    catch
        return nothing
    end
    ispath(path) || return nothing
    return (st.size, st.mtime)
end

"""
    sky_cache_key(tele, mjd, expnum, ar1Dfname, almanacFile; skyZcut)

Hex digest over the FULL invalidation set:

| component            | why it is in the key                                          |
|----------------------|---------------------------------------------------------------|
| `SKY_CACHE_SCHEMA`   | format + guard-semantics version (manual bump, CI-pinned)      |
| `tele`               | the SAME (mjd, expnum) exists at BOTH observatories --- 31,307 |
|                      | such pairs measured in 57618-61230; without this they collide  |
| `mjd`, `expnum`      | identify the exposure within a telescope                       |
| ar1Dunical identity  | (size, mtime) --- a re-reduction changes the sky fiber data    |
| almanac identity     | (size, mtime) --- a re-built almanac changes which fibers are  |
|                      | sky, and the fiber ordering                                    |
| `skyZcut`            | the only tunable that changes the cached guard verdict         |

NOT in the key, deliberately: the prior set (`ARM_SKY_PRIOR_DIR`, chebmsk, starCont,
starLines), `sky_obs_thresh` and `min_fibers`. None of them is read before the
`sky_decomp` loop, so none of them can change a cached value. `min_fibers` is applied
by the caller AFTER the cache lookup, so it stays correct even if it changes.

Returns `nothing` when either input file cannot be stat'ed (-> forced miss).
"""
function sky_cache_key(tele, mjd, expnum, ar1Dfname, almanacFile; skyZcut)
    ar_id = file_identity(ar1Dfname)
    alm_id = file_identity(almanacFile)
    (isnothing(ar_id) || isnothing(alm_id)) && return nothing
    payload = join(string.((
        "arM-skybundle", SKY_CACHE_SCHEMA, tele, mjd, expnum,
        basename(ar1Dfname), ar_id[1], ar_id[2],
        basename(almanacFile), alm_id[1], alm_id[2],
        skyZcut,
    )), "\0")
    return bytes2hex(sha256(payload))[1:16]
end

"""
    sky_cache_path(root, tele, mjd, expnum, key)

`<root>/skybundle/v<SCHEMA>/<tele>/<mjd>/skybundle_<tele>_<mjd>_<expnum>_<key>.h5`

`tele` appears in BOTH the directory and the file name, so the measured
(mjd, expnum)-at-both-observatories collisions cannot alias. Sharding by
`<tele>/<mjd>` bounds a directory to one night's exposures at one telescope
(<= ~120 files measured), giving ~5k directories and ~154k files for the whole
57618-61230 range --- well inside the enforced ceph inode quotas.
"""
function sky_cache_path(root, tele, mjd, expnum, key)
    joinpath(root, "skybundle", "v$(SKY_CACHE_SCHEMA)", string(tele), string(mjd),
        "skybundle_$(tele)_$(mjd)_$(lpad(expnum, 4, "0"))_$(key).h5")
end

"""
    read_sky_bundle(path, key)

Read a cached bundle, or `nothing` on ANY problem (missing, truncated, wrong schema,
wrong key, inconsistent shapes). Fail-open by construction: the caller recomputes.
"""
function read_sky_bundle(path, key)
    isfile(path) || return nothing
    try
        return h5open(path, "r") do f
            (read(f["schema"]) == SKY_CACHE_SCHEMA) || return nothing
            # the key is also in the file name, but a stray/relinked file must not be
            # trusted on its name alone
            (read(f["key"]) == key) || return nothing
            skyfibIndxs = read(f["skyfibIndxs"])::Vector{Int}
            skyFibBits = read(f["skyFibBits"])::Vector{Int}
            mskSky = convert.(Bool, read(f["mskSky"]))
            skyBit = read(f["skyBit"])::Int
            nSkyFibers = read(f["nSkyFibers"])::Int
            ncand = length(skyfibIndxs)
            (length(skyFibBits) == ncand && length(mskSky) == ncand) || return nothing
            (count(mskSky) == nSkyFibers) || return nothing
            if nSkyFibers == 0
                npix = length(logUniWaveAPOGEE)
                return (skyfibIndxs=skyfibIndxs, skyFibBits=skyFibBits, mskSky=mskSky,
                    skyBit=skyBit, nSkyFibers=0,
                    survspec=zeros(npix, 0), survivar=zeros(npix, 0), survmsk=falses(npix, 0))
            end
            survspec = read(f["survspec"])::Matrix{Float64}
            survivar = read(f["survivar"])::Matrix{Float64}
            survmsk = convert.(Bool, read(f["survmsk"]))
            (size(survspec, 2) == nSkyFibers && size(survivar) == size(survspec) &&
             size(survmsk) == size(survspec)) || return nothing
            (size(survspec, 1) == length(logUniWaveAPOGEE)) || return nothing
            return (skyfibIndxs=skyfibIndxs, skyFibBits=skyFibBits, mskSky=mskSky,
                skyBit=skyBit, nSkyFibers=nSkyFibers,
                survspec=survspec, survivar=survivar, survmsk=survmsk)
        end
    catch e
        sky_cache_warn("unreadable entry $path ($(sprint(showerror, e)))")
        return nothing
    end
end

"""
    write_sky_bundle!(path, key, bundle)

Publish a bundle ATOMICALLY: write to a unique temporary name IN THE SAME DIRECTORY
(so the rename stays within one filesystem and is therefore atomic), then `rename(2)`
it into place. A reader can only ever observe the complete file or no file at all.

Concurrency model: duplicate COMPUTATION is accepted; there is no claim/lock. The
cached value is a pure function of the key inputs, so every concurrent writer produces
identical content and "last rename wins" is not a correctness hazard. The duplicated
work is ~0.5 s; a claim protocol would need TTLs and a reaper for claims orphaned by
preempted or OOM-killed workers, and a stale claim is a far worse failure than a
duplicated half-second. See DESIGN.md.

Returns `true` on publish, `false` on any failure (never throws).
"""
function write_sky_bundle!(path, key, bundle)
    SKY_CACHE_WRITE_DISABLED[] && return false
    dir = dirname(path)
    tmp = ""
    try
        mkpath(dir)
        tmp = joinpath(dir, "." * basename(path) * ".tmp." * string(getpid()) * "." * randstring(12))
        h5open(tmp, "w") do f
            f["schema"] = SKY_CACHE_SCHEMA
            f["key"] = key
            f["skyfibIndxs"] = collect(Int, bundle.skyfibIndxs)
            f["skyFibBits"] = collect(Int, bundle.skyFibBits)
            f["mskSky"] = UInt8.(bundle.mskSky)
            f["skyBit"] = Int(bundle.skyBit)
            f["nSkyFibers"] = Int(bundle.nSkyFibers)
            if bundle.nSkyFibers > 0
                f["survspec"] = bundle.survspec
                f["survivar"] = bundle.survivar
                f["survmsk"] = UInt8.(bundle.survmsk)
            end
        end
        # raw rename(2), NOT Base.mv(...; force=true): mv unlinks the destination first,
        # which would leave a window in which the entry does not exist. rename(2)
        # replaces the destination in one atomic step.
        Base.Filesystem.rename(tmp, path)
        return true
    catch e
        if !isempty(tmp)
            try
                isfile(tmp) && rm(tmp; force = true)
            catch
            end
        end
        SKY_CACHE_WRITE_DISABLED[] = true
        sky_cache_warn("cannot write $path ($(sprint(showerror, e))); caching disabled for this worker")
        return false
    end
end

"""
    compute_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile; skyZcut=10)

The per-EXPOSURE half of `getSky4visit`, with no prior dependence: almanac sky-fiber
lookup, the ar1Dunical read, and the M-SKY guard chain, reduced to the surviving sky
columns. This is the only thing the cache ever stores.
"""
function compute_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile; skyZcut = 10)
    skyfibIndxs = get_sky_fiber_indices(almanacFile, tele, mjd, expnum)
    npix = length(logUniWaveAPOGEE)
    if length(skyfibIndxs) == 0
        return (skyfibIndxs = Int[], skyFibBits = Int[], mskSky = Bool[],
            skyBit = SKY_NO_FIBERS_BIT | SKY_TOO_FEW_FIBERS_BIT, nSkyFibers = 0,
            survspec = zeros(npix, 0), survivar = zeros(npix, 0), survmsk = falses(npix, 0))
    end

    ar1Dfname = get_1Duni_name(reduxBase, tele, mjd, expnum)
    f = jldopen(ar1Dfname)
    skyspec = f["flux_1d"][:, skyfibIndxs]
    skyivar = f["ivar_1d"][:, skyfibIndxs]
    skymskmat = f["mask_1d"][:, skyfibIndxs]
    # AR's per-fiber throughput verdict for the SAME candidate sky columns, in the same
    # order. `nothing` when the reduction predates the flag, which DISABLES the
    # pre-filter rather than pretending every fiber is healthy.
    sky_bitmsk_relthrpt = if haskey(f, "relthrpt") && haskey(f, "bitmsk_relthrpt")
        bits = f["bitmsk_relthrpt"]   # (N_CHIPS, N_FIBERS)
        thrpt = f["relthrpt"]
        [begin
             acc = reduce(|, Int.(bits[:, ix]); init = 0)
             any(.!isfinite.(thrpt[:, ix])) && (acc |= AR_RELTHRPT_NOTFINITE_BIT)
             acc
         end
         for ix in skyfibIndxs]
    else
        nothing
    end
    close(f)

    mskSky, nSkyFibers, skyBit,
    skyFibBits = select_sky_fibers(skyspec, skyivar, skymskmat;
        skyZcut = skyZcut, bitmsk_relthrpt = sky_bitmsk_relthrpt)
    if any(b -> (b & SKYFIB_RELTHRPT_BIT) != 0, skyFibBits)
        skyBit |= SKY_RELTHRPT_FIBER_BIT
    end

    # NOTHING IS PRINTED HERE, ON PURPOSE (see `getSky4visit`). `skyBit` and
    # `skyFibBits` ARE the verdict, they are returned to the caller, and `skyBit` is
    # written out as a per-spectrum column of every batch product. The log is no
    # longer the record; `test/regression/arm_census.jl` in ApogeeReduction reads the
    # products instead.
    surv = findall(mskSky)
    return (skyfibIndxs = collect(Int, skyfibIndxs), skyFibBits = collect(Int, skyFibBits),
        mskSky = collect(Bool, mskSky), skyBit = Int(skyBit), nSkyFibers = Int(nSkyFibers),
        survspec = skyspec[:, surv], survivar = skyivar[:, surv],
        survmsk = Bool.(skymskmat[:, surv]))
end

"""
    get_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile; skyZcut=10, cache_root=sky_cache_root())

Cached `compute_sky_bundle`. Fail-open at every step: a disabled/unset cache root, a
missing or corrupt entry, an unwritable directory, or an unstattable input all fall
back to recomputation with (at most) one warning per worker.
"""
function get_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile;
        skyZcut = 10, cache_root = sky_cache_root())
    isnothing(cache_root) && return compute_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile; skyZcut = skyZcut)

    ar1Dfname = get_1Duni_name(reduxBase, tele, mjd, expnum)
    key = sky_cache_key(tele, mjd, expnum, ar1Dfname, almanacFile; skyZcut = skyZcut)
    if isnothing(key)
        # a missing ar1Dunical/almanac is a real error; let the normal path raise it
        return compute_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile; skyZcut = skyZcut)
    end
    path = sky_cache_path(cache_root, tele, mjd, expnum, key)

    hit = read_sky_bundle(path, key)
    isnothing(hit) || return hit

    bundle = compute_sky_bundle(reduxBase, tele, mjd, expnum, almanacFile; skyZcut = skyZcut)
    write_sky_bundle!(path, key, bundle)
    return bundle
end
