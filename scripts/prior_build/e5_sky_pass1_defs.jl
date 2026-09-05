## E5 pass-1 sky-prior retraining definitions (DR21 testbed corpus).
#
# Implements the AKS-decided (2026-09-04) pooling design for ALL 600 fibers:
#
#   POOLING RULE (AKS-decided 2026-09-04): the sample list for target fiber f
#   (native fiberindx 1:300, per telescope) pools every validated sky sample from
#   fibers g with |g - f| <= E5_WIND (= 5, the DR17 `wind = 5` window from
#   apMADGICS src/verify_drp.jl at Utah), HARD-CLIPPED to f's own MTP block —
#   contiguous 30-fiber v-groove blocks, fiberindx 1-30, 31-60, ..., 271-300
#   (mtp_block(x) = fld(x-1,30)+1; same mapping the production wavelength solution
#   uses, ApogeeReduction.jl/src/wavecal.jl "mtp_inds"). Pools NEVER cross block
#   or telescope boundaries. No other exclusion applies (DR17 had none: the
#   apparent "skips" in the DR17 lists, e.g. fiberindx 214/220 absent from file
#   519, are pure availability — those fibers had zero plate-era sky assignments,
#   dr17_dr17_skycounts_raw.txt; DR17 itself did NOT clip at MTP boundaries,
#   e.g. its file 030 pools fiberindx 25-35 across the 30|31 boundary — the
#   clipping here is new, AKS-specified behavior).
#
# Pooling is implemented in the LIST-GENERATION step (e5_make_pooled_lists), as
# DR17 did at Utah (verify_drp.jl wrote per-fiber pooled sky_input lists that
# sample_sky.jl consumed verbatim) — not as a patch on the sampler.
#
# Execution avoids the naive ~10M redundant exposure re-reads by splitting the
# merged getSkyRough(fibindxoi=...) path into two provably-equivalent stages:
#   1. e5_extract_packs: per unique (tele,mjd,expnum), run the EXACT merged guard
#      chain (combine_sky_fibers: validate_sky_fiber bits + skyZcut) ONCE and
#      cache all candidate sky columns + guard verdicts in a "pack" HDF5.
#      combine_sky_fibers is deterministic and depends only on the exposure's full
#      candidate set (never on fibindxoi), so pack contents reproduce what
#      getSkyRough would hand every caller.
#   2. e5_gather_source_samples: per SOURCE fiber, collect its surviving columns
#      (mskSky[j] && fiberindx match — the same condition getSkyRough applies
#      before returning a column) into per-source files, in deterministic
#      (mjd, expnum) order.
#   3. e5_assemble_pooled: per TARGET fiber, hcat the pool members' per-source
#      samples (minus screened outliers) into the standard sampler output files
#      (skyflux_NNN.h5 / skyivar_NNN.h5 / skymsk_NNN.h5 + meta_* datasets), in
#      deterministic (mjd, expnum, fiberindx) order. Downstream, the MERGED
#      decomposition (get_sky_samples' second stage, target fiber's starCont
#      basis on all pooled samples — DR17-faithful) and the MERGED builders
#      (build_skyCont / build_skyLines) run unmodified on these files.
#
# Outlier screens (AKS-agreed, E6 SVD-amplifies-rare-outliers lesson) are DESIGNED
# from per-source-fiber sample statistics (e5_source_stats) and APPLIED at
# assembly time (a screened sample is dropped from every pool it appears in).
#
# Requires (from includer): HDF5, LinearAlgebra, StatsBase, Distributed helpers,
# src/utils.jl, src/fileNameHandling.jl, src/ingest.jl,
# scripts/prior_build/prior_utils.jl (nanzeromedian).
# Author - Andrew Saydjari (E5 pass 1)

const E5_WIND = 5

"MTP (slit v-groove) block id for a native fiberindx 1:300 (blocks 1-30,31-60,...)."
mtp_block(fiberindx::Int) = fld(fiberindx - 1, 30) + 1

"In-block +/-E5_WIND pool member range for a native fiberindx (AKS rule above)."
function e5_pool_range(fiberindx::Int; wind::Int=E5_WIND)
    blo = 30 * fld(fiberindx - 1, 30) + 1
    bhi = blo + 29
    return max(fiberindx - wind, blo):min(fiberindx + wind, bhi)
end

"""
    e5_collect_sky_runlist(almanacFile; runlist_parallel=true)

All exact-fiber sky entries in the corpus: vcat of the merged
get_telemjd_runlist_from_almanac(accepted_fibtypes=["sky"]) over every
(tele, mjd) in the almanac. The DR21 almanac stores telescopes under a `raw/`
root, which the installed ApogeeReduction readers (read_almanac_exp_df,
get_fibTargDict) expect natively — the file is passed straight through.
Entries are the merged NamedTuples (tele, mjd, expnum, adjfiberindx, sdss_id).
"""
function e5_collect_sky_runlist(almanacFile; runlist_parallel=true)
    tele_mjd_pairs = Tuple{String,String}[]
    h5open(almanacFile, "r") do f
        root = haskey(f, "raw") ? f["raw"] : f
        for tele in ("apo", "lco")
            haskey(root, tele) || continue
            for mjd in keys(root[tele])
                push!(tele_mjd_pairs, (tele, mjd))
            end
        end
    end
    get_runlist_partial(argtup) = get_telemjd_runlist_from_almanac(
        almanacFile, argtup[1], argtup[2], accepted_fibtypes=["sky"])
    run_lsts = if runlist_parallel
        pmap(get_runlist_partial, tele_mjd_pairs)
    else
        map(get_runlist_partial, tele_mjd_pairs)
    end
    return vcat(run_lsts...)
end

# deterministic sort key for sample ordering
e5_sortkey(e) = (e.tele, parse(Int, e.mjd), e.expnum, e.adjfiberindx)

"""
    e5_make_pooled_lists(almanacFile, out_dir; wind=E5_WIND, runlist_parallel=true)

List-generation step (the pooling lives HERE, per AKS). Writes:
  out_dir/e5_sky_pool_lst_NNN.h5  (NNN = adjfiberindx 001-600) with datasets
    tele (String), mjd (String), expnum (Int), fiberindx (Int, native 1:300
    SOURCE fiber), sdss_id (Int) — sorted by (mjd, expnum, fiberindx) —
    and attributes documenting the rule;
  out_dir/e5_pool_summary.h5 with per-target pool membership and counts.
Returns the per-target counts vector (600).
"""
function e5_make_pooled_lists(almanacFile, out_dir; wind=E5_WIND, runlist_parallel=true)
    mkpath(out_dir)
    run_lst = e5_collect_sky_runlist(almanacFile; runlist_parallel=runlist_parallel)
    sort!(run_lst, by=e5_sortkey)

    # index entries by (teleind, native fiberindx)
    by_fiber = [NamedTuple[] for _ in 1:600]
    for e in run_lst
        push!(by_fiber[e.adjfiberindx], e)
    end

    counts = zeros(Int, 600)
    members_lo = zeros(Int, 600); members_hi = zeros(Int, 600)
    for adjfib in 1:600
        offset = adjfib > 300 ? 300 : 0
        f = adjfib - offset
        pr = e5_pool_range(f; wind=wind)
        members_lo[adjfib] = first(pr); members_hi[adjfib] = last(pr)
        pooled = reduce(vcat, [by_fiber[g + offset] for g in pr]; init=NamedTuple[])
        sort!(pooled, by=e5_sortkey)
        counts[adjfib] = length(pooled)
        fname = joinpath(out_dir, "e5_sky_pool_lst_" * lpad(adjfib, 3, "0") * ".h5")
        h5open(fname, "w") do fh
            fh["tele"] = String[e.tele for e in pooled]
            fh["mjd"] = String[e.mjd for e in pooled]
            fh["expnum"] = Int[e.expnum for e in pooled]
            fh["fiberindx"] = Int[e.adjfiberindx - offset for e in pooled]
            fh["sdss_id"] = Int[e.sdss_id for e in pooled]
            attrs(fh)["pooling_rule"] = "AKS 2026-09-04: |g-f|<=$(wind) AND mtp_block(g)==mtp_block(f) (blocks fld(x-1,30)+1, 1-30/31-60/.../271-300), same telescope; no other exclusions. DR17 wind=5 window (verify_drp.jl) with NEW hard MTP-block clipping."
            attrs(fh)["target_adjfiberindx"] = adjfib
            attrs(fh)["pool_members_native"] = collect(pr)
        end
    end
    h5open(joinpath(out_dir, "e5_pool_summary.h5"), "w") do fh
        fh["adjfiberindx"] = collect(1:600)
        fh["pool_lo"] = members_lo
        fh["pool_hi"] = members_hi
        fh["n_pooled_candidates"] = counts
    end
    return counts
end

"""
    e5_pack_path(pack_dir, tele, mjd, expnum)

Per-exposure guard-cache pack file path.
"""
e5_pack_path(pack_dir, tele, mjd, expnum) =
    joinpath(pack_dir, tele, "skypack_" * tele * "_" * mjd * "_" * lpad(expnum, 4, "0") * ".h5")

"""
    e5_extract_pack(reduxBase, almanacFile, tele, mjd, expnum; pack_dir)

Run the merged guard chain once for one exposure and cache the result.
Mirrors getSkyRough(...; fibindxoi) exactly: same fibtargDict -> skyfibIndxs
selection, same column reads, same combine_sky_fibers(skyZcut=10) verdicts.
Stores ALL candidate columns plus mskSky/skyFibBits so any (target, source)
lookup reproduces the merged return. Checkpointed on file existence.
"""
function e5_extract_pack(reduxBase, almanacFile, tele, mjd, expnum; pack_dir)
    pname = e5_pack_path(pack_dir, tele, mjd, expnum)
    isfile(pname) && return :exists
    f = h5open(almanacFile)
    fibtargDict, _ = get_fibTargDict(f, tele, mjd, expnum)
    close(f)
    fibtypelist = map(x -> fibtargDict[x], 1:300)
    skyfibIndxs = findall(map(x -> x[1:3] == "sky", fibtypelist))
    isempty(skyfibIndxs) && return :nosky # cannot happen for exposures on the pooled lists

    ar1Dfname = get_1Duni_name(reduxBase, tele, mjd, expnum)
    fj = jldopen(ar1Dfname)
    skyspec = fj["flux_1d"][:, skyfibIndxs]
    skyivar = fj["ivar_1d"][:, skyfibIndxs]
    skymsk = fj["mask_1d"][:, skyfibIndxs]
    close(fj)

    nSkyFibers, _, _, _, skyBit, mskSky, skyFibBits =
        combine_sky_fibers(skyspec, skyivar, skymsk; skyZcut=10)

    # per-candidate scale (for screen design; nanzeromedian as in the guard)
    scales = [nanzeromedian(view(skyspec, :, j)) for j in eachindex(skyfibIndxs)]

    mkpath(dirname(pname))
    tmpname = pname * ".tmp"
    h5open(tmpname, "w") do fh
        fh["flux"] = skyspec
        fh["ivar"] = skyivar
        # strict Bool conversion (as in the merged sampler's zeros(Bool) assignment),
        # materialized as Array: JLD2/broadcast hand back BitArray, which HDF5.jl
        # cannot write (MethodError: strides)
        fh["msk"] = Matrix{Bool}(skymsk)
        fh["fiberindx"] = skyfibIndxs
        fh["mskSky"] = Vector{Bool}(mskSky)
        fh["skyFibBits"] = Int.(skyFibBits)
        fh["scale"] = convert.(Float64, scales)
        attrs(fh)["skyBit"] = skyBit
        attrs(fh)["nSkyFibers"] = nSkyFibers
        attrs(fh)["tele"] = tele; attrs(fh)["mjd"] = mjd; attrs(fh)["expnum"] = expnum
    end
    mv(tmpname, pname; force=true)
    return :done
end

"""
    e5_unique_exposures(list_dir)

Unique (tele, mjd, expnum) across all 600 pooled lists (= across the exact-fiber
corpus; pooling adds no new exposures).
"""
function e5_unique_exposures(list_dir)
    seen = Set{Tuple{String,String,Int}}()
    for adjfib in 1:600
        fname = joinpath(list_dir, "e5_sky_pool_lst_" * lpad(adjfib, 3, "0") * ".h5")
        h5open(fname, "r") do fh
            tele = read(fh["tele"]); mjd = read(fh["mjd"]); expnum = read(fh["expnum"])
            for i in eachindex(tele)
                push!(seen, (tele[i], mjd[i], expnum[i]))
            end
        end
    end
    out = collect(seen)
    sort!(out, by=x -> (x[1], parse(Int, x[2]), x[3]))
    return out
end

"""
    e5_gather_source_samples(adjfib, list_dir, pack_dir, src_dir)

Per-SOURCE-fiber sample gather: for native source fiber f = adjfib (mod 300),
read every pack where f was a sky candidate (from fiber f's own pooled list
restricted to fiberindx == f) and keep the columns where the guard passed
(mskSky). Writes src_dir/skysrc_NNN.h5 with flux/ivar/msk + meta (mjd, expnum,
scale) in (mjd, expnum) order. Checkpointed. Returns (n_candidates, n_kept).
"""
function e5_gather_source_samples(adjfib, list_dir, pack_dir, src_dir)
    offset = adjfib > 300 ? 300 : 0
    fsrc = adjfib - offset
    sname = joinpath(src_dir, "skysrc_" * lpad(adjfib, 3, "0") * ".h5")
    lname = joinpath(list_dir, "e5_sky_pool_lst_" * lpad(adjfib, 3, "0") * ".h5")
    tele_l, mjd_l, exp_l, fib_l = h5open(lname, "r") do fh
        read(fh["tele"]), read(fh["mjd"]), read(fh["expnum"]), read(fh["fiberindx"])
    end
    own = findall(fib_l .== fsrc)
    if isfile(sname)
        nkept = h5open(sname, "r") do fh; size(fh["flux"], 2); end
        return (length(own), nkept)
    end
    npix = length(wavetarg)
    flux = zeros(npix, length(own)); ivar = zeros(npix, length(own))
    msk = zeros(Bool, npix, length(own))
    mjds = zeros(Int, length(own)); exps = zeros(Int, length(own))
    scales = zeros(length(own))
    keep = falses(length(own))
    for (k, i) in enumerate(own)
        pname = e5_pack_path(pack_dir, tele_l[i], mjd_l[i], exp_l[i])
        isfile(pname) || error("e5_gather_source_samples: missing pack $pname")
        h5open(pname, "r") do fh
            fibs = read(fh["fiberindx"])
            j = findfirst(fibs .== fsrc)
            isnothing(j) && error("e5_gather_source_samples: fiber $fsrc not a candidate in $pname (list/pack mismatch)")
            if read(fh["mskSky"])[j]
                flux[:, k] .= fh["flux"][:, j]
                ivar[:, k] .= fh["ivar"][:, j]
                msk[:, k] .= fh["msk"][:, j]
                mjds[k] = parse(Int, mjd_l[i]); exps[k] = exp_l[i]
                scales[k] = read(fh["scale"])[j]
                keep[k] = true
            end
        end
    end
    kidx = findall(keep)
    tmpname = sname * ".tmp"
    h5open(tmpname, "w") do fh
        fh["flux"] = flux[:, kidx]
        fh["ivar"] = ivar[:, kidx]
        fh["msk"] = msk[:, kidx]
        fh["mjd"] = mjds[kidx]
        fh["expnum"] = exps[kidx]
        fh["scale"] = scales[kidx]
        attrs(fh)["adjfiberindx"] = adjfib
        attrs(fh)["n_candidates"] = length(own)
        attrs(fh)["n_kept_guard"] = length(kidx)
    end
    mv(tmpname, sname; force=true)
    return (length(own), length(kidx))
end

"""
    e5_source_stats(adjfib, src_dir)

Per-sample screening statistics for one source fiber, from its skysrc file:
  scale       nanzeromedian flux (airglow brightness proxy; bright-outlier screen)
  goodfrac    fraction of pixels with msk & ivar>0
  chi2r_full  mean over good pixels of ivar*(flux - scale*ref)^2, ref = per-pixel
              median of flux/scale over the fiber's samples (chi2-class screen)
  chi2r_cont  same restricted to continuum pixels (ref below its 60th percentile;
              targets continuum-shape anomalies like the E4 light-leak case
              without penalizing natural airglow-line variability)
Returns (mjd, expnum, scale, goodfrac, chi2r_full, chi2r_cont, nsamp).
"""
function e5_source_stats(adjfib, src_dir)
    sname = joinpath(src_dir, "skysrc_" * lpad(adjfib, 3, "0") * ".h5")
    flux, ivar, msk, mjd, expnum, scale = h5open(sname, "r") do fh
        read(fh["flux"]), read(fh["ivar"]), read(fh["msk"]), read(fh["mjd"]),
        read(fh["expnum"]), read(fh["scale"])
    end
    npix, n = size(flux)
    goodfrac = zeros(n); chi2r_full = fill(NaN, n); chi2r_cont = fill(NaN, n)
    if n >= 10
        # reference shape: per-pixel median of scale-normalized flux over samples
        good = msk .& (ivar .> 0)
        normed = flux ./ reshape(max.(scale, 1e-3), 1, :)
        normed[.!good] .= NaN
        ref = [nanzeromedian(view(normed, i, :)) for i in 1:npix]
        refv = replace(ref, NaN => 0.0)
        refq = quantile(filter(x -> isfinite(x) && x != 0, refv), 0.6)
        contpix = isfinite.(ref) .& (refv .< refq) .& (refv .!= 0)
        for k in 1:n
            g = view(good, :, k)
            goodfrac[k] = count(g) / npix
            if count(g) > 100
                r2 = (view(flux, :, k) .- scale[k] .* refv) .^ 2 .* view(ivar, :, k)
                chi2r_full[k] = mean(r2[g])
                gc = g .& contpix
                chi2r_cont[k] = count(gc) > 100 ? mean(r2[gc]) : NaN
            end
        end
    end
    return (mjd=mjd, expnum=expnum, scale=scale, goodfrac=goodfrac,
        chi2r_full=chi2r_full, chi2r_cont=chi2r_cont, nsamp=n)
end

"""
    e5_apply_screens(stats; scale_zmax, chi2cont_max, chi2full_max)

Screen verdict per sample from e5_source_stats output for ONE source fiber.
bright-outlier: robust z of log10(scale) within the fiber's own distribution
(z = (x - median)/(1.4826*MAD)) above scale_zmax OR non-positive scale with
|scale| large is left to the guard (negscale is recorded-only per M-SKY).
chi2-class: chi2r_cont > chi2cont_max or chi2r_full > chi2full_max.
Returns (drop::BitVector, reason::Vector{Int}) with bit 1 = bright, bit 2 = chi2.
"""
function e5_apply_screens(stats; scale_zmax=6.0, chi2cont_max=25.0, chi2full_max=200.0)
    n = stats.nsamp
    drop = falses(n); reason = zeros(Int, n)
    pos = stats.scale .> 0
    if count(pos) >= 10
        lx = log10.(stats.scale[pos])
        med = median(lx); madv = 1.4826 * median(abs.(lx .- med))
        madv = max(madv, 1e-3)
        z = fill(0.0, n)
        z[pos] .= (log10.(stats.scale[pos]) .- med) ./ madv
        for k in 1:n
            if pos[k] && (z[k] > scale_zmax)
                drop[k] = true; reason[k] |= 1
            end
        end
    end
    for k in 1:n
        if (isfinite(stats.chi2r_cont[k]) && stats.chi2r_cont[k] > chi2cont_max) ||
           (isfinite(stats.chi2r_full[k]) && stats.chi2r_full[k] > chi2full_max)
            drop[k] = true; reason[k] |= 2
        end
    end
    return drop, reason
end

"""
    e5_assemble_pooled(adjfib, list_dir, src_dir, sample_dir; drop_dict=Dict(), variant_unscreened=false)

Assemble the per-TARGET pooled sample files consumed by the merged decomposition
and builders: sample_dir/{skyflux,skyivar,skymsk}_NNN.h5 (standard datasets) plus
meta_{fiberindx,mjd,expnum} datasets in the skyflux file. Pool membership from
the target's own pooled list; per-source samples from skysrc files; screened
samples (drop_dict[srcadjfib] :: BitVector over that source's kept samples)
dropped unless variant_unscreened. Deterministic (mjd, expnum, fiberindx) order.
Checkpointed on skyflux existence. Returns number of assembled samples.
"""
function e5_assemble_pooled(adjfib, list_dir, src_dir, sample_dir;
        drop_dict=Dict{Int,BitVector}(), variant_unscreened=false)
    savefluxname = joinpath(sample_dir, "skyflux_" * lpad(adjfib, 3, "0") * ".h5")
    saveivarname = joinpath(sample_dir, "skyivar_" * lpad(adjfib, 3, "0") * ".h5")
    savemskname = joinpath(sample_dir, "skymsk_" * lpad(adjfib, 3, "0") * ".h5")
    if isfile(savefluxname) && isfile(saveivarname) && isfile(savemskname)
        return h5open(savefluxname, "r") do fh; size(fh["skyflux"], 2); end
    end
    mkpath(sample_dir)
    offset = adjfib > 300 ? 300 : 0
    lname = joinpath(list_dir, "e5_sky_pool_lst_" * lpad(adjfib, 3, "0") * ".h5")
    pool = h5open(lname, "r") do fh; attrs(fh)["pool_members_native"][:]; end

    blocks = [] # (fiberindx, mjd, expnum, flux, ivar, msk, keepmask)
    for g in pool
        srcadj = g + offset
        sname = joinpath(src_dir, "skysrc_" * lpad(srcadj, 3, "0") * ".h5")
        isfile(sname) || error("e5_assemble_pooled: missing source file $sname")
        h5open(sname, "r") do fh
            n = size(fh["flux"], 2)
            keep = trues(n)
            if !variant_unscreened && haskey(drop_dict, srcadj)
                dd = drop_dict[srcadj]
                length(dd) == n || error("e5_assemble_pooled: drop mask length $(length(dd)) != nsamp $n for src $srcadj")
                keep .= .!dd
            end
            kidx = findall(keep)
            if !isempty(kidx)
                # read full then column-index (HDF5 fancy indexing avoided)
                fl = read(fh["flux"]); iv = read(fh["ivar"]); mk = read(fh["msk"])
                push!(blocks, (g, read(fh["mjd"])[kidx], read(fh["expnum"])[kidx],
                    fl[:, kidx], iv[:, kidx], mk[:, kidx]))
            end
        end
    end
    isempty(blocks) && error("e5_assemble_pooled: 0 samples for target adjfib=$adjfib")
    fibv = reduce(vcat, [fill(b[1], length(b[2])) for b in blocks])
    mjdv = reduce(vcat, [b[2] for b in blocks])
    expv = reduce(vcat, [b[3] for b in blocks])
    flux = reduce(hcat, [b[4] for b in blocks])
    ivar = reduce(hcat, [b[5] for b in blocks])
    msk = reduce(hcat, [b[6] for b in blocks])
    ord = sortperm(collect(zip(mjdv, expv, fibv)))
    h5open(savefluxname * ".tmp", "w") do fh
        fh["skyflux"] = flux[:, ord]
        fh["meta_fiberindx"] = fibv[ord]
        fh["meta_mjd"] = mjdv[ord]
        fh["meta_expnum"] = expv[ord]
        attrs(fh)["screened"] = !variant_unscreened
    end
    h5open(saveivarname * ".tmp", "w") do fh; fh["skyivar"] = ivar[:, ord]; end
    h5open(savemskname * ".tmp", "w") do fh; fh["skymsk"] = msk[:, ord]; end
    for nm in (savefluxname, saveivarname, savemskname)
        mv(nm * ".tmp", nm; force=true)
    end
    return length(ord)
end
