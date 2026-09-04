# E6 step 5b: rebuild adjfib 595's NEW prior with the leak-affected sample columns
# dropped (columns whose Tfunindx came from lco 59160/0018 or 60291/0009).
# Replicates build_starCont.jl's math exactly (mask, zero-sum col drop, grand-mean
# norm, C=VV'/nsamp, svd, nsub=60), on the SAME new sample file, minus 3 columns.
using HDF5, Serialization, LinearAlgebra, Statistics

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
const FIB = 595
const DROPCOLS = [668, 7106, 8208]  # from leak_identify.jl (RNG reproduction)
const NSUB = 60

samples = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/prior_outputs/starCont_20260903/tell_prior_disk/starCont_" *
          lpad(FIB, 3, "0") * ".jdat"
starcont = deserialize(samples)
println("loaded ", samples, " ", size(starcont))

chipgapmsk = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5", "r") do f
    Bool.(read(f[FIB > 300 ? "lco" : "apo"]))
end
wavelen = size(starcont, 1)

function build(starcont, chipgapmsk, savename)
    specsum = dropdims(sum(starcont, dims=1), dims=1)
    Vred = starcont[chipgapmsk, specsum .> 0]
    mnorm = mean(filter(!iszero, Vred))
    Vred ./= mnorm
    nsamp = size(Vred, 2)
    Cexp = Vred * transpose(Vred)
    Cexp ./= nsamp
    SF = svd(Cexp)
    EVEC = zeros(wavelen, size(SF.U, 2))
    EVEC[chipgapmsk, :] .= SF.U
    rm(savename; force=true)
    h5write(savename, "Vmat", EVEC[:, 1:NSUB] * Diagonal(sqrt.(SF.S[1:NSUB])))
    h5write(savename, "λv", SF.S[1:NSUB])
    h5write(savename, "chipgapmsk", collect(Bool, chipgapmsk))
    println("wrote ", savename, "  (nsamp used = ", nsamp, ")")
end

keepcols = setdiff(1:size(starcont, 2), DROPCOLS)
build(starcont[:, keepcols], chipgapmsk, joinpath(RUND, "built_new", "APOGEE_starcont_svd_60_f595_dropleak.h5"))
println("DONE")
