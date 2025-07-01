### File Name Handlers

chipRaw2Redux = Dict(
    "a" => "R",
    "b" => "G",
    "c" => "B",
)

chipRedux2Raw = Dict(
    "R" => "a",
    "G" => "b",
    "B" => "c",
)

## AR.jl
function get_1Duni_name(reduxBase, tele, mjd, expnum; fnameType = "ar1Dunical")
    return joinpath(reduxBase, "apred", mjd, join([fnameType, tele, mjd, lpad(expnum,4,"0"), "OBJECT.h5"], "_"))
end

function adjfiberindx2fiberindx(adjfiberindx::Int)
    if adjfiberindx >300
        return adjfiberindx-300
    else
        return adjfiberindx
    end
end

function get_almanac_file(reduxBase, mjd)
    return joinpath(reduxBase, "almanac", join(["objects", mjd*".h5"], "_"))
end