# Tests for the AllExposure builder.   julia --project=scripts/release_tables/allexposure test/runtests.jl
#
# Unit tests run anywhere (no cluster data). The full-table regression runs only when
#   ALLEXPOSURE_REFERENCE=<released h5>  ALLEXPOSURE_CANDIDATE=<built h5>
# are set; it checks every dataset bitwise (the authoritative check is
# ../../validation/compare_h5.py, which also checks attribute string flavours and layouts).
using Test, HDF5, AllExposure
const AE = AllExposure

hexf(s) = parse(Float64, s)
function fixture()
    lines = filter(l -> !startswith(l, "#"), readlines(joinpath(@__DIR__, "fixture_rows_040.csv")))
    hdr = split(lines[1], ",")
    rows = [split(l, ",") for l in lines[2:end]]
    Dict(String(h) => (h == "row" ? [parse(Int, r[j]) for r in rows] : [hexf(r[j]) for r in rows])
         for (j, h) in enumerate(hdr))
end

@testset "AllExposure" begin
    @testset "schema" begin
        names = [c.name for c in AE.SCHEMA]
        @test length(unique(names)) == length(names) == 52
        @test all(c.layout in (:gz, :gz_nosh, :gz13045, :contig) for c in AE.SCHEMA)
        @test AE.chunk_rows(Float64) == 13045 && AE.chunk_rows(Int16) == 26090 && AE.chunk_rows(UInt8) == 52179
    end

    @testset "numpy casting semantics" begin
        @test AE._conv(Int32, true) === Int32(1)
        @test AE._conv(Int32, Int64(2)^32 + 5) === Int32(5)          # wraps like numpy
        @test AE._conv(Float32, 0.1) === Float32(0.1)
        @test AE._strs(["ab\0\0", "c"]) == ["ab", "c"]
    end

    @testset "categorical first-seen order" begin
        ct = AE.CodeTable()
        @test AE.code!(ct, ["science", "sky_apogee", "science"]) == UInt8[1, 2, 1]
        @test AE.code!(ct, ["standard_apogee", "science"]) == UInt8[3, 1]
        @test ct.order == ["", "science", "sky_apogee", "standard_apogee"]
    end

    @testset "fiber index convention" begin
        @test AE.fiber_index_almanac(Int8(1), Int16(1)) == 300
        @test AE.fiber_index_almanac(Int8(2), Int16(600)) == 1
    end

    fx = fixture()
    @testset "Galactic: bit-identical to astropy (0.4.0 fixture rows)" begin
        gl, gb = AE.galactic(fx["ra"], fx["dec"])
        @test all(isequal.(gl, fx["glon"]))
        @test all(isequal.(gb, fx["glat"]))
    end

    @testset "libm cos/sin used by LegacyNumerics" begin
        x = fx["ra"] .* AE.DEG2RAD_NP
        @test all(AE._cos(AE.LegacyNumerics(), v) == ccall(:cos, Float64, (Float64,), v) for v in x)
    end

    @testset "velocities: native numerics agree with 0.4.0 to 1e-9 km/s (not bitwise)" begin
        pw = AE.pow10(AE.NativeNumerics(), fx["RV_pixoff_final"], "")
        vrel, vrad, ev = AE.velocities(pw, fx["RV_pix_var"], fx["v_barycentric_correction"])
        ok = isfinite.(fx["v_rel"])
        @test maximum(abs.(vrel[ok] .- fx["v_rel"][ok])) < 1e-9
        @test maximum(abs.(vrad[ok] .- fx["v_rad"][ok])) < 1e-9
    end

    @testset "attribute string flavours" begin
        p = tempname() * ".h5"
        h5open(p, "w") do f
            AE.put_attr!(f, "fixed", AE.FixedStr("abc"))
            AE.put_attr!(f, "vlen", AE.VlenStr("abc"))
            AE.put_attr!(f, "vvec", AE.VlenStrVec(["", "x"]))
        end
        h5open(p, "r") do f
            a = attributes(f)
            @test read(a["fixed"]) == "abc" && read(a["vlen"]) == "abc" && read(a["vvec"]) == ["", "x"]
            dt(n) = HDF5.API.h5a_get_type(a[n])
            @test HDF5.API.h5t_is_variable_str(dt("fixed")) == false
            @test HDF5.API.h5t_is_variable_str(dt("vlen")) == true
        end
        rm(p)
    end

    ref = get(ENV, "ALLEXPOSURE_REFERENCE", ""); cand = get(ENV, "ALLEXPOSURE_CANDIDATE", "")
    if !isempty(ref) && !isempty(cand)
        @testset "full-table regression vs $(basename(ref))" begin
            h5open(ref, "r") do r; h5open(cand, "r") do c
                @test sort(collect(keys(r))) == sort(collect(keys(c)))
                for k in keys(r)
                    a = read(r[k]); b = read(c[k])
                    @test eltype(a) == eltype(b)
                    @test length(a) == length(b) && all(reinterpret(UInt8, a) .== reinterpret(UInt8, b))
                end
            end; end
        end
    end
end
