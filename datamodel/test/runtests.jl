# julia --project=. test/runtests.jl
using Test
using DataModel
using DataModel: normalize_dtype, shape_string, h5py_shape, sanitize, harvest_markdown_tables,
                 harvest_columns_block, is_schema_change, is_doc_problem, banner_line, md
using HDF5, TOML

snapcol(d, s) = Dict{String,Any}("dtype" => d, "shape" => s)
snap(prods) = Dict{String,Any}("release" => "x", "products" =>
    Dict{String,Any}(p => Dict{String,Any}("kind" => "hdf5_table", "n_rows" => 10, "columns" => cols)
                     for (p, cols) in prods))

@testset "dtype and shape normalization" begin
    @test normalize_dtype(Int16) == "int16"
    @test normalize_dtype(Float32) == "float32"
    @test normalize_dtype(Union{Missing,Int64}) == "int64?"
    @test normalize_dtype(HDF5.FixedString{3,0}) == "S3"
    @test normalize_dtype(String) == "string"
    @test shape_string((8700, 26715540), 26715540) == "(8700, N)"
    @test shape_string((10,), 10) == "(N,)"
    @test shape_string((), 10) == "()"
    @test h5py_shape("(8700, N)") == "(N, 8700)"
    @test h5py_shape("(N,)") == "(N,)"
end

@testset "diff_snapshots" begin
    prev = snap(["A" => Dict("x" => snapcol("int64", "(N,)"), "gone" => snapcol("float64", "(N,)"),
                             "s" => snapcol("float64", "(8700, N)")),
                 "Old" => Dict("q" => snapcol("int8", "(N,)"))])
    cur = snap(["A" => Dict("x" => snapcol("int32", "(N,)"), "new" => snapcol("float32", "(N,)"),
                            "s" => snapcol("float64", "(8575, N)")),
                "Fresh" => Dict("r" => snapcol("uint8", "(N,)"))])
    ch = diff_snapshots(prev, cur)
    kinds = Set((c.kind, c.product, c.column) for c in ch)
    @test (:removed_column, "A", "gone") in kinds
    @test (:new_column, "A", "new") in kinds
    @test (:dtype_changed, "A", "x") in kinds
    @test (:shape_changed, "A", "s") in kinds
    @test (:removed_product, "Old", "") in kinds
    @test (:removed_column, "Old", "q") in kinds
    @test (:new_product, "Fresh", "") in kinds
    @test (:new_column, "Fresh", "r") in kinds
    @test all(is_schema_change, ch)
    @test isempty(diff_snapshots(cur, cur))
    # row-count changes are not schema changes
    cur2 = deepcopy(cur); cur2["products"]["A"]["n_rows"] = 99
    @test isempty(diff_snapshots(cur, cur2))
    # the loud wording the release process relies on
    rm = only(filter(c -> c.kind === :removed_column && c.column == "gone", ch))
    @test banner_line(rm) == "COLUMN REMOVED: A.gone — are we sure we don't want this column?"
    nw = DataModel.Change(:new_column, "A", "new"; new = "float32", note = "datamodel/registry/A.toml")
    @test banner_line(nw) == "NEW COLUMN: A.new (float32) — needs a description in datamodel/registry/A.toml"
end

@testset "find_previous_snapshot" begin
    mktempdir() do d
        for v in ("0.3.0", "0.3.10", "0.4.0", "0.5.0")
            write(joinpath(d, v * ".toml"), "release = \"$v\"\n")
        end
        write(joinpath(d, "notes.toml"), "")
        @test basename(find_previous_snapshot(d, "0.4.0")) == "0.3.10.toml"
        @test basename(find_previous_snapshot(d, "0.5.0")) == "0.4.0.toml"
        @test find_previous_snapshot(d, "0.3.0") === nothing
    end
end

@testset "harvesters" begin
    t = """
    | column | type | what it is |
    |---|---|---|
    | `v_rad_mean` | Float64 | unweighted mean, km/s |
    | `ra`, `dec` | Float64 | ICRS position, deg |
    | not a column | x | y |
    """
    h = harvest_markdown_tables(t)
    @test h["v_rad_mean"] == "unweighted mean, km/s"
    @test h["ra"] == "ICRS position, deg" && h["dec"] == "ICRS position, deg"
    @test !haskey(h, "not")
    b = "intro\nCOLUMNS\n  tele, fiber      telescope and fiber\n                   index 1..300\n  era              1 = plate\n\nafter"
    hb = harvest_columns_block(b)
    @test hb["tele"] == "telescope and fiber index 1..300"
    @test hb["fiber"] == hb["tele"]
    @test hb["era"] == "1 = plate"
end

@testset "registry resolution and checks" begin
    p = ProductSchema("P", "hdf5_table", 5, 0, [
            ColumnSchema("a", "int64", "(N,)"),
            ColumnSchema("b", "float64", "(N,)", Dict{String,Any}("description" => "from attr"), Dict{String,Any}()),
            ColumnSchema("c", "float64", "(N,)"),
            ColumnSchema("d", "uint8", "(N,)", Dict{String,Any}("categories" => ["", "x", "y"]), Dict{String,Any}()),
        ], Dict{String,Any}(), Dict{String,Any}())
    reg = Dict{String,Any}("product" => Dict(), "columns" => Dict(
        "a" => Dict("description" => "curated", "units" => "km/s"),
        "d" => Dict("description" => "coded"),
        "zombie" => Dict("description" => "no longer in the file")))
    harv = Dict("c" => "from docs", "a" => "ignored")
    @test resolve_column("a", reg, p.columns[1], harv).source === :registry
    @test resolve_column("a", reg, p.columns[1], harv).description == "curated"
    @test resolve_column("b", reg, p.columns[2], harv).source === :attribute
    @test resolve_column("c", reg, p.columns[3], harv).source === :harvested
    @test resolve_column("c", reg, p.columns[3], Dict{String,String}()).source === :undocumented
    vals = resolve_column("d", reg, p.columns[4], harv).values
    @test [v["value"] for v in vals] == [0, 1, 2] && vals[2]["meaning"] == "x"
    ch = registry_check(p, reg, Dict{String,String}(); registry_file = "datamodel/registry/P.toml")
    @test Set((c.kind, c.column) for c in ch) == Set([(:undocumented, "c"), (:registry_orphan, "zombie")])
    @test all(is_doc_problem, ch)
end

@testset "sanitize and markdown" begin
    pm = Dict("/mnt/ceph/users/sdssv/work/scifest/0.4.0" => "data/work/mwm/0.4.0")
    @test sanitize("see /mnt/ceph/users/sdssv/work/scifest/0.4.0/public/x.h5", pm) == "see data/work/mwm/0.4.0/public/x.h5"
    @test sanitize("raw at /mnt/ceph/users/someone/raw_1/; done", pm) == "raw at [internal path]; done"
    h = md("a `b` <c>\n\n- one\n- two"; path_map = pm)
    @test occursin("<code>b</code>", h) && occursin("&lt;c&gt;", h) && occursin("<li>one</li>", h)
end

# End-to-end on a synthetic release: initial run, then a release that drops one column,
# adds one and retypes one. Exit codes must make unacknowledged changes impossible to miss.
@testset "run_release end to end" begin
    mktempdir() do d
        h5 = joinpath(d, "t.h5")
        h5open(h5, "w") do f
            f["keep"] = collect(1:4)
            f["drop_me"] = rand(4)
            f["retype"] = Int16[1, 2, 3, 4]
        end
        regdir = joinpath(d, "registry"); mkpath(regdir)
        write(joinpath(regdir, "T.toml"), """
        [product]
        general_description = "test"
        [columns.keep]
        description = "kept"
        [columns.drop_me]
        description = "dropped later"
        [columns.retype]
        description = "retyped later"
        """)
        cfg(v) = (p = joinpath(d, "r$v.toml");
                  write(p, "release = \"$v\"\n[[products]]\nname = \"T\"\nkind = \"hdf5_table\"\nsource = \"$h5\"\nregistry = \"T.toml\"\n");
                  p)
        snapdir = joinpath(d, "snaps")
        io = IOBuffer()
        code, ch = run_release(cfg("1.0.0"); outdir = joinpath(d, "html1"), snapshot_dir = snapdir,
                               registry_dir = regdir, io = io)
        @test code == 0
        @test isfile(joinpath(snapdir, "1.0.0.toml"))
        @test isfile(joinpath(d, "html1", "index.html")) && isfile(joinpath(d, "html1", "T.html"))

        h5open(h5, "w") do f
            f["keep"] = collect(1:5)
            f["retype"] = Int32[1, 2, 3, 4, 5]
            f["brand_new"] = rand(5)
        end
        io = IOBuffer()
        code, ch = run_release(cfg("1.1.0"); outdir = joinpath(d, "html2"), snapshot_dir = snapdir,
                               registry_dir = regdir, io = io)
        out = String(take!(io))
        @test code == 2
        @test occursin("COLUMN REMOVED: T.drop_me — are we sure we don't want this column?", out)
        @test occursin("NEW COLUMN: T.brand_new (float64) — needs a description in datamodel/registry/T.toml", out)
        @test occursin("DTYPE CHANGED: T.retype  int16 -> int32", out)
        @test occursin("UNDOCUMENTED COLUMN: T.brand_new", out)
        @test occursin("REGISTRY ENTRY WITHOUT COLUMN: T.drop_me", out)
        html = read(joinpath(d, "html2", "index.html"), String)
        @test occursin("alert red", html) && occursin("COLUMN REMOVED:", html)
        @test occursin("UNDOCUMENTED", read(joinpath(d, "html2", "T.html"), String))

        # acknowledging schema changes still leaves the documentation gap -> 3
        code, _ = run_release(cfg("1.1.0"); outdir = joinpath(d, "html3"), snapshot_dir = snapdir,
                              registry_dir = regdir, accept_changes = true, io = IOBuffer())
        @test code == 3
        code, _ = run_release(cfg("1.1.0"); outdir = joinpath(d, "html4"), snapshot_dir = snapdir,
                              registry_dir = regdir, accept_changes = true, allow_undocumented = true, io = IOBuffer())
        @test code == 0
    end
end
