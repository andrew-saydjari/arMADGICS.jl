## Schema types shared by introspection, snapshots and diffs

"""
    ColumnSchema

One column (HDF5 dataset or CSV column) as found in the file.

- `dtype`: normalized type name (`int64`, `float32`, `uint8`, `S3`, `string`, `bool`, ...);
  CSV columns that contain empty cells carry a `?` suffix (`float64?`).
- `shape`: Julia (column-major) shape with the row axis written as `N`, e.g. `(8700, N)`.
  h5py reports the reverse order, `(N, 8700)`.
- `attrs`: HDF5 dataset attributes (empty for CSV).
- `extra`: per-column facts that are not schema (file size, presence in the release dir).
"""
struct ColumnSchema
    name::String
    dtype::String
    shape::String
    attrs::Dict{String,Any}
    extra::Dict{String,Any}
end
ColumnSchema(name, dtype, shape) = ColumnSchema(name, dtype, shape, Dict{String,Any}(), Dict{String,Any}())

"""
    ProductSchema

A released product (one file, or a family of per-key files) as found on disk.
`kind` is one of `hdf5_table`, `hdf5_perkey`, `csv`, `reference`.
"""
struct ProductSchema
    name::String
    kind::String
    n_rows::Int
    size_bytes::Int
    columns::Vector{ColumnSchema}
    file_attrs::Dict{String,Any}
    extra::Dict{String,Any}
end

column_names(p::ProductSchema) = [c.name for c in p.columns]

function getcolumn(p::ProductSchema, name::AbstractString)
    i = findfirst(c -> c.name == name, p.columns)
    return i === nothing ? nothing : p.columns[i]
end

"""
    Change

One finding from a diff. `kind` is one of

- `:removed_product`, `:new_product`
- `:removed_column`, `:new_column`, `:dtype_changed`, `:shape_changed`  (schema changes vs previous release)
- `:undocumented`        column present in the files with no description from any source
- `:registry_orphan`     registry entry for a column that is not in the files
- `:initial`             no previous snapshot; informational
"""
struct Change
    kind::Symbol
    product::String
    column::String
    old::String
    new::String
    note::String
end
Change(kind, product, column; old = "", new = "", note = "") =
    Change(kind, String(product), String(column), String(old), String(new), String(note))

const SCHEMA_CHANGE_KINDS = (:removed_product, :new_product, :removed_column, :new_column,
                             :dtype_changed, :shape_changed)
const DOC_CHANGE_KINDS = (:undocumented, :registry_orphan)

is_schema_change(c::Change) = c.kind in SCHEMA_CHANGE_KINDS
is_doc_problem(c::Change) = c.kind in DOC_CHANGE_KINDS
