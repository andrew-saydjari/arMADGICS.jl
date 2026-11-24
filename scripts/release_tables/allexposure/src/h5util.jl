# h5util.jl -- small HDF5 helpers: layout-faithful dataset writes and typed string attributes.

chunk_rows(::Type{T}) where {T} = sizeof(T) >= 4 ? 13045 : sizeof(T) == 2 ? 26090 : 52179

"Write vector `v` as dataset `name` with one of the schema layout classes."
function write_column(f::HDF5.File, name::AbstractString, v::AbstractVector{T}, layout::Symbol) where {T}
    n = length(v)
    if layout === :contig
        f[name] = v
        return f[name]
    end
    c = layout === :gz13045 ? 13045 : chunk_rows(T)
    c = min(c, max(n, 1))
    d = if layout === :gz_nosh
        create_dataset(f, name, datatype(T), dataspace((n,)); chunk = (c,), deflate = 1)
    else
        create_dataset(f, name, datatype(T), dataspace((n,)); chunk = (c,), shuffle = true, deflate = 1)
    end
    n > 0 && write(d, v)
    return d
end

"Fixed-length NULLPAD ASCII string dataset, as numpy's `|S<n>` arrays are stored by h5py."
function write_fixed_strings(f, name::AbstractString, v::Vector{String})
    n = max(1, maximum(ncodeunits, v; init = 0))
    dt = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
    HDF5.API.h5t_set_size(dt, n)
    HDF5.API.h5t_set_strpad(dt, HDF5.API.H5T_STR_NULLPAD)
    buf = zeros(UInt8, n * length(v))
    for (i, s) in enumerate(v)
        cu = codeunits(s); buf[(i-1)*n+1:(i-1)*n+length(cu)] .= cu
    end
    d = create_dataset(f, name, dt, dataspace((length(v),)))
    write_dataset(d, dt, buf)
    close(d); close(dt)
    return nothing
end

function _vlen_utf8()
    dt = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
    HDF5.API.h5t_set_size(dt, HDF5.API.H5T_VARIABLE)
    HDF5.API.h5t_set_cset(dt, HDF5.API.H5T_CSET_UTF8)
    return dt
end

function _put_vlen!(obj, name, strs::Vector{String}, scalar::Bool)
    dt = _vlen_utf8()
    sp = scalar ? dataspace(()) : dataspace((length(strs),))
    a = create_attribute(obj, name, dt, sp)
    GC.@preserve strs begin
        ptrs = [Base.unsafe_convert(Cstring, s) for s in strs]
        HDF5.API.h5a_write(a, dt, ptrs)
    end
    close(a); close(dt)
    return nothing
end

"Write an attribute with the exact HDF5 string flavour the 0.4.0 file used."
put_attr!(obj, name, x::FixedStr) = (attributes(obj)[name] = x.s; nothing)
put_attr!(obj, name, x::VlenStr) = _put_vlen!(obj, name, [x.s], true)
put_attr!(obj, name, x::VlenStrVec) = _put_vlen!(obj, name, x.v, false)
put_attr!(obj, name, x::Integer) = (attributes(obj)[name] = Int64(x); nothing)
