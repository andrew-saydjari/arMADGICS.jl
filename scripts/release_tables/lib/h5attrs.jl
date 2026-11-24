# h5attrs.jl -- small shared HDF5 helpers for the release-table helper scripts (include()d, not a
# package): attributes with an exact HDF5 string flavour, a raw type-preserving attribute copy, and
# the attribute "flavour" record that validation/compare_h5.{py,jl} compare.
#
# Why: h5py writes a Python `str` attribute as a SCALAR, VARIABLE-length, UTF-8, NULLTERM string,
# whereas HDF5.jl's `attributes(obj)[name] = "..."` writes a fixed-length string. The release
# files carry h5py's flavour, so the Julia ports write it explicitly. Same approach as
# allexposure/src/h5util.jl (_vlen_utf8 / _put_vlen!), copied here so that the helper scripts do
# not depend on the AllExposure package.
using HDF5

"Variable-length, UTF-8, NULLTERM C string type (what h5py uses for `str`). Caller closes it."
function vlen_utf8_type()
    dt = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
    HDF5.API.h5t_set_size(dt, HDF5.API.H5T_VARIABLE)
    HDF5.API.h5t_set_cset(dt, HDF5.API.H5T_CSET_UTF8)
    return dt
end

"Write `strs` as a variable-length UTF-8 string attribute: scalar (one string) or 1-D array."
function put_vlen_attr!(obj, name::AbstractString, strs::Vector{String}; scalar::Bool = false)
    scalar && length(strs) != 1 && error("scalar attribute needs exactly one string")
    HDF5.API.h5a_exists(obj, name) && HDF5.API.h5a_delete(obj, name)   # h5py overwrites silently
    dt = vlen_utf8_type()
    sp = scalar ? dataspace(()) : dataspace((length(strs),))
    a = create_attribute(obj, name, dt, sp)
    GC.@preserve strs begin
        ptrs = [Base.unsafe_convert(Cstring, s) for s in strs]
        HDF5.API.h5a_write(a, dt, ptrs)
    end
    close(a); close(sp); close(dt)
    return nothing
end

"h5py's `obj.attrs[name] = <python str>`."
put_str_attr!(obj, name::AbstractString, s::AbstractString) = put_vlen_attr!(obj, name, [String(s)]; scalar = true)

"""
    copy_attr_raw!(src, name, dst, dstname = name)

Copy attribute `name` of `src` to `dst` with the SAME HDF5 datatype and dataspace (fixed or
variable-length strings, numeric, any shape): the bytes are read in the file type and written back
in it, then any variable-length memory HDF5 allocated is reclaimed.
"""
function copy_attr_raw!(src, name::AbstractString, dst, dstname::AbstractString = name)
    A = HDF5.API
    aid = A.h5a_open(src, name, A.H5P_DEFAULT)
    tid = A.h5a_get_type(aid)
    sid = A.h5a_get_space(aid)
    npts = prod(Int.(A.h5s_get_simple_extent_dims(sid, nothing)); init = 1)
    buf = zeros(UInt8, max(1, npts) * A.h5t_get_size(tid))
    A.h5a_read(aid, tid, buf)
    A.h5a_exists(dst, dstname) && A.h5a_delete(dst, dstname)
    bid = A.h5a_create(dst, dstname, tid, sid, A.H5P_DEFAULT, A.H5P_DEFAULT)
    A.h5a_write(bid, tid, buf)
    A.h5a_close(bid)
    if A.h5t_get_class(tid) == A.H5T_VLEN ||
       (A.h5t_get_class(tid) == A.H5T_STRING && A.h5t_is_variable_str(tid))
        A.h5d_vlen_reclaim(tid, sid, A.H5P_DEFAULT, buf)
    end
    A.h5s_close(sid); A.h5t_close(tid); A.h5a_close(aid)
    return nothing
end

"""
    attr_flavour(obj, name) -> Vector{Pair{String,Any}}

The record compare_h5.py's `attr_flavour` builds, in the same key order: HDF5 type class and
dataspace shape (C order), then for strings vlen / strpad / cset (and size if fixed), else the
type size.
"""
function attr_flavour(obj, name::AbstractString)
    A = HDF5.API
    aid = A.h5a_open(obj, name, A.H5P_DEFAULT)
    tid = A.h5a_get_type(aid)
    sid = A.h5a_get_space(aid)
    shape = A.h5s_get_simple_extent_type(sid) == A.H5S_NULL ? Int[] :
            Int.(A.h5s_get_simple_extent_dims(sid, nothing))
    d = Pair{String,Any}["class" => Int(A.h5t_get_class(tid)), "shape" => shape]
    if A.h5t_get_class(tid) == A.H5T_STRING
        vl = A.h5t_is_variable_str(tid)
        push!(d, "vlen" => vl, "pad" => Int(A.h5t_get_strpad(tid)), "cset" => Int(A.h5t_get_cset(tid)))
        vl || push!(d, "size" => Int(A.h5t_get_size(tid)))
    else
        push!(d, "size" => Int(A.h5t_get_size(tid)))
    end
    A.h5s_close(sid); A.h5t_close(tid); A.h5a_close(aid)
    return d
end
