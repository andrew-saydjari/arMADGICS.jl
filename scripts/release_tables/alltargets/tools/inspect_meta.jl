# Print the datasets of /meta (or another group) with dtype, shape and attributes. Read-only.
# Usage: julia --project=. tools/inspect_meta.jl FILE [GROUP]
using HDF5
path = ARGS[1]; grp = length(ARGS) > 1 ? ARGS[2] : "meta"
h5open(path, "r") do f
    println("top-level: ", keys(f))
    println("file attrs: ", Dict(k => read(attributes(f)[k]) for k in keys(attributes(f))))
    g = f[grp]
    println("group attrs: ", Dict(k => (v = read(attributes(g)[k]); v isa AbstractArray && length(v) > 10 ? "$(typeof(v)) len $(length(v))" : v) for k in keys(attributes(g))))
    for (i, k) in enumerate(keys(g))   # keys(): HDF5 iteration order (creation order if tracked, else name)
        d = g[k]
        d isa HDF5.Dataset || (println(i, " ", k, " GROUP"); continue)
        dt = datatype(d)
        cls = HDF5.API.h5t_get_class(dt)
        vl = cls == HDF5.API.H5T_STRING ? (HDF5.API.h5t_is_variable_str(dt) ? "vlen-str" : "fixed-str$(HDF5.API.h5t_get_size(dt))") : string(cls)
        a = Dict(kk => read(attributes(d)[kk]) for kk in keys(attributes(d)))
        println(i, " ", k, " ", vl, " size=", sizeof(dt), " shape=", size(d), " attrs=", a)
    end
end
