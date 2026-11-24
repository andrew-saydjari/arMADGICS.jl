# 03_write_allvisitrv.jl -- write the 21-column allVisitRV file from the columns computed by
# 02_build_visits.jl: chunk 65,536 rows, SHUFFLE + gzip 9 on every dataset, `tele` as fixed S3.
#
# CHAIN POSITION: allVisitRV step 3 of 4 (then 04_attrs_allvisitrv.py; the fiber-offset and
#                 offtarget steps later add 4 more columns).
# INPUTS:  ARGS[1]  v2_cols.jls from step 2
# OUTPUTS: ARGS[2]  allVisitRV h5 (21 datasets, no attributes yet)
# SOURCE:  work/2026_09_16/allvisitrv/s10_write_v2.jl (0.4.0). The 0.4.0 script wrote four
#          candidate files (deflate 4/6/9, with/without shuffle), kept the smallest and deleted
#          the rest; the one kept was deflate 9 + shuffle, so this copy writes only that one
#          (writefile() itself unchanged). Paths became ARGS.
using HDF5, Printf, Serialization, Statistics
say(s) = (println(s); flush(stdout))

const C = deserialize(ARGS[1])

# column order: the released file's logical grouping, with v_rad_scatter_jack inserted next to
# the plain v_rad_scatter it debiases.  21 datasets.
cols() = Any[
    ("sdss_id",            C.O_sid),
    ("mjd",                C.O_mjd),
    ("mjd_mid",            C.O_mid),
    ("v_rad_mean",         C.O_v),
    ("n_exp",              C.O_n),
    ("snr_mean",           C.O_snrm),
    ("snr_visit",          C.O_snrq),
    ("v_rad_scatter",      C.O_sc),
    ("v_rad_scatter_jack", C.O_jk),
    ("e_v_rad_visit",      C.O_ev),
    ("ra",                 C.O_ra),
    ("dec",                C.O_dec),
    ("glon",               C.O_glon),
    ("glat",               C.O_glat),
    ("fiber",              C.O_fib),
    ("era",                C.O_era),
    ("field_id",           C.O_fld),
    ("plate_id",           C.O_plate),
    ("config_id",          C.O_config),
    ("category",           C.O_cat),
]

function writefile(path, level, shuffle)
    isfile(path) && rm(path)
    n = length(C.O_sid); ck = (min(1 << 16, n),)
    h5open(path, "w") do g
        for (nm, v) in cols()
            g[nm, chunk = ck, shuffle = shuffle, deflate = level] = v
        end
        dt = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
        HDF5.API.h5t_set_size(dt, 3)
        HDF5.API.h5t_set_cset(dt, HDF5.API.H5T_CSET_ASCII)
        buf = Vector{UInt8}(undef, 3n)
        for i in 1:n
            b = codeunits(C.telestr[i]); buf[3i-2] = b[1]; buf[3i-1] = b[2]; buf[3i] = b[3]
        end
        ds = create_dataset(g, "tele", dt, dataspace((n,)); chunk = ck, shuffle = shuffle, deflate = level)
        write_dataset(ds, dt, buf)
    end
    chmod(path, 0o644)
    return filesize(path)
end

function main()
    out = ARGS[2]
    t = time()
    sz = writefile(out, 9, true)
    @printf("wrote %s  deflate=9 shuffle=true : %10.2f MB   (%.0f s)\n", out, sz/1e6, time()-t)
end

main()
