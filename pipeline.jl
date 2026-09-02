## This is the main pipeline that will batch over APOGEE files
# Author - Andrew Saydjari, CfA

import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
Pkg.update("ApogeeReduction");
Pkg.instantiate();
Pkg.precompile();
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Package activation took $dt");
t_then = t_now;
flush(stdout);
using BLISBLAS
using Distributed, ArgParse, SlurmClusterManager, Suppressor, DataFrames, DelimitedFiles
using ApogeeReduction: read_almanac_exp_df, get_fibTargDict, check_type_for_jld2

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--redux_base"
        required = true
        help = "redux base directory"
        arg_type = String
        "--almanac_file"
        required = true
        help = "almanac file"
        arg_type = String
        "--outdir"
        required = false
        help = "output directory"
        arg_type = String
        default = "../outdir/arMADGICS/raw/"
        "--local_cache_dir"
        required = false
        help = "local cache directory"
        arg_type = String
        default = "../local_cache/"
    end
    return parse_args(s)
end
parg = parse_commandline()

proj_path = dirname(Base.active_project()) * "/"
if "SLURM_NTASKS" in keys(ENV)
    using SlurmClusterManager
    addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
else
    addprocs(32, exeflags = ["--project=$proj_path"]) # change to a workers per node variable
end
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker allocation took $dt");
t_then = t_now;
flush(stdout);
println("Running Main on ", gethostname());
flush(stdout);

@everywhere begin
    using BLISBLAS
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, JLD2, FileIO
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using ThreadPinning, ApogeeReduction, DataFrames
    using StatsBase, ProgressMeter
    using ApogeeReduction: check_type_for_jld2
end
@passobj 1 workers() parg
@passobj 1 workers() proj_path
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker loading took $dt");
t_then = t_now;
flush(stdout);
println(BLAS.get_config());
flush(stdout);

@everywhere begin
    prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
    src_dir = "$proj_path"
    include(joinpath(src_dir, "src/utils.jl"))
    include(joinpath(src_dir, "src/gridSearch.jl"))
    include(joinpath(src_dir, "src/componentAndPosteriors.jl"))
    include(joinpath(src_dir, "src/fileNameHandling.jl"))
    include(joinpath(src_dir, "src/ingest.jl"))
    include(joinpath(src_dir, "src/lowRankPrescription.jl"))
    include(joinpath(src_dir, "src/marginalizeEW.jl"))
    include(joinpath(src_dir, "src/spectraInterpolation.jl"))
    include(joinpath(src_dir, "src/chi2Wrappers.jl"))
    include(joinpath(src_dir, "src/pipelineCore.jl"))
end

using LibGit2;
println(proj_path)
git_branch, git_commit, git_clean = initalize_git(proj_path);
@passobj 1 workers() git_branch;
@passobj 1 workers() git_commit;
@passobj 1 workers() git_clean;

# These global allocations for the injest are messy... but we plan on changing the ingest
# relatively soon... so don't worry for now.
@everywhere begin
    refine_iters = 5
    ddstaronly = false
    checkpoint_mode = "commit_exists"
    runlist_range = 1:600 # 295, 245, 335, 101
    batchsize = 100 #10

    # Step Size for Chi2 Surface Error Bars
    RV_err_step = 4
    # DIB_pix_err_step = 3 # consider increasing to 4 (self consistency + LSF test)
    # DIB_sig_err_step = 3
    # Flux marginalize region
    # sigMarg0 = -50//100:10//100:50//100
    # svalMarg0 = -0//10:1//10:0//10;

    out_dir = parg["outdir"]
    cache_dir = parg["local_cache_dir"]
    inject_cache_dir = replace(cache_dir, "local_cache/" => "inject_local_cache/")

    reduxBase = parg["redux_base"]
    almanacFile = parg["almanac_file"]

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Sky Priors
    prior_dict["skycont"] = prior_dir * "2025_07_31/prior_dump/APOGEE_skycont_svd_30_f"
    prior_dict["skyLines_bright"] = prior_dir * "2025_07_31/prior_dump/sky_priors/APOGEE_skyline_bright_GSPICE_svd_120_f"
    prior_dict["skyLines_faint"] = prior_dir * "2025_07_31/prior_dump/sky_priors/APOGEE_skyline_faint_GSPICE_svd_120_f"

    # Star Priors
    # prior_dict["starCont"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/star_priors/APOGEE_starcont_svd_60_f"
    prior_dict["chebmsk"] = prior_dir * "2025_07_31/prior_dump/chebmsk_exp.h5"
    prior_dict["starCont"] = prior_dir * "2025_07_31/prior_dump/APOGEE_starcont_svd_60_rough.h5"
    prior_dict["starLines_refLSF"] = prior_dir * "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5"
    # prior_dict["starLines_LSF"] = prior_dir*"2024_03_16/arMADGICS.jl/src/prior_build/starLine_priors_norm94_dd/APOGEE_starCor_svd_50_subpix_f" # DD Version
    # prior_dict["starLines_LSF"] = prior_dir*"2024_02_21/arMADGICS.jl/src/prior_build/starLine_priors_norm94/APOGEE_stellar_kry_50_subpix_f" # TH Version

    # # DIB Priors
    # dib_waves = [15273, 15672]
    # for dib in dib_waves
    #     # prior_dict["DIB_noLSF_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_1_$(dib)_analyticDeriv_stiff.h5"
    #     scratch1/working/
    #     prior_dict["DIB_noLSF_soft_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_3_$(dib)_analyticDeriv_soft.h5"
    #     prior_dict["DIB_LSF_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_1_$(dib)_analyticDerivLSF_stiff_"
    #     prior_dict["DIB_LSF_soft_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_3_$(dib)_analyticDerivLSF_soft_"
    # end

    # # maps dib scans to dib_waves
    # dib_ind_prior = Dict{Int,Int}()
    # dib_ind_prior[1] = 1
    # dib_ind_prior[2] = 1
    # dib_ind_prior[3] = 2
    # dib_ind_prior[4] = 2
end

# it would be great to move this into a parameter file that is read for each run
# I should revisit the error bars in the context of chi2 versus frame number trends
@everywhere begin
    # minimum step sizes of the priors (might want to store in the file structures of those priors and read in)
    minRVres = 1 // 10
    # minDIBvelres = 1//10
    # minDIBsigres = 1//100

    ## sigrng is somewhat counterintuitively defined in the chi2Wrappers.jl file
    minSigval, maxSigval = extrema(sigrng)

    # Star Wave
    lvl1 = -70:1//2:70
    lvl2 = -8:2//10:8
    lvl3 = -3:1//10:3
    slvl_tuple = (lvl1, lvl2, lvl3)
    # tuple1dprint(slvl_tuple)

    # # (Wave, Sig) DIB
    # dib_center_lst = map(x->dib_waves[dib_ind_prior[x]],1:length(dib_ind_prior)) # not clear we need this anymore since we don't scanOffset
    # lvl1d_15273wide = ((-137:4:150),(18//10:18//10))
    # lvl1d_15672wide = ((-54:4:150),(18//10:18//10))
    # lvl1d_narrow = ((-54:4:54),(18//10:18//10))
    # lvl1d = ((-150:4:150),(18//10:18//10))
    # lvl2d = ((0:0), (-7//5:4//100:11//5))
    # lvl3d = ((-18:2//10:18), (0:0))
    # lvl4d = ((0:0), (-90//100:2//100:90//100))
    # lvl5d = ((-1:2//10:1), (0:0))
    # lvl6d = ((0:0), (-10//100:2//100:10//100))
    # lvl7d = ((-6//10:2//10:6//10), (0:0))
    # lvl8d = ((0:0), (-6//100:2//100:6//100))
    # lvl9d = ((-4//10:1//10:4//10), (-4//100:1//100:4//100));
    # lvltuple = (lvl1d, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_15273wide = (lvl1d_15273wide, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_15672wide = (lvl1d_15672wide, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_narrow = (lvl1d_narrow, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_lst = [lvltuple_15273wide, lvltuple_narrow, lvltuple_15672wide, lvltuple_narrow]
    # # tuple2dprint(lvltuple)
end

@everywhere begin
    logUniWaveAPOGEE = 10 .^ range((start = 4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)
    minw, maxw = extrema(logUniWaveAPOGEE)

    c = 299792.458 # in km/s
    delLog = 6e-6

    Xd_stack = zeros(3 * 2048)
    Xd_std_stack = zeros(3 * 2048)
    waveobs_stack = zeros(3 * 2048)
    waveobs_stack_old = zeros(3 * 2048)
    pixmsk_stack = zeros(Int, 3 * 2048)
    fullBit = zeros(Int, 3 * 2048)
    outvec = zeros(length(logUniWaveAPOGEE))
    outvar = zeros(length(logUniWaveAPOGEE))
    cntvec = zeros(Int, length(logUniWaveAPOGEE))
end


@everywhere begin
    function multi_spectra_batch(indsubset; out_dir=out_dir, ddstaronly=ddstaronly, checkpoint_mode=checkpoint_mode)
        ### Set up
        out = []
        startind = indsubset[1].linear_index
        adjfiberindx = indsubset[1].adjfiberindx

        ### Save and cache restart handling
        savename = join([out_dir, lpad(adjfiberindx, 3, "0"), "arMADGICS_fiber_" * lpad(adjfiberindx, 3, "0") * "_batch_" * lpad(startind, 7, "0") * ".h5"], "/")
        dirName = splitdir(savename)[1]
        if !ispath(dirName)
            mkpath(dirName)
        end
        # probably should shift this to a check_file function like in ApogeeReduction.jl
        # if check_file(savename, mode = checkpoint_mode)
        if !isfile(savename)
            # We are loading the priors EVERY time, so there is no benefit to ordering
            # This is not optimal, but reduces scope confusion
            # These first two could be globals, but load them here for consistency
            # V_dib_noLSF_soft_lst = []
            # for dib in dib_waves
            #     local f = h5open(prior_dict["DIB_noLSF_soft_$(dib)"])
            #     push!(V_dib_noLSF_soft_lst, read(f["Vmat"]))
            #     close(f)
            # end

            f = h5open(prior_dict["starLines_refLSF"])
            V_starlines_refLSF = read(f["Vmat"])
            close(f)

            ### Need to load the priors here
            # f = h5open(prior_dict["skycont"]*lpad(adjfiberindx ,3,"0")*".h5")
            # V_skycont = read(f["Vmat"])
            # chebmsk_exp = convert.(Bool,read(f["chebmsk_exp"]))
            # close(f)

            # f = h5open(prior_dict["skyLines_bright"]*lpad(adjfiberindx ,3,"0")*".h5")
            # V_skyline_bright = read(f["Vmat"])
            # submsk_bright = convert.(Bool,read(f["submsk"]))
            # close(f)

            # f = h5open(prior_dict["skyLines_faint"]*lpad(adjfiberindx ,3,"0")*".h5")
            # V_skyline_faint = read(f["Vmat"])
            # submsk_faint = convert.(Bool,read(f["submsk"]))
            # close(f)

            # skymsk_bright = chebmsk_exp .& submsk_bright #
            # skymsk_faint = chebmsk_exp .& submsk_faint
            # # global skymsk = chebmsk_exp .& (submsk_bright .| submsk_faint)
            # skymsk = chebmsk_exp .& submsk_faint # completely masking all bright lines b/c detector response is nonlinear;

            f = h5open(prior_dict["chebmsk"])
            chebmsk_exp = read(f["chebmsk_exp"])
            close(f)
            skymsk_bright = chebmsk_exp #.& submsk_bright #
            skymsk_faint = chebmsk_exp #.& submsk_faint
            skymsk = chebmsk_exp #.& submsk_faint # completely masking all bright lines b/c detector response is nonlinear;

            # f = h5open(prior_dict["starCont"]*lpad(adjfiberindx ,3,"0")*".h5")
            f = h5open(prior_dict["starCont"])
            V_starcont = read(f["Vmat"])
            close(f)

            V_starlines = V_starlines_refLSF #hack
            # f = h5open(prior_dict["starLines_LSF"]*lpad(adjfiberindx ,3,"0")*".h5")
            # V_starlines = read(f["Vmat"])
            if ddstaronly
                V_starlines_refLSF = V_starlines
                msk_starCor = convert.(Bool, read(f["msk_starCor"]))
            else
                msk_starCor = ones(Bool, length(chebmsk_exp))
            end
            # close(f)

            # V_dib_lst = []
            # for dib in dib_waves
            #     local f = h5open(prior_dict["DIB_LSF_$(dib)"]*lpad(adjfiberindx ,3,"0")*".h5")
            #     push!(V_dib_lst,read(f["Vmat"]))
            #     close(f)
            # end

            # V_dib_soft_lst = []
            # for dib in dib_waves
            #     local f = h5open(prior_dict["DIB_LSF_soft_$(dib)"]*lpad(adjfiberindx ,3,"0")*".h5")
            #     push!(V_dib_soft_lst,read(f["Vmat"]))
            #     close(f)
            # end

            ### Single spectrum loop
            # prior_vec = (V_skycont,chebmsk_exp,V_skyline_bright,V_skyline_faint,skymsk_bright,skymsk_faint,skymsk,V_starcont,V_starlines_refLSF,V_starlines,msk_starCor,V_dib_lst, V_dib_soft_lst,V_dib_noLSF_soft_lst)
            prior_vec = (chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont, V_starlines_refLSF, V_starlines, msk_starCor)
            pipeline_single_spectra_bind(argtup) = pipeline_single_spectra(argtup, prior_vec; ddstaronly=ddstaronly)
            for (ind, indval) in enumerate(indsubset)
                push!(out, pipeline_single_spectra_bind(indval))
            end

            ### Save Exporting
            metai = 1
            RVind, RVchi, RVcom, strpo = 2, 3, 4, 5
            # DIBind, EWind, DIBchi, DIBcom = 6, 7, 8, 9
            # dibsavesz = 4

            ## RV Block
            RVextract = [
                # meta info
                (x -> x[metai][1], "data_pix_cnt"),
                (x -> x[metai][2], "starscale"),
                (x -> x[metai][3], "skyscale0"),
                # (x->x[metai][4],                        "frame_counts"),
                # (x->x[metai][5],                        "chip_midtimes"),
                # (x->x[metai][6],                        "a_relFlux"),
                # (x->x[metai][7],                        "b_relFlux"),
                # (x->x[metai][8],                        "c_relFlux"),
                # (x->x[metai][9],                        "cartVisit"),
                # (x->x[metai][10],                       "ingestBit"),
                (x -> x[metai][4], "flux"),
                (x -> x[metai][5], "fluxivar"),
                (x -> x[metai][6], "flux_nans"),
                (x -> x[metai][7], "fluxerr2_nans"),
                (x -> convert(Vector{Int}, x[metai][8]), "simplemsk"),
                (x -> x[metai][9], "nSkyFibers"),
                (x -> x[metai][10], "snr"),
                (x -> x[metai][11], "ingestBit"), # M2: per-spectrum ingest/failure code (bits in src/ingest.jl)
                (x -> adjfiberindx, "adjfiberindx"),
                (x -> Float64.(x[RVind][1][1]), "RV_pixoff_final"),
                (x -> Float64.(x[RVind][1][3]), "RV_pixoff_disc_final"),
                (x -> x[RVind][1][2], "RV_minchi2_final"),
                (x -> x[RVind][1][6], "RV_flag"),
                (x -> x[RVind][1][7], "RV_pix_var"), (x -> x[RVchi][1], "RVchi2_residuals"),
                # (x->x[RVchi][2],                        "RVchi2_residuals_flux_scaled"),
                (x -> x[RVchi][2], "avg_flux_conservation"),
                (x -> x[RVchi][3], "final_pix_cnt"),
                (x -> x[RVchi][4], "starscale1"),
                (x -> x[RVchi][5], "skyscale1"), (x -> x[RVind][2][1][3], "RV_p5delchi2_lvl1"),
                (x -> x[RVind][2][2][3], "RV_p5delchi2_lvl2"),
                (x -> x[RVind][2][3][3], "RV_p5delchi2_lvl3"), (x -> x[RVcom][1], "x_residuals_z_v0"),
                (x -> x[RVcom][2], "x_residuals_v0"),
                # (x->x[RVcom][3],                        "x_skyLines_bright_v0"),
                (x -> x[RVcom][3], "x_skyLines_faint_v0"),
                (x -> x[RVcom][4], "x_skyContinuum_v0"),
                (x -> x[RVcom][5], "x_starContinuum_v0"),
                # (x->x[RVcom][6],                        "x_starLineCor_v0"),
                (x -> x[RVcom][6], "x_starLines_v0"),
                (x -> x[RVcom][7], "x_starLineCof_v0"),
                (x -> x[RVcom][8], "tot_p5chi2_v0"),
                (x -> x[RVcom][9], "apVisit_v0"),
                (x -> Int.(x[RVcom][10]), "finalmsk"),
                (x -> x[RVcom][11], "x_starLines_restFrame_v0"), (x -> x[strpo], "x_starLines_err_v0"),
            ]
            # ## DIB Block
            # DIBextract = []
            # for (dibindx,dibw) in enumerate(dib_center_lst)
            #     dib = string(round(Int,dibw))
            #     dibind = string(dibindx)
            #     push!(DIBextract,[
            #     # Further chi2 refinement does not have fixed sizing because can hit grid edge
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][1][1]),        "DIB_pixoff_final_$(dibind)_$(dib)"),
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][1][2]),        "DIB_sigval_final_$(dibind)_$(dib)"),
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][3][1]),        "DIB_pixoff_disc_final_$(dibind)_$(dib)"),
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][3][2]),        "DIB_sigval_disc_final_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][1][2],                     "DIB_minchi2_final_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][1][6],                     "DIB_flag_$(dibind)_$(dib)"),
            #     (x->[x[DIBind+dibsavesz*(dibindx-1)][1][7:11]...],             "DIB_hess_var_$(dibind)_$(dib)"),

            #     (x->x[DIBind+dibsavesz*(dibindx-1)][2][1][3],                  "DIB_p5delchi2_lvl1_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][2][2][3],                  "DIB_p5delchi2_lvl2_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][2][3][3],                  "DIB_p5delchi2_lvl3_$(dibind)_$(dib)"),

            #     (x->x[EWind+dibsavesz*(dibindx-1)][1],                         "EW_dib_$(dibind)_$(dib)"),
            #     (x->x[EWind+dibsavesz*(dibindx-1)][2],                         "EW_dib_err_$(dibind)_$(dib)"),

            #     (x->x[DIBchi+dibsavesz*(dibindx-1)][1],                        "DIBchi2_residuals_$(dibind)_$(dib)"),
            #     # (x->x[DIBchi+dibsavesz*(dibindx-1)][2],                        "DIBchi2_residuals_flux_scaled_$(dibind)_$(dib)"),

            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][1],                        "x_residuals_z_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][2],                        "x_residuals_v1_$(dibind)_$(dib)"),
            #     # (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_bright_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_faint_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][4],                        "x_skyContinuum_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][5],                        "x_starContinuum_v1_$(dibind)_$(dib)"),
            #     # (x->x[DIBcom+dibsavesz*(dibindx-1)][6],                        "x_starLineCor_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][6],                        "x_dib_flux_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][7],                        "x_starLines_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][8],                        "x_dib_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][9],                        "tot_p5chi2_v1_$(dibind)_$(dib)"),
            #     ])
            # end
            # extractlst = vcat(RVextract...,DIBextract...)
            extractlst = vcat(RVextract...)

            for elelst in extractlst
                extractor(out, elelst[1], elelst[2], savename)
            end
            
            hdr_dict = Dict(
                "pipeline" => "arMADGICS.jl",
                "git_branch" => git_branch,
                "git_commit" => git_commit,
                "git_clean" => git_clean,
            )
            h5write(savename, "hdr", "This is only a header")
            h5writeattr(savename, "hdr", hdr_dict)
            # # could merge into safe_jldsave-like handling from AR.jl
            # # add the git info to the metadata
            # metadata = Dict{String, Any}()
            # metadata["pipeline"] = "arMADGICS.jl"
            # metadata["git_branch"] = git_branch
            # metadata["git_commit"] = git_commit
            # metadata["git_clean"] = git_clean

            # # add metadata group to the file
            # h5open(filename, "r+") do f
            #     g = create_group(f, "metadata")
            #     for (k, v) in metadata
            #         g[k] = check_type_for_jld2(v)
            #     end
            # end
        end
        return 0
    end

    function extractor(x, elemap, elename, savename)
        len = length(x)
        exobj = elemap(x[1])
        outmat = zeros(eltype(exobj), size(exobj)..., len)
        for i = 1:len
            flush(stdout)
            outmat[.., i] .= elemap(x[i])
        end
        h5write(savename, elename, outmat)
    end
end

# Collect all (tele, mjd) pairs first
f = h5open(almanacFile)
tele_mjd_pairs = []
if haskey(f, "raw")
    for tele in keys(f["raw"])
        for mjd in keys(f["raw"][tele])
            push!(tele_mjd_pairs, (tele, mjd))
        end
    end
else
    for tele in keys(f)
        for mjd in keys(f[tele])
            push!(tele_mjd_pairs, (tele, mjd))
        end
    end
end
close(f)
@everywhere get_telemjd_runlist_from_almanac_partial(argtup) = get_telemjd_runlist_from_almanac(almanacFile, argtup[1], argtup[2])
run_lsts = pmap(get_telemjd_runlist_from_almanac_partial, tele_mjd_pairs)
run_lst = vcat(run_lsts...)

# this takes like 5 min to run (speed up)
iterlst = []
iterlist_full = []
Base.length(f::Iterators.Flatten) = sum(length, f.it)
for adjfiberindx in 1:600
    subiter = filter(x -> x[:adjfiberindx] .== adjfiberindx, run_lst)
    push!(iterlist_full, subiter)
    indexed_tuples = map(enumerate(subiter)) do (idx, named_tuple)
        merge((linear_index=idx,), named_tuple)
    end
    subiterpart = Iterators.partition(indexed_tuples, batchsize)
    push!(iterlst, subiterpart)
end
ittot = Iterators.flatten(iterlst)
lenargs = length(ittot)
nwork = length(workers())
println("Batches to Do: $lenargs, number of workers: $nwork")
flush(stdout)

# Write flattened iterlist_full to a jld2 file (clean this up later)
sdss_id_lst = []
tele_lst = []
mjd_lst = []
expnum_lst = []
adjfiberindx_lst = []
for iter in iterlist_full
    push!(sdss_id_lst, map(x -> x.sdss_id, iter))
    push!(tele_lst, map(x -> x.tele, iter))
    push!(mjd_lst, map(x -> x.mjd, iter))
    push!(expnum_lst, map(x -> x.expnum, iter))
    push!(adjfiberindx_lst, map(x -> x.adjfiberindx, iter))
end
sdss_id_all = vcat(sdss_id_lst...)
tele_all = vcat(tele_lst...)
mjd_all = vcat(mjd_lst...)
expnum_all = vcat(expnum_lst...)
adjfiberindx_all = vcat(adjfiberindx_lst...);

full_list_info_file = joinpath(out_dir, "full_list_info.h5")
dirName = splitdir(full_list_info_file)[1]
if !ispath(dirName)
    mkpath(dirName)
end
jldsave(full_list_info_file, sdss_id=sdss_id_all, tele=tele_all, mjd=mjd_all, expnum=expnum_all, adjfiberindx=adjfiberindx_all)

# Write the batch information to a simple text file for easy parsing
println("Writing batch information to file...")
batch_info_file = joinpath(out_dir, "batch_info.txt")

open(batch_info_file, "w") do io
    println(io, "# Batch information for arMADGICS pipeline")
    println(io, "# Format: linear_index, tele, mjd, expnum, adjfiberindx ")
    println(io, "# Total batches: $lenargs")
    println(io, "# Generated on: $(now())")
    println(io, "")

    for batch in ittot
        for item in batch
            println(io, "$(item.linear_index), $(item.tele), $(item.mjd), $(item.expnum), $(item.adjfiberindx )")
        end
    end
end
println("Batch information written to: $batch_info_file")
flush(stdout)

# M2: batch-level failures (prior load, HDF5 write, ...) must not kill the whole
# pmap; per-spectrum failures are already handled inside pipeline_single_spectra.
# Failed batches are logged and marked with 1 in pout_arMADGICS.txt (success = 0).
batch_on_error = ex -> begin
    println("multi_spectra_batch failed (batch skipped, marked 1 in pout): ", sprint(showerror, ex))
    flush(stdout)
    1
end
pout = @showprogress pmap(multi_spectra_batch, ittot; on_error = batch_on_error)
writedlm(out_dir * "pout_arMADGICS.txt", pout)
rmprocs(workers())

t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t0));
println("Total script runtime: $dt");
t_then = t_now;
flush(stdout);