
using JLD
using DSP # For conv function in Psi
using Dates

mrb = include("../../Tools/WFMR_bs.jl")
mr  = include("../../Tools/WFMR.jl")
at  = include("../../Tools/AnalysisToolbox.jl")

load_pred = true

# Load Old Data
geni = 1
gen = "lin1e5_r$geni"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

println("Getting preds")
save_file_pred = "../../../../data/KSE_Data/ks_pred_$gen.jld" 
if load_pred
    pred = load(save_file_pred,"pred")
    println("loaded pred from: " * save_file_pred)
else
    sol_file = "../../../../data/KSE_Data/ks_sol_$gen.jld" 
    println("Sol load location: " * sol_file)

    @time vv = load(sol_file,"dat_vv");

    # get observations
    d = 5
    h = 0.1
    # collect observations
    obs_gap = 1
    signal = vv[2:d+1,1:obs_gap:end]

    Psi = mr.get_Psi_2017(h)

    pred = @timev complex(mr.get_pred(signal, Psi)[:,1:end-1])
    sig = complex(signal[:,2:end]);
    save(save_file_pred,"pred",pred,"sig",sig)
    println("saved sig and pred in: " * save_file_pred)
end

# Get Whitening filter
println("Computing WFBS")
@timev h_wf_bs = mrb.get_wf_bs(sig, pred, M_out = 50)

savefile_wf_bs = "../../../../data/KSE_Data/ks_wf_bs_$gen.jld"
save(savefile_wf_bs,"h_whf",h_wf_bs)
println("Whitening filter saved at "*savefile_wf_bs) 