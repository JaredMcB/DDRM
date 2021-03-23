
using JLD
using DSP # For conv function in Psi
using Dates

mrb = include("../../Tools/WFMR_bs.jl")
at = include("../../Tools/AnalysisToolbox.jl")

# Load Old Data

for i = 1:2
    gen = "lin1e5_r$i"     # this is just a reference designation it shows up in the
                    # output file. I think of generatrion.

    server = startswith(pwd(), "/u5/jaredm") ? true : false
    println("on server = $server")

    sol_file = server ? "../../../data/KSE_Data/ks_sol_$gen.jld" :
       "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"
    println("Sol load location: " * sol_file)

    @time vv = load(sol_file,"dat_vv")


    ## Get Reduced model #########################################################
    # Model reductrion parameters

    d = 5
    h = 0.1
    # collect observations
    obs_gap = 1
    V_obs = vv[2:d+1,1:obs_gap:end]
    vv = []

    # Build PSI
    Psi = mrb.get_Psi_2017(h)

    # Get Wiener filter
    #@time h_wf = get_wf(V_obs,Psi, M_out = M_out,PI = true)
    signal = V_obs
    V_obs = []
    M_out = 20
    n = 3; p = 1500; par = 1500
    ty = "bin"
    xspec_est = "old"
    nfft = 0
    rl = true
    Preds = false
    N_ckms = 3000
    PI = false
    rtol = 1e-6
    Verb = false
    tm = now()

    for M_out in [20, 60, 80, 100, 120]

       paramaters = Dict(
           "M_out" => M_out,
           "n" => n,
           "p" => p,
           "par" => par,
           "ty" => ty,
           "xspec_est" => xspec_est,
           "nfft" => nfft,
           "rl" => rl,
           "Preds" => Preds,
           "N_ckms" => N_ckms,
           "rtol" => rtol,
           "tm" => tm
       )


       h_wf_bs = @time mrb.get_wf_bs(signal, Psi; M_out)

       # Save Wienerfilter
       dat = Dict("dat_h_wf" => h_wf_bs)
       Data = merge(paramaters,dat)
       # save("../data/KSE_Data/KSE_sol_wienerfilter.jld",Data)

       wf_file = server ? "../../../data/KSE_Data/ks_wf_bs$M_out-$gen.jld" :
          "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_wf_bs$M_out-$gen.jld"
       save(wf_file,Data)
       end

       println("Wiener filters $i saved")
end

# # Load old Wiener Filter
# LD = load("Data\\KSE_sol_wienerfilter.jld")
# h_wf = LD["dat_h_wf"]
# println("Wiener filter load
