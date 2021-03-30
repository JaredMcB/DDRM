using JLD
using DSP
using FFTW
using PyPlot
using Statistics

mr = include("../../../Tools/WFMR.jl") # Now includes the whitening filter
at = include("../../../Tools/AnalysisToolbox.jl")

# Load Old Data
gen = "lin1e5_r4"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")

sol_file = server ? "../../../../data/KSE_Data/ks_sol_$gen.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"
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

# Load Whitening filter
h_whf = load("../../../../data/KSE_Data/ks_whf_$gen.jld","h_whf")

## truncate it
M_out = 10000
h = h_whf[:,:,1:M_out];

@timev y = at.my_filt(h, pred)

# Save the whitened process.
save("../../../../data/KSE_Data/ks_whsol_$gen.jld","wh_pred",y,"M_out",M_out)