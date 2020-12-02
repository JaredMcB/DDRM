################################################################################
#
# File: KSE_wf.jl
# Author: Jared McBride (11-23-2020)
#
# This is an experiment script originally for Exper_11_23_2020_1
#
# This file is simply designed to run KSE data (or load it) and us the gerenal
# Wiener filtering routine ("Tools/Model_Reduction_Dev.jl") to produce the wiener
# filter.
#
# It forms the psi function based on that of Lin, Lu, and Chorin (2017).
#
################################################################################

using JLD
using DSP: conv # For conv function in Psi
using Dates

mr  = include("../../Tools/Model_Reduction_Dev.jl")
kse = include("Model_KSE.jl")
kmr = include("KSE_modredTools.jl")

Exp = "12_02_20_4"

## Parameters
# run parameters
T = 10^5              # Length (in seconds) of time of run
T_disc = T ÷ 2               # Length (in seconds) of time discarded
P = 2π / sqrt(0.085)    # Period
N = 96                # Number of fourier modes used
h = 1e-3              # Timestep
g = x -> cos(x) * (1 + sin.(x))

q = 2π / P * (0:N-1)      # Wve numbers (derivative)

# Observation parameters
obs_gap = 100
d = 5 # No. of lowest modes taken in reduced model

# Wiener filtering parameters
M_out = 20
nfft = 2^10
par = 100
xspec_est = "old" # Default
short = true
loadsol = true

paramaters = Dict(
   "Exp" => Exp,
   "T" => T,
   "T_disc" => T_disc,
   "P" => P,
   "N" => N,
   "h" => h,
   "g" => "x -> cos(x) * (1 + sin.(x))",
   "q" => q,
   "d" => d,
   "tm" => now(),
   "M_out" => M_out,
   "nfft" => nfft,
   "par" => par,
   "obs_gap" => obs_gap,
   "d" => d, # No. of lowest modes taken in reduced model
   "M_out" => M_out,
   "short" => short,
   "loadsol" => loadsol,
)
println("the Parameters ===================")
for x in keys(paramaters)
   println(x * " : ", paramaters[x])
end
println("==================================")

## Get full run
# decide on save location based of platform
server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")
sol_file = server ? "../../../data/KSE_Data/KSE_sol$Exp.jld" :
   "Examples/KSE/Data/KSE_sol$Exp.jld"
println("Sol save location: " * sol_file)
wf_file = server ? "../../../data/KSE_Data/KSE_wf$Exp-Mo$M_out.jld" :
   "Examples/KSE/Data/KSE_wf$Exp-Mo$M_out.jld"

# When I want the standard lin et al. (2017) data.
sol_file = server ? "../../../data/KSE_Data/KSE_sol_lin1.jld" :
   "Examples/KSE/Data/KSE_sol_lin.jld"

if loadsol
   # Load Old Data
   @time Data = load(sol_file)
   print("Data Loaded")
   uu = Data["dat_uu"]
   vv = Data["dat_vv"]
   tt = Data["dat_tt"]
else
   uu, vv, tt = kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)

   dat = Dict("dat_uu" => uu, "dat_vv" => vv, "dat_tt" => tt)
   Data = merge(paramaters, dat)
   save(sol_file, Data)
   println("data saved")
end

## Get Observations
X = vv[2:d+1, 1:end]

## Build Psi
Psi(x) = kmr.PSI(x)
## Get Wiener filter

print("Get_wf computation time: ")
@time h_wf = mr.get_wf(X, Psi; M_out, par, nfft, rl = false, PI = false)

# Save Wienerfilter
dat = Dict("dat_h_wf" => h_wf)
Data = merge(paramaters,dat)
save(wf_file, Data)
println("Wiener filter saved")
