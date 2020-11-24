################################################################################
#
# File: KSE_data_gen.jl
# Author: Jared McBride (11-23-2020)
#
# This just generates and saves KSE model runs
#
################################################################################

using JLD
using DSP: conv # For conv function in Psi
using Dates
using Random

include("Model_KSE.jl")

T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2 # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100
seed     = 2020

Random.seed!(seed)
uu, vv, tt =  my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)

paramaters = Dict(
   "T" => T,
   "T_disc" => T_disc,
   "P" => P,
   "N" => N,
   "h" => h,
   "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
   "obs_gap" => obs_gap,
   "tm" => now()
   )

sol_file = "Examples/KSE/Data/KSE_sol_lin.jld"

dat = Dict("dat_uu" => uu, "dat_vv" => vv, "dat_tt" => tt)
Data = merge(paramaters, dat)
save(sol_file, Data)
println("data saved")
