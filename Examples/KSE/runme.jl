using Dates
using JLD
include("Model_KSE.jl")

## Parameters for KSE model
T = 10^4 # Length (in seconds) of time of run
T_disc = Int(T/2) # Length (in seconds) of time discarded
P = 21.55  # Period
N = 96  # Number of fourier modes used
h = 1e-3 # Timestep
g = x -> randn()*cos(2π*x/P)*(randn() + sin.(2π*x/P))
obs_gap = 100

uu, vv, tt =  my_KSE_solver(T,
    T_disc  = T_disc,
    P = P,
    N = N,
    h = h,
    g = g,
    n_gap = obs_gap)

# set save destinations
sol_file = "../data/KSE_Data/KSE_sol_Lin.jld"

paramaters = Dict(
      "T" => T,
      "T_disc" => T_disc,
      "P" => P,
      "N" => N,
      "h" => h,
      "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
      "obs_gap" => obs_gap,
      "tm" => now())
dat = Dict(
     "dat_uu" => uu,
     "dat_vv" => vv,
     "dat_tt" => tt)
Data = merge(paramaters,dat)
save(sol_file,Data)
