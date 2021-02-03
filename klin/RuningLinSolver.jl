################################################################################
#
# File: RunningLinSolver.jl
# Author: Jared McBride (Feb-2-2021)
#
# This just generates and saves KSE model runs using
# Dr. Lin's solver.
#
################################################################################

using JLD
using Dates

ks = include("ks.jl")

gen = "lin"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

T        = 10^3 # Length (in seconds) of time of run
T_disc   = T ÷ 2 # Length (in seconds) of time discarded (taken from T)
P        = 2π/sqrt(0.085)  # Period
n        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100


steps   = ceil(Int,T/h)
discard = ceil(Int,T_disc/h)

@time vv = ks.run_ks(ks.fei_init(),steps, .1;
             nsubsteps = obs_gap,
             verbose = false, L = P)

paramaters = Dict(
   "gen" => gen,
    "T" => T,
   "T_disc" => T_disc,
   "P" => P,
   "n" => n,
   "h" => h,
   "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
   "obs_gap" => obs_gap,
   "seed" => seed,
   "tm" => now()
   )

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")
sol_file = server ? "../../../data/KSE_Data/ks_sol_$gen.jld" :
   "Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"
println("Sol save location: " * sol_file)

dat = Dict("dat_vv" => vv)
Data = merge(paramaters, dat)
save(sol_file, Data)
println("data saved")
