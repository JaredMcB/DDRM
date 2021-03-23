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

kse = include("myKSE_solver.jl")

gen = "lin1e5"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

T        = 10^5 # Length (in seconds) of time of run
T_disc   = T ÷ 2 # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
n        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
obs_gap  = 100

for i = 1:2

   gen = "lin1e5_r$i"
   r = randn(6)
   g(x) = sum(r[n+1]*cos(2π/P*n*x) for n = 0:5)

   vv = @timev kse.my_KSE_solver(T; T_disc, P, n, h, g, n_gap = obs_gap)

   paramaters = Dict(
      "gen" => gen,
       "T" => T,
      "T_disc" => T_disc,
      "P" => P,
      "n" => n,
      "h" => h,
      "g" => "x -> sum(r[n+1]*cos(2π/P*n*x) for n = 0:5)",
      "r" => r,
      "obs_gap" => obs_gap,
      "tm" => now()
      )

   server = startswith(pwd(), "/u5/jaredm") ? true : false
   println("on server = $server")
   sol_file = server ? "../../../data/KSE_Data/ks_sol_$gen.jld" :
      "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"

   println("Sol save location: " * sol_file)

   dat = Dict("dat_vv" => vv)
   Data = merge(paramaters, dat)
   save(sol_file, Data)
   println("data $i saved")

end
