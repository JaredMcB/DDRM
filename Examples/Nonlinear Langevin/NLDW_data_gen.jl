################################################################################
#
# File: NLDW_data_gen.jl
# Author: Jared McBride (11-23-2020)
#
# This just generates and saves KSE model runs
#
################################################################################

using JLD
using DSP: conv # For conv function in Psi
using Dates

dg = include("DataGenDWOL.jl") # This has many packages in it's preamble

gen = "lin1e6"     # this is just a reference designation it shows up in the
                   # output file. I think of generatrion.
steps = 10^5 + 1
scheme = "FE"
h = .1
discard = 10^6
sig_init = [1.5]
sigma = [0.5]
V_prime = x -> -x.*(x.^2 .- 1)
ObsNoise = false
d = 1
gap = 1
#e = randn(d,steps + discard)

# Get full model run
Random.seed!(2014)
X = dg.DataGen_DWOL(;
  steps,
  scheme, h, discard,
  sig_init, sigma, V_prime, ObsNoise, gap
  )


paramaters = Dict(
   "gen" => gen,
    "steps" => steps,
   "scheme" => scheme,
   "h" => h,
   "discard" => discard,
   "sig_init" => sig_init,
   "V_prime" => "x -> -x.*(x.^2 .- 1)",
   "gap" => gap,
   "d" => d,
   "tm" => now()
   )

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")
sol_file = server ? "../../../data/NDWL_Data/ks_sol_$gen.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"

println("Sol save location: " * sol_file)

dat = Dict("dat_vv" => vv)
Data = merge(paramaters, dat)
save(sol_file, Data)
println("data saved")
