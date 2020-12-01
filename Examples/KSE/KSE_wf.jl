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

Exp = "12_01_20_1"

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
nfft = 2^14
par = 1500
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
sol_file = server ? "../../../data/KSE_Data/KSE_sol_lin.jld" :
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
function InvBurgRK4_1step(x)
   lx = length(x)
   function F(x)
      𝑥 = [conj(reverse(x, dims = 1)); 0; x]
      -im / 2 * (2π / P * (1:lx) / N) .* conv(𝑥, 𝑥)[2*lx+2:3*lx+1]
   end

   Δt = h * obs_gap

   k1 = F(x)
   k2 = F(x .+ Δt * k1 / 2)
   k3 = F(x .+ Δt * k2 / 2)
   k4 = F(x .+ Δt * k3)
   A = @. x + Δt / 6 * (k1 + 2k2 + 2k3 + k4)
end

function Inertialman_part(x)
   lx = length(x)
   𝑥(j) = (j <= lx ? x[j] : im * sum(x[l] * x[j-l] for l = j-lx:lx))

   L = complex(zeros(lx^2))
   for j = 1:lx
      for k = 1:lx
         L[(j-1)*lx+k] = 𝑥(j + lx) * 𝑥(j + lx - k)
      end
   end
   L
end

function Inertialman_part_short(x)
   lx = length(x)
   𝑥(j) = (j <= lx ? x[j] : im * sum(x[l] * x[j-l] for l = j-lx:lx))

   L = complex(zeros(binomial(lx + 1, 2)))
   i = 1
   for j = 1:lx
      for k = j:lx # k should normaly go from 1 to lx but i changed it to go from j to lx.
         L[i] = 𝑥(j + lx) * 𝑥(j + lx - k)
         i += 1
      end
   end
   L
end

Psi(x) = short ? [x; InvBurgRK4_1step(x); Inertialman_part_short(x)] :
   [x; InvBurgRK4_1step(x); Inertialman_part(x)]


## Get Wiener filter

print("Get_wf computation time: ")
@time h_wf = mr.get_wf(X, Psi; M_out, par, nfft, rl = false, PI = false)

# Save Wienerfilter
dat = Dict("dat_h_wf" => h_wf)
Data = merge(paramaters,dat)
save(wf_file, Data)
println("Wiener filter saved")
