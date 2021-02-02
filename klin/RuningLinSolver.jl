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

T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2 # Length (in seconds) of time discarded (taken from T)
P        = 2π/sqrt(0.085)  # Period
n        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100


steps   = ceil(Int,T/h)
discard = ceil(Int,T_disc/h)

vv = ks.run_ks(steps, .1;
             nsubsteps = obs_gap,
             verbose = false, L = P)

paramaters = Dict(
   "gen" => gen,
    "T" => T,
   "T_disc" => T_disc,
   "P" => P,
   "N" => N,
   "h" => h,
   "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
   "obs_gap" => obs_gap,
   "seed" => seed,
   "tm" => now()
   )

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")
sol_file = server ? "../../../data/KSE_Data/KSE_sol_$gen.jld" :
   "Desktop/DDMR/Examples/KSE/Data/KSE_sol_$gen.jld"
println("Sol save location: " * sol_file)

dat = Dict("dat_uu" => uu, "dat_vv" => vv, "dat_tt" => tt)
Data = merge(paramaters, dat)
save(sol_file, Data)
println("data saved")




using PyPlot
using FFTW

ks = include("ks.jl")
at = include("../Tools/AnalysisToolbox.jl")


P = 32π
steps = 100
N = 128

# Trefethen
vv = ks.run_ks(steps,1.5;
             nsubsteps = 6,
             verbose = false, L = P)

# get solution
uu = real(N*ifft([zeros(1, steps+1); vv; reverse(conj(vv),dims = 1)],1))
H1 = imshow(reverse(reverse(uu',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")


VV = [vv_lin, vv_myNA, vv_myA]
UU = [uu_my,uu_lin, uu]

sty = ["-","--",":"]
labels = ["vv_lin", "vv_myNA", "vv_myA"]
for i = 1:5
    for j = 1 : length(VV)
        subplot(5,1,i)
        plot(1.5*(1:300),real(VV[j][i,1:300]),linestyle = sty[j])
    end
    title("model $i")
end
suptitle("Comparison of L's nonlinear part (solid) with M's nonlinear part (dealiased, dashed;aliased, dotted) in L's ODE solver")






v = ks.kt_init()

vc = real(N*ifft([0; v[:];0;reverse(conj(v[:]))]))

vi = ks.kt_init(real = true)

maximum(abs.(vi - vc))

EE = abs2.(vv)

semilogy(mean(EE,dims = 2))

# 2017
steps = 10000
vv = ks.run_ks(ks.fei_init(), steps,.1;
             nsubsteps = 100,
             verbose = false)

A = at.my_autocor(vv[2,:],-150:150)

plot(A)
