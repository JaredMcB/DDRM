using FFTW
using PyPlot


at    = include("../../Tools/AnalysisToolbox.jl")
mykse = include("myKSE_solver.jl")
# mykse_NA = include("myKSE_solver_NA.jl")
# mykse_A = include("myKSE_solver_A.jl")

kse   =  include("Model_kse.jl")

## Kassam and Trefethen paramaters
T       = 150
P       = 32π
n       = 64
h       = 1/4
g       = x -> cos(x/16)*(1 + sin(x/16))
T_disc  = 0
n_gap   = 6
aliasing= true

# My New Solver
vv = mykse.my_KSE_solver(T;
                        P,
                        n,
                        h,
                        g,
                        T_disc,
                        n_gap)

uu = real((2n+1)*ifft(vv,1))
ender = findfirst(isnan, sum(uu,dims = 1))[2]-10
figure()
H1 = imshow(reverse(reverse(uu[:,1:end]',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")
title("Aliasing = $aliasing")
EE = abs2.(vv)
figure()
semilogy(mean(EE[:,1:ender],dims=2),label = "me (Julia) $i")

lags = -50:50
A = at.my_autocor(vv[3,1:ender],lags)
plot(A)
# MyOld Solver
uu_o, vv_o, tt = kse.my_KSE_solver(T;
                        P,
                        N = 2n+1,
                        h,
                        g,
                        T_disc,
                        n_gap,
                        aliasing)

ender = findfirst(isnan, sum(uu_o,dims = 1))[2]-10
figure()
H1 = imshow(reverse(uu_o[:,1:ender]',dims = 1), extent=[0,P,0,150], aspect="auto")

EE_o = abs2.(vv_o)
semilogy(mean(EE_o[:,1:ender],dims=2),label = "me (Julia) $i")

A_o = at.my_autocor(vv_o[3,1:ender],lags)
plot(A_o)




## 2017 parameters
T        = 10^3 # Length (in seconds) of time of run
T_disc   = 0    # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
n        = 96  # Number of fourier modes used
h        = .001  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
n_gap    = 100
aliasing = false

# My New Solver
vv = mykse.my_KSE_solver(T;
                        P,
                        n,
                        h,
                        g,
                        T_disc,
                        n_gap,
                        aliasing)

uu = real((2n+1)*ifft(vv,1))
ender = findfirst(isnan, sum(uu,dims = 1))[2]-1000
H1 = imshow(reverse(uu[:,1:ender]',dims = 1), extent=[0,P,0,1000], aspect="auto")

EE = abs2.(vv)
figure()
semilogy(mean(EE[:,1:ender],dims=2),label = "me (Julia) $i")

lags = -50:50
A = at.my_autocor(vv[3,1:ender],lags)
plot(A)
# MyOld Solver
uu_o, vv_o, tt = kse.my_KSE_solver(T;
                        P,
                        N = 2n+1,
                        h,
                        g,
                        T_disc,
                        n_gap,
                        aliasing)

ender = findfirst(isnan, sum(uu_o,dims = 1))[2]-10
figure()
H1 = imshow(reverse(uu_o[:,1:ender]',dims = 1), extent=[0,P,0,150], aspect="auto")

EE_o = abs2.(vv_o)
semilogy(mean(EE_o[:,1:ender],dims=2),label = "me (Julia) $i")

A_o = at.my_autocor(vv_o[3,1:ender],lags)
plot(A_o)
