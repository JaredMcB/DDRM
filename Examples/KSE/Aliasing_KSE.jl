using PyPlot
using Statistics: mean, var
using FFTW


at   = include("../../Tools/AnalysisToolbox.jl")
kse  = include("Model_KSE.jl")
ksed = include("Model_KSE_Dev.jl")

# 2017 parameters
T        = 10000 # Length (in seconds) of time of run
T_disc   = 0    # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = .001  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
obs_gap  = 60 #floor(Int, T/h/100)

Δt = h*obs_gap

uu_a, vv_a, tt   =  @time kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap);

# Trefethen parameters. 
uu, vv, tt = kse.my_KSE_solver(1500;T_disc = 0, n_gap = 60)
uu
uu_an = uu_a



H1 = imshow(reverse(reverse(uu,dims = 2),dims = 1)', extent=[0,P,0,150], aspect="auto")

uu_a2n, vv_a, tt   =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap);

figure()
H2 = imshow(reverse(reverse(uu_a2n[:,:],dims = 2),dims = 1)', extent=[0,P,0,150], aspect="auto")

uu_a2n

uu_a2n_org, vv_a, tt   =  @time kse.my_KSE_solver(T; T_disc, P, N = 2N+1, h, g, n_gap = obs_gap);
figure()
H1 = imshow(reverse(reverse(uu_a2n_org[:,:],dims = 2),dims = 1)', extent=[0,P,0,150], aspect="auto")
