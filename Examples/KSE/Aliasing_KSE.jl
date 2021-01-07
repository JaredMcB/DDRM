using PyPlot
using Statistics: mean, var
using FFTW


at    = include("../../Tools/AnalysisToolbox.jl")
kse   = include("Model_KSE.jl")
ksed2 = include("Model_KSE_Dev2.jl")

## 2017 parameters
T        = 10^3 # Length (in seconds) of time of run
T_disc   = 0    # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = .001  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
obs_gap  = 100 #floor(Int, T/h/100)

Δt = h*obs_gap

uu_a, vv_a, tt   = @time kse.my_KSE_solver(T; P, n = N, h, g,
    n_gap = obs_gap,
    aliasing = false);

uu_o, vv_o, tt   = @time ksed2.my_KSE_solver(T; P, N, h, g, n_gap = obs_gap);




## Trefethen parameters.
uu, vv, tt = kse.my_KSE_solver(1500;
                               N = 128,
                               T_disc = 0,
                               n_gap = 6,
                               aliasing = false)
uu

plot(32π*(0:127)/128,uu[:,end], label = "me (Julia)")
title("KSE solution at T = 150 (both with aliasing)")

UU = uu
VV = vv


ender = findfirst(isnan, sum(UU,dims = 1))[2]-5
EE = abs2.(VV)

N = 128
i += 1
semilogy(mean(EE[:,1:ender],dims=2)*2N^2,label = "me (Julia) $i")
title("Energy spectrum with aliasing")
legend()
for i in 1:20:200
    semilogy(EE[:,i],label = i,lw = 1)
    legend()
end



lags = -150:150
C2 = at.my_autocor(vv[3,:],lags)

plot(lags,[C2 c2_ml],".-")
grid()
H1 = imshow(reverse(reverse(uu,dims = 2),dims = 1)', extent=[0,P,0,150], aspect="auto")

uu_a2n, vv_a, tt   =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap);

figure()
H2 = imshow(reverse(reverse(uu_a2n[:,:],dims = 2),dims = 1)', extent=[0,P,0,150], aspect="auto")

uu_a2n

uu_a2n_org, vv_a, tt   =  @time kse.my_KSE_solver(T; T_disc, P, N = 2N+1, h, g, n_gap = obs_gap);
figure()
H1 = imshow(reverse(reverse(uu_a2n_org[:,:],dims = 2),dims = 1)', extent=[0,P,0,150], aspect="auto")
