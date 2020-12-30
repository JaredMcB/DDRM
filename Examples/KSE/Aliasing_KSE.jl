using PyPlot
using Statistics: mean, var
using FFTW


at   = include("../../Tools/AnalysisToolbox.jl")
kse  = include("Model_KSE.jl")
ksed2 = include("Model_KSE_Dev2.jl")

# 2017 parameters
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




# Trefethen parameters.
uu, vv, tt = kse.my_KSE_solver(150;n = 64,T_disc = 0, n_gap = 6,aliasing = false)
uu

UU = uu
VV = vv

findfirst(isnan, sum(UU,dims = 1))
EE = abs2.(VV)

semilogy(mean(EE,dims=2)*2N^2,label = "aliasing
controled")
legend()
for i in 1:10:100
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

EEE = [1.81e-29, 549.7216426, 399.4952178,
    132.8961573,
    145.9358596,
    186.0774939,
    160.3824272,
    180.2053082,
    449.8552998,
    314.8072577,
    2390.572761,
    778.0807514,
    887.2342493,
    472.4222544,
    492.439317,
    828.1559251,
    412.6699241,
    275.8383614,
    515.1858937,
    125.7039713,
    129.5540942,
    118.3811256,
    69.82428727,
    33.25648266,
    43.19265574,
    27.49228091,
    23.89401988,
    12.06031186,
    16.03049699,
    6.838577858,
    6.46509257,
    2.705236466,
    3.022531191,
    1.323293487,
    0.77401524,
    0.608287976,
    0.484772079,
    0.2226381,
    0.294299149,
    0.147621629,
    0.148894384,
    0.069739826,
    0.062648917,
    0.02702955,
    0.022809037,
    0.012599292,
    0.008948572,
    0.005721031,
    0.004197851,
    0.00215609,
    0.001782571,
    0.001140648,
    0.001004431,
    0.000540758,
    0.000457465,
    0.000214173,
    0.000149708,
    8.59e-05,
    7.70e-05,
    3.12e-05,
    3.66e-05,
    1.26e-05,
    2.05e-05,
    8.47e-06,
    8.03e-30,
    8.47e-06,
    2.05e-05,
    1.26e-05,
    3.66e-05,
    3.12e-05,
    7.70e-05,
    8.59e-05,
    0.000149708,
    0.000214173,
    0.000457465,
    0.000540758,
    0.001004431,
    0.001140648,
    0.001782571,
    0.00215609,
    0.004197851,
    0.005721031,
    0.008948572,
    0.012599292,
    0.022809037,
    0.02702955,
    0.062648917,
    0.069739826,
    0.148894384,
    0.147621629,
    0.294299149,
    0.2226381,
    0.484772079,
    0.608287976,
    0.77401524,
    1.323293487,
    3.022531191,
    2.705236466,
    6.46509257,
    6.838577858,
    16.03049699,
    12.06031186,
    23.89401988,
    27.49228091,
    43.19265574,
    33.25648266,
    69.82428727,
    118.3811256,
    129.5540942,
    125.7039713,
    515.1858937,
    275.8383614,
    412.6699241,
    828.1559251,
    492.439317,
    472.4222544,
    887.2342493,
    778.0807514,
    2390.572761,
    314.8072577,
    449.8552998,
    180.2053082,
    160.3824272,
    186.0774939,
    145.9358596,
    132.8961573,
    399.4952178,
    549.7216426]

semilogy(EEE, label = "Matlab")
legend()
