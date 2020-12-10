using PyPlot
using Statistics: mean, var

kse = include("../Model_KSE.jl")
ksed = include("../Model_KSE_Dev.jl")

T        = 150 # Length (in seconds) of time of run
T_disc   = 0    # Length (in seconds) of time discarded
P        = 32π  # Period
N        = 128  # Number of fourier modes used
h        = 1/4  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
obs_gap  = floor(Int, T/h/100)

Δt = h*obs_gap
uu, vv, tt    =  @time kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)
uud, vvd, ttd =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)

t_start = 0
t_stop = 150
ind_start = floor(Int,t_start/Δt)+1
ind_stop =floor(Int,t_stop/Δt)
H1 = imshow(reverse(uu,dims = 2)', extent=[0,P,0,150], aspect="auto")

mean(uu, dims = 2)

u = rand(128)
using FFTW

u - fft(ifft(u))

sum(abs.(uu - uud).^2)
uu
