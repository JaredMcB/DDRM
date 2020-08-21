"""
Today we will investigate:

Lorenz 63

"""

using Statistics
using StatsBase
using Plots
using StatsPlots
using DataFrames
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
include("..\\DataGen.jl")
include("..\\RedModRun.jl")

## Preference parameters
t_start = 0
t_stop = 10^4
steps = 10^6 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

d = 3
sig_init = randn(3,1)
sigma = I + zeros(d,d)
sigma_v = sigma

Nen = 10

Psi(x) = [x;
          x[1]*x[3];
          x[1]*x[2] ]

M_out = 10

# Derivitive parameters
steps_tot = steps + discard
nu = size(Psi(zeros(3,1)),1)
dt = (t_stop - t_start)/(steps)
tim = range(t_start,t_stop,length = steps)

signal = DataGen_Lorenz_63(
    steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    SM1 = true,
    Obs_noise = false
    )

wind = 1000:2000
plot(tim[wind],signal[:, wind]')

plot(sig[1,wind],sig[2,wind],sig[3,wind],
    seriestype = :scatter)

sig = signal[:,2:end]

h_wf= get_wf(signal,Psi)

h_wf_ens = Run_and_get_WF(Psi,
    simulator = DataGen_Lorenz_63,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = 3,
    V_prime = (x -> -x.*(x.^2 .- 1)),
    Nen = 10,
    M_out = 10,
    rl = true # if ture output made REAL
    )

h_mean = analyse_h_ens(h_wf_ens,plt = false)

h_var = h_mean[2]
h_m = h_mean[1]

sig_rm = redmodrun(
    h_wf, Psi;
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard
    )

plot(sig_rm[1,wind],sig_rm[2,wind],sig_rm[3,wind],
        seriestype = :scatter)

### Compare ACF's
L = 1000; lags = 0:L
A = zeros(L+1,2*d)
A[:,1] = autocov(sig[1,:],lags)
A[:,2] = autocov(sig[2,:],lags)
A[:,3] = autocov(sig[3,:],lags)
A[:,4] = autocov(sig_rm[1,:],lags)
A[:,5] = autocov(sig_rm[2,:],lags)
A[:,6] = autocov(sig_rm[3,:],lags)

plot(lags, A,
    line = (2,[:solid :solid :solid :dot :dot :dot]),
    color = [:red :green :blue])
