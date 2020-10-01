"""
Today we will investigate:
1. The effect zeroing out very small coefficients of the wiener filter.

"""

using Statistics
using StatsBase
using Plots
using StatsPlots
using DataFrames
#ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


# include("..\\..\\..\\Tools\\Wiener Filtering\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
include("..\\DataGen.jl")
include("../Wiener_Filter_analyzer.jl")
include("../../../Tools/Model_Reduction_Dev.jl")
# include("..\\RedModRun.jl")

## Preference parameters
t_start = 0
t_stop = 10^4
steps = 10^5 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [.45]
sigma_v = sigma
d = 1

Nen = 100

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)
Psi(x) = [x; (x).^3]

# simulator = DataGen_Langevin_FE

M_out = 20

# Derivitive parameters
steps_tot = steps + discard
nu = size(Psi(0),1)
dt = (t_stop - t_start)/(steps)
tim = range(t_start,t_stop,length = steps)

## Code begins

Nen = 10

h_ens = Run_and_get_WF_DWOL(
    Psi;
    steps,
    t_start,
    t_stop,
    discard,
    sig_init,
    sigma,
    d = d,
    V_prime = dVdx,
    Nen = Nen,
    M_out = 20)

h_mean = analyse_h_ens(h_ens)

h_mean[3]
h_var = h_mean[2]
h_m = h_mean[1]

x= -3:.01:3

h_m[1,:,1]'*Psi(5)

plot(x,[x .- dt*(x.^3 .- x) map(x -> h_m[1,:,1]'*Psi(x),x)])

# signal = simulator(steps,
#     t_start = t_start,
#     t_stop = t_stop,
#     discard = discard,
#     sig_init = sig_init,
#     sigma = sigma,
#     d = d,
#     V_prime = dVdx,
#     SM1 = false,
#     Obs_noise = true)

SM1 = false
Obs_noise = true
d = 1
e = randn(d,steps + discard)

signal = DataGen_DWOL(
    steps;
    t_start, t_stop, discard,
    sig_init , sigma, V_prime = dVdx,
    SM1, Obs_noise, d, e
    )

obs_noise = signal[2]
sig = signal[1]

pred = get_pred(sig,Psi)

timeseries_plot(tim,sig[1,:],
    title = "Signal")

h_wf= get_wf(sig,Psi)

h_ideal = zeros(1,2,20); h_ideal[1,:,1] = [1+dt -dt]
h_trunc_1 = h_wf[:,:,1:1]; M_trunc_1 = 1
h_trunc_5 = h_wf[:,:,1:5]; M_trunc_5 = 5

## Repreduced Model Information

sig_rm_ideal = redmodrun(h_ideal, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

# P = timeseries_plot(tim,sig_rm_m[1,:],
#     title = "Signal")
# gui()

sig_rm_trunc_1 = redmodrun(h_trunc_1, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

sig_rm_trunc_5 = redmodrun(h_trunc_5, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

# timeseries_plot(tim,sig_rm_trunc_5[1,:],
#     title = "trunc_5 filter Reproduced Signal")
# gui()

sig_rm_wf = redmodrun(h_wf, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

# timeseries_plot(tim,sig_rm_wf[1,:],
#     title = "Run-specific Filter Reproduced Signal")
# gui()

wind = 1:100
plot(tim[wind], [sig[1,wind] sig_rm_trunc_1[1,
    wind] sig_rm_trunc_5[1,wind] sig_rm_wf[1,wind]],
    line = (2,[:solid :dash :dot :dot]),
    label = ["sig" "sig_rm_trunc_1" "sig_rm_trunc_5" "sig_rm_wf"],
    title = "Double-well Signal and sine model reduction")

lags = 0:2000
A_o = autocov(sig[1,:],lags)
A_wf = autocov(sig_rm_wf[1,:],lags)
A_trunc_1 = autocov(sig_rm_trunc_1[1,:],lags)
A_trunc_5 = autocov(sig_rm_trunc_5[1,:],lags)

plot([A_o A_wf A_trunc_1 A_trunc_5],
    xlabel ="lags",
    labels = ["origninal" "20 filter" "trunc_1" "trunc_5"],
    line =([:solid :dash :dot :dot],3),
    title = "ACF Double-well Signal and sine model reduction")
gui()
