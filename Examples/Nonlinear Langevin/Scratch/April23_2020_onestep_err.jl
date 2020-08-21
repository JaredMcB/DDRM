"""
Today we will investigate:
1. The effect zeroing out very small coefficients of the wiener filter.

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
steps = 10^5 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [.4]
sigma_v = sigma
d = 1

Nen = 500

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)
Psi(x) = [x; x.^3]

simulator = DataGen_Langevin_FE

M_out = 20

# Derivitive parameters
steps_tot = steps + discard
nu = size(Psi(0),1)
dt = (t_stop - t_start)/(steps)
tim = range(t_start,t_stop,length = steps)

## Code begins

signal = simulator(steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx,
    SM1 = false,
    Obs_noise = true)

obs_noise = signal[2]
signal = signal[1]

sig = signal[:,2:end]
pred = get_pred(sig,Psi)

h_wf= get_wf(signal,Psi)

h_ideal = zeros(1,2,20); h_ideal[1,:,1] = [1+dt -dt]
h_trunc_1 = h_wf[:,:,1:1]; M_trunc_1 = 1
h_trunc_5 = h_wf[:,:,1:5]; M_trunc_5 = 5


## prep the one step predictors
sig_hat_ideal = zeros(d,steps); sig_hat_ideal[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_ideal[:,i] = sum([real.(h_ideal)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

sig_hat_f = zeros(d,steps); sig_hat_f[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_f[:,i] = sum([real.(h_wf)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

sig_hat_trunc_1 = zeros(d,steps); sig_hat_trunc_1[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_trunc_1[:,i] = sum([real.(h_trunc_1)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_trunc_1)])
end

# Plot estimates
win = 10080:10100
plot([sig[1,win .- 1] sig_hat_f[1,win] sig_hat_trunc_1[1,win] sig_hat_ideal[1,win]],
    label = ["sig" "sig_hat_f" "sig_hat_trunc_1" "sig_hat_ideal"])
