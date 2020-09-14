using Statistics
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
t_stop = 10^3
steps = 10^5 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [.35]
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
h = (t_stop - t_start)/(steps)

tim = range(t_start,t_stop,length = steps)

## Get Wiener Filter Ensemble

h_wf_ens = Run_and_get_WF(
    Psi;
    simulator = simulator,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx,
    Nen = Nen,
    M_out = M_out
    )
# Times
# RK4 10  = 1:45:50 - 1:47:05
# RK4 100 = 1:49:52 - 1:50:54
# FE 100  = 1:57:22 - 1:57:53

# Analyze ensamble to get the mean filter
A2 = analyse_h_ens(h_wf_ens)
A1[3]

T = h_wf_ens[1,1,:,:]'
T[:,1] .-= 1

S = h_wf_ens[1,2,:,:]'

dfT, dfS = DataFrame(T), DataFrame(S)
p1 = @df dfT boxplot(T,
    marker=(0.3,:orange,stroke(.5)),
    alpha=0.75,
    leg = :none,
    title = "First coefficients")
H = h_wf[1,1,:]; H[1] -= 1
scatter!(1:20,H,marker = (:star5,7))


p2 = @df dfS boxplot(S,
    marker=(0.3,:orange,stroke(.5)),
    alpha=0.75,
    leg = :none,
    title = "Second coefficients")
H = h_wf[1,2,:];
scatter!(1:20,H,marker = (:star5,7))

P = plot(p1,p2, layout = (2,1))
gui(P)

h_m = A2[1]
A2[3]
hva = A2[2]

signal = simulator(steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx,
    SM1 = true,
    Obs_noise = true)

obs_noise = signal[2]
signal = signal[1]
sig = signal[:,2:end]


# Notice since SM1 is true the output series will include one
# presample. Thus effecting an offset by one in the index.
# Signal is one time step behind, the desired sig.

timeseries_plot(tim,sig[1,:], title = "Original Signal")

gui()

pred = zeros(nu, steps)
for n = 1:steps
    pred[:,n] = Psi(signal[:,n])
end

h_wf = real.(vector_wiener_filter_fft(pred, sig, M_out))

h_ideal = zeros(1,2,20); h_ideal[1,:,1] = [1+h -h]

sig_hat_ideal = zeros(d,steps); sig_hat_ideal[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_ideal[:,i] = sum([real.(h_ideal)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

sig_hat_f = zeros(d,steps); sig_hat_f[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_f[:,i] = sum([real.(h_wf)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

sig_hat_m = zeros(d,steps); sig_hat_m[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_m[:,i] = sum([real.(h_m)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

win = 10080:10100

plot([sig[1,win .- 1] sig_hat_f[1,win] sig_hat_m[1,win] sig_hat_ideal[1,win]],
    label = ["sig" "sig_hat_f" "sig_hat_m" "sig_hat_ideal"])


sig_rm_ideal = redmodrun(
    h_ideal, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

P = timeseries_plot(tim,sig_rm_ideal[1,:],
    title = "Idea Reproduced Signal")
gui()



sig_rm_m = redmodrun(
    h_m, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

timeseries_plot(tim,sig_rm_m[1,:], title = "Mean Filter 500 samples Reproduced Signal")
gui()

sig_rm_wf = redmodrun(
    h_wf, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    obs_noise = obs_noise
    )

timeseries_plot(tim,sig_rm_wf[1,:], title = "Run-specific Filter Reproduced Signal")
gui()
