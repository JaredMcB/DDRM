# ##########################################################################
#  Title : sensitiviy_test.jl
#  Author: Jared McBride (04/09/2020)
#
#  Description: Here we investigate the sensitivty the Wiener filter may have
#  to a particular realization. We generate
#     N realizations of model
#     steps long.
#  for each realization a Wiener filter is computed of length
#     h_m we then analyse the ensamble of wiener filters
#
# ##########################################################################

# Packages used
using Statistics
using Plots
pyplot()

include("..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
include("DataGen.jl")
include("RedModRun.jl")
include("my_Histogram.jl")

## All the parameters
t_start = 0
t_stop = 10^3
steps = 10^5 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [0.2]
sigma_v = sigma

d = 1

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)

M_out = 20

## For investigating PSI

Psi1(x) = [x; x.^3 .- x]
Psi2(x) = [x; x.^3]
PSI = [Psi1, Psi2]

## For investigating variance sensitivity
h_wf_results = []


for i = 1: length(PSI)

    h_wf_ens = Run_and_get_WF(
        PSI[i],
        simulator = DataGen_Langevin_FE,
        steps = steps,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard,
        sig_init = sig_init,
        sigma = sigma,
        d = d,
        V_prime = dVdx,
        Nen = 100,
        M_out = M_out
        )
    ##

    append!(h_wf_results, [h_wf_ens])
end

h_wf_results_means = []
for H in h_wf_results
    append!(h_wf_results_means,[mean(H,dims = 4)])
end


plot(real.([h_wf_results_means[1][1,1,1:20] h_wf_results_means[2][1,1,1:20]]),
    leg = false,
    line = 1,
    marker = (3,[:h :h]))

plot(real.([h_wf_results_means[1][1,2,1:20] h_wf_results_means[2][1,2,1:20]]),
    leg = false,
    line = 1,
    marker = (3,[:h :h]))

# print(real.(h_wf_results[1,1,1:5,:]))

# h_wf = h_wf_mean[:,:,:,1]

sig_rm = zeros(d,steps,2)
for i = 1:2
    h_wf = real.(h_wf_results_means[i][:,:,1:20])
    sig_rm[:,:,i] = redmodrun(
        h_wf, PSI[i],
        sigma_v = [1],
        sig_init = [1.5],
        steps = steps,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard
        )
end
plot([sig_rm[1,:,1] sig_rm[1,:,1]])


signal = DataGen_Langevin(10^5 + 1,
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [.2],
    d = 1,
    V_prime = dVdx,
    SM1 = true)
# Notice since SM! is true the output series will include one
# presample. Thus effecting an offset by one in the index.
# Signal is one time step behind, the desired sig.

sig = signal[:,2:end]

pred = zeros(nu, steps)
for n = 1:steps
    pred[:,n] = Psi(signal[n])
end

h_wf = vector_wiener_filter_fft(pred, sig, M_out)

sig_hat = zeros(1,steps)
for i = 1:steps
    sig_hat[:,i] = sum([real.(h_wf[:,:,k+1])*pred[:,i-k] for k = 0:min(i - 1,M_out - 1)])
end

view = 100
start = steps - 10 - view
tim = start:start + view;

sig_hat_mean = sig_hat

plot(tim,[sig[1,tim ] sig_hat[1,tim] sig_hat_mean[1,tim]],
    color=[:red :black :orange],
    line=(1,[:solid :dash :dot]),
    label=["sig" "sig_hat" "sig_hat_mean"])
