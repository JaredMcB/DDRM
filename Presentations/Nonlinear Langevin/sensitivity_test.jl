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
sigma = [1]
sigma_v = sigma

d = 1

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)
Psi(x) = [x; x.^3]

M_out = 20
##

nu = size(Psi(0),1)
h = (t_stop - t_start)/(steps)
Nen = 100 # The number of series in the ensamble
h_wf_ens = complex(zeros(d,nu,M_out,Nen))

for k = 1:Nen
    signal = DataGen_Langevin(10^5 + 1,
        t_start = 0,
        t_stop = 10^3,
        discard = 10^4,
        sig_init = [1.5],
        sigma = [1],
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
    M_out = 20
    h_wf = vector_wiener_filter_fft(pred, sig, M_out)
    h_wf_ens[:,:,:,k] = h_wf
end

h_wf_mean = mean(h_wf_ens,dims = 4)
h_wf_var = var(h_wf_ens,dims = 4)

plot(real.(h_wf_ens[1,2,2:20,:]),
    leg = false,
    line = 1)

h_wf = h_wf_mean[:,:,:,1]

sig_rm = redmodrun(
    h_wf, Psi,
    sigma_v = [1],
    sig_init = [1.5],
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard
    )

plot(sig_rm)


signal = DataGen_Langevin(10^5 + 1,
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [1],
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
