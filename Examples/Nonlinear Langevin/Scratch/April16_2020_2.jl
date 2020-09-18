using Statistics
using Plots
pyplot()

include("C:\\Users\\jared\\Desktop\\Github Repos\\DDMR\\Tools\\Wiener Filtering\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
include("..\\DataGen.jl")
include("..\\RedModRun.jl")
include("C:\\Users\\jared\\Desktop\\Github Repos\\DDMR\\Tools\\Model_Reduction_Dev.jl")
## Preference parameters
t_start = 0
t_stop = 10^3
steps = 10^6 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [.35]
sigma_v = sigma
d = 1

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)
Psi(x) = [x; x.^3]

M_out = 25


## Derivitive parameters
steps_tot = steps + discard

nu = size(Psi(0),1)
h = (t_stop - t_start)/(steps)


## Generate signal and predictors
signal = DataGen_Langevin_FE(steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx,
    SM1 = true)
# Notice since SM! is true the output series will include one
# presample. Thus effecting an offset by one in the index.
# Signal is one time step behind, the desired sig.

sig = signal[:,2:end] # True signal


plot(sig[1,:])

pred = zeros(nu, steps)
for n = 1:steps
    pred[:,n] = Psi(signal[:,n])
end

h_ideal = zeros(1,2,20);
h_ideal[1,:,1] = [1+h -h]

h_wf = vector_wiener_filter_fft(pred, sig, M_out)
h_wp = vector_wiener_predictor_fft(pred, sig, M_out)
h_old = h_wf

errh = sum(abs.(h_wf .- h_old),dims = 3)[:,:,1]
sum(h_wp[1,1,:])
sum(h_wf[1,1,:])


p1 = plot(real.([h_wf[1,1,:] h_old[1,1,:]]),
        marker = (2,[:h]),
        label = ["filter" "Predictor"])
p2 = plot(real.([h_wf[2,1,:] h_old[2,1,:]]),
        marker = (2,[:h]),
        label = ["filter" "Predictor"])

plot(p1,p2,layout = (2,1))

sig_hat_ideal = zeros(d,steps)
sig_hat_ideal[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_ideal[:,i] = sum([real.(h_ideal)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

sig_hat_f = zeros(d,steps)
sig_hat_f[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_f[:,i] = sum([real.(h_wf)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

sig_hat_p = zeros(d,steps)
sig_hat_p[1,1] = sig[1,1]
for i = 2 : steps
    sig_hat_p[:,i] = sum([real.(h_wp)[:,:,k]*pred[:,i-k] for k = 1:min(i - 1,M_out)])
end

win = 10080:10100

plot([sig[1,win] sig_hat_f[1,win] sig_hat_p[1,win] sig_hat_ideal[1,win]],
    label = ["sig" "sig_hat_f" "sig_hat_p" "sig_hat_ideal"])

err_ideal = sig .- sig_hat_ideal
err_f = sig .- sig_hat_f
err_p = sig .- sig_hat_p

lags = -15:10

C_ideal = Crosscov(err_ideal,pred,lags)
C_f = Crosscov(err_f,pred,lags)
C_p = Crosscov(err_p,pred,lags)

plot(lags, [C_ideal[1,1,:] C_f[1,1,:] C_p[1,1,:]],
    color=[:red :black :grey],
    line=(2,[:solid :dash :dot]),
    label=["sig_hat_ideal" "sig_hat_f" "sig_hat_p"],
    marker=([:hex :d],6))
