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

Psi(x) = [diagm(x) [0; x[1]*x[3];0] [0;0;x[1]*x[2]]]

M_out = 10

# Derivitive parameters
steps_tot = steps + discard
nu = size(Psi(zeros(d)),2)
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

function get_pred_M(sig, Psi)
    d, steps = size(sig)
    d, nu = size(Psi(zeros(d)))

    pred = zeros(d, nu, steps)
    for n = 1:steps
        pred[:,:,n] = Psi(sig[:,n])
    end
    pred
end

pred = get_pred_M(sig,Psi)
Nex = 2^10
win = "Par"

d, stepsy = size(sig)
dp, nu, stepsx = size(pred)

d == dp || error("first dimension of pred should equal dimension of sig.")

stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])

Nexh = Int(floor(Nex/2))

L = 50
lags = 0:L;

# Smoothed viewing window
if win == "Bar"
    lam = 1 .- (0:L)/L
elseif win == "Tuk"
    lam = .5*(1 .+ cos.(pi/L*(0:L)))
elseif win == "Par"
    LL = Int(floor(L/2))
    lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
    lam2 = 2*(1 .- (LL+1:L)/L).^3
    lam = [lam1; lam2]
else
    lam = ones(L+1)
end

R_pred_smoothed = zeros(nu,nu,length(lags))
for i = 1 : nu
    for j = 1 : nu
        R_pred_smoothed[i,j,:] = lam .* sum([my_crosscov(pred[k,i,:],pred[k,j,:],lags) for k = 1:d])
    end
end

l = Matrix_CKMS_c(R_pred_smoothed);

l_pad_minus = cat(dims = 3,l,zeros(nu,nu,Nex - L - 1))

z_spect_pred_minus_num_fft = fft(l_pad_minus,3)
z_spect_pred_plus_num_fft =complex(zeros(nu,nu,Nex))
for i = 1 : Nex
    z_spect_pred_plus_num_fft[:,:,i] = z_spect_pred_minus_num_fft[:,:,i]'
end

z_crossspect_predsig_num_fft = z_crossspect_fft(pred,sig)

# This computes the impule response (coefficeints of z) for S_{yx}{S_x^+}^{-1}
S_predsig_overS_plus_fft_num = complex(zeros(nu,Nex))
for i = 1: Nex
    S_predsig_overS_plus_fft_num[:,i] = z_spect_pred_plus_num_fft[:,:,i]\z_crossspect_predsig_num_fft[:,i]
end

S_predsig_overS_plus_fft_num_fft = ifft(S_predsig_overS_plus_fft_num,2)

# Extracts causal part coefficinets of S_{yx}{S_x^+}^{-1}, {S_{yx}{S_x^+}^{-1}}_+
S_predsig_overS_plus_fft_plus_num_fft = cat(dims = 2,S_predsig_overS_plus_fft_num_fft[:,1: Nexh], zeros(nu,Nex - Nexh))

# Computes causal part of S_{yx}/S_x^+, {S_{yx}/S_x^+}_+
S_predsig_overS_plus_plus_num_fft = fft(S_predsig_overS_plus_fft_plus_num_fft,2);

# Obtain transfer function H by dividing {S_{yx}/S_x^+}_+ by S_x^-

H_num = complex(zeros(nu,Nex))
for i = 1: Nex
    H_num[:,i] = z_spect_pred_minus_num_fft[:,:,i]\S_predsig_overS_plus_plus_num_fft[:,i]
end

# Extrct tranferfunction coeffifcients (impulse responce of Weiner filter)
h_num_raw = ifft(H_num,2)

# Truncate
h_num_fft = h_num_raw[:,1:M_out]

h_num_fft[:,1]
h_wf = h_num_fft

##Test it
v = randn(d,steps_tot)
signal_rm = zeros(d,steps_tot)
signal_rm[:,1] = sig_init
PSI = zeros(d,nu,steps_tot)
for i = 1: steps_tot-1
    PSI[:,:,i] = Psi(signal_rm[:,i])
    signal_rm[:,i+1] = sum([PSI[:,:,i-k]*real.(h_wf)[:,k+1]
        for k = 0:min(i - 1,M_out - 1)]) + sqrt(dt)*sigma_v*v[:,i+1]
end

sig_rm = signal_rm[:,discard+1:steps_tot]

plot(tim,sig_rm[:, :]')

plot(sig[1,wind],sig[2,wind],sig[3,wind],
    seriestype = :scatter)
