using Statistics
using StatsBase
using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
include("..\\DataGen.jl")

## Preference parameters
t_start = 0
t_stop = 10^4
steps = 10^6 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = randn(3,1)
sigma = I
sigma_v = sigma
d = 3

Nen = 100

Psi(x) = [x;
          0; x[1]*x[3]; 0;
          0; 0; x[1]*x[2] ]

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
    d = d,
    SM1 = true,
    Obs_noise = false
    )

sig = signal[:,2:end]
d, steps = size(sig)
nu = size(Psi(zeros(d,1)),1)

pred = zeros(nu, steps)
for n = 1:steps
    pred[:,n] = Psi(signal[:,n])
end

function Autocov_sm(pred; L = 55, win = "Par")

    lags = 0:L

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
            R_pred_smoothed[i,j,:] = lam .* my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
        end
    end
    R_pred_smoothed
end

P = Autocov_sm(pred)

function Autocov_sm_Alt(pred, lags; win = "Par")
    pred_L = size(pred,2) # I flatten pred to (nuxd) x steps
    m = length(lags)

    L = maximum(lags)
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

    pred .+= mean(pred,dims = 2)

    R_pred_smoothed = zeros(nu,nu,m)
    for k = 1:m

    for i = 1 : nu
        for j = 1 : nu
            R_pred_smoothed[i,j,:] = lam .* my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
        end
    end
    R_pred_smoothed
end


eigen(P[:,:,1])

Rr = P[:,:,1]

cond(Rr)

mean(Pred,dims = 2)
Pred = pred
Pred .-= mean(Pred,dims = 2)

D = Pred[:,1]*Pred[:,1]'

Pred[:,1]
