using Printf
using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


include("..\\DataGen.jl")
include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")



t_start = 0
t_stop = 10000

sig_init = [1.5]
sigma = [.3]
sigma_v = sigma
d = 1

# gap = 10^2

dt = 10^-3
# Δt = gap*dt

Time = range(t_start,t_stop, step = dt)
# N_grid = range(t_start,t_stop, step = Δt)

# discard_N = 10^3
discard_T = 10^7 # gap*discard_N

T = length(Time)
# N = length(N_grid)

dW = randn(1,T + discard_T)
W = cumsum(dW,dims = 2)
# ΔW = zeros(d,N + discard_N)
# for i = 1:(N + discard_N - 1)
#     ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,gap*i] - W[:,gap*(i-1)+1])
# end

ANS = @timed DataGen_DWOL(T,
    scheme = "FE",
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    sigma = sigma,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = dW) # 60sec

signal_T = ANS[1]

data[i,:] = signal_T[1,:]

clock_time = ANS[2]
