
using LinearAlgebra
using DSP
using StatsBase
using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()

include("..\\DataGen.jl")


t_start = 0
t_stop = 7.2 * 10^5

sig_init = [1.5]
sigma = [.3]
sigma_v = sigma
d = 1

dt = 2^-5
Time = range(t_start,t_stop, step = dt)
discard_T = 508*10^3 #gap*discard_N
T = length(Time)

dW = randn(1,T + discard_T)

signal_T = DataGen_DWOL(T,
    scheme = "FE",
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    sigma = sigma,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = dW)

lags = 0:500000

A = @timed autocov(signal_T[1,:],lags)
@time A_con = autocov_con(signal_T[1,:],lags)
plot!(lags,A_con,label = "A_con")

x = signal_T[1,:]
lx = size(x,1)
A = conv(x,reverse(x))/lx
A = [A[k + lx] for k in lags]

plot(A)
plot!(A)

function autocov_con(x::AbstractVector{<:Real},lags::UnitRange{Int})
    lx = size(x,1)
    x .-= mean(x)
    A = conv(x,reverse(x))/lx
    A = [A[k + lx] for k in lags]
end


x = rand(100)
x[1:95] .+= x[6:100]
x
A = autocov(x)

plot(A)

A_con = autocov_con(x,0:20)

plot!(A_con)
