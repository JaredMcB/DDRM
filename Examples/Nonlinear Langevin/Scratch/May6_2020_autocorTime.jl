"""
The purpose of this script is to investigate order mismatch.
perticularly here we investigate the autocorrelation time and use it
to be sure we have enough samples.
"""

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

gap = 2^5

dt = 2^-10
Δt = gap*dt

Time = range(t_start,t_stop, step = dt)
N_grid = range(t_start,t_stop, step = Δt)

discard_N = 10^3
discard_T = 508*10^3 #gap*discard_N

T = length(Time)
N = length(N_grid)

dW = randn(1,T + discard_T)
W = cumsum(dW,dims = 2)
ΔW = zeros(d,N + discard_N)
for i = 1:(N + discard_N - 1)
    ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,gap*i] - W[:,gap*(i-1)+1])
end

@time signal_T = DataGen_DWOL(T,
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

@time τ_exp, τ_int, P = auto_times(signal_T[1,:],plt = true)


## Finding exponetial autocorrelation time
lags = 0:10^6
A = autocov_con(signal_T[1,:],lags)
plot(A)

end_ind = Int64(round( (findall(A.<0)[1] - 1)/2 ))

plot(log.(A[1:end_ind]))

end_ind = 60000

@time A_mat = [ones(end_ind,1) reshape(1:end_ind,end_ind,1)]
@time b = inv(A_mat'*A_mat)*A_mat'*log.(A[1:end_ind])

b[2]
τ_exp = -1/b[2]
plot!(1:end_ind,A_mat*b)
N_disc = 20*τ_exp

## Integrated autocorrelation time

plot(A)
end_ind_int = findall(A.<0)[1] - 1
plot(A[1:end_ind])


A ./= A[1]

τ_int = .5 + sum(A[2:end_ind_int])
N_eff = T/τ_int
