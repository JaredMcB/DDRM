"""
The purpose of this script is to investigate order mismatch.
"""

using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


include("..\\DataGen.jl")
include("C:\\Users\\jared\\Desktop\\Github Repos\\DDMR\\Tools\\Wiener Filtering\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")

t_start = 0
t_stop = 100

sig_init = [1.5]
sigma = [.5]
sigma_v = sigma
d = 1

dt = 2^-10
Δt = 2^-5
scale = Int( Δt/dt )

Time = range(t_start,t_stop, step = dt)
N_grid = range(t_start,t_stop, step = Δt)

discard_N = 10^3
discard_T = scale*discard_N


T = length(Time)
N = length(N_grid)

dW = randn(1,T + discard_T)
W = cumsum(dW,dims = 2)

signal_T = DataGen_DWOL(T,
    scheme = "FE",
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    sigma = [.5],
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = dW)


ΔW = zeros(d,N + discard_N)
for i = 1:(N + discard_N - 1)
    ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,scale*i] - W[:,scale*(i-1)+1])
end
signal_N = DataGen_DWOL(N,
    scheme = "RK4",
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_N,
    sig_init = sig_init,
    sigma = [.5],
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = ΔW)

plot(Time,signal_T[1,:])
plot!(N_grid,signal_N[1,:])


## Forward Euler

Psi1(x) = [x ; x.^3]


plot(Time[1:2^10*5],signal_T[1,1:2^10*5])
scatter!(N_grid[1:2^5*5],signal_T[1,1:scale:2^10*5], marker = (4,:d))

h_wf = get_wf(signal_T[:,1:scale:end], Psi1)

h_ideal = zeros(1,2,1)

sig_hat_m = redmodrun(
    , Psi1;
    sigma_v = sigma_v,
    steps = T,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    obs_noise = dW
    )

plot!(Time[1:2^10*5],sig_hat_m[1,1:2^10*5],
    line = (3,:dot))


plot!(N_grid[1:2^5*5],sig_hat_m_trunc[1,1:2^5*5],
    marker = (4,:h),
    line = (3,:dot))

h_ens = Run_and_get_WF_DWOL(
    Psi1,
    scheme = "FE",
    steps = N*10^3,
    t_start = t_start,
    t_stop = t_stop*10^3,
    discard = discard_N,
    sig_init = sig_init,
    sigma = sigma,
    d = 1,
    V_prime = (x -> -x.*(x.^2 .- 1)),
    Nen = 100)

ANA = analyse_h_ens(h_ens; plt = true)
h_m = ANA[1]

## this is not done right see May 5th for a better experiment.
