"""
The purpose of this script is to investigate order mismatch.
"""


using Printf
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


include("..\\DataGen.jl")
include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")



t_start = 0
t_stop = 100000

sig_init = [1.5]
sigma = [.35]
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

Plts = Dict()
N_itr = 10
data = zeros(N_itr,4)
data_sig = zeros(N_itr,T)
@printf("\n\ni | t_stop  | tau_exp | tau_int |  N_disc  | N_eff \n")
@printf(    "-------------------------------------------------")
@time for i = 1:N_itr

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

    data_sig[i,:] = signal_T[1,:]

    clock_time = ANS[2]
    # The goal is to approximate this times series

    # Set the parameters of the batch run to
    # get an ensemble of wiener filters.


    # Check Autocor times
    τ_exp, τ_int, P = auto_times(signal_T[1,:], plt = true) # 7 sec

    Plts[i] = P

    N_disc = 20*τ_exp
    N_eff = T/τ_int
    @printf("%i | %.1f | %.3f | %.3f | %.0f | %.0f \n",i,t_stop,τ_exp*dt, τ_int*dt,N_disc,N_eff)
    plot(data_sig[i,1:50:T],title = "run "*string(i))
    gui()
    data[i,:] = [τ_exp*dt, τ_int*dt,N_disc,N_eff]
end
data_m = mean(data,dims = 1)
@printf("mean: %.3f | %.3f | %.0f | %.0f \n", data_m[1,1], data_m[1,2], data_m[1,3], data_m[1,4])


Psi(x) = [x; x.^3]
nu = 2
M_out = 20
Nen = 10
scheme = "FE"

h_wf_ens = complex(zeros(d,nu,M_out,Nen))
@time for k = 1:Nen
    signal = DataGen_DWOL(T,
        scheme = scheme,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard_T,
        sig_init = sig_init,
        sigma = sigma,
        V_prime = x -> -x.*(x.^2 .- 1),
        SM1 = true)
    print("done with run ",k,"/n")
    h_wf_ens[:,:,:,k] = get_wf(signal[:,1:gap:end],Psi)
end   #764 sec
h_wf_ens = real.(h_wf_ens)

ANA = analyse_h_ens(h_wf_ens; plt = true)
h_m = ANA[1]
# Here I extra the mean filter
ANA[3]

@time sig_hat_m = redmodrun(
    h_m, Psi;
    sigma_v = sigma_v,
    steps = N,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_N,
    sig_init = sig_init,
    obs_noise = ΔW
    )# 65 sec
t_end = 50
plot(Time[1:2^5*t_end],signal_T[1,1:2^5*t_end],
    color = :red,
    line = :solid,
    xlabel = "time",
    label = "Full Model",
    title = "Preformance of Low Resolution Reduced model")
# scatter!(N_grid[1:2^5*t_end],signal_T[1,1:gap:2^10*t_end],
#     marker = (3,:d),
#     color = :black)
plot!(N_grid[1:2^2*t_end],sig_hat_m[1,1:2^2*t_end],
    marker = (2,:h),
    label = "Reduced Model",
    color = :green,
    line = :dash)


H_true = my_hist(signal_T[1,1:gap:end],50)
H_rm = my_hist(sig_hat_m[1,:],50)
plot!(H_true,title = "Histogram of full full-model run")
plot!(H_rm,title = "Histogram of full reduced-model run")

plot(H_true,H_rm)

plot(signal_T[1,:])
