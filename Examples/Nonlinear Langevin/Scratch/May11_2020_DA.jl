"""
The purpose of this script is to investigate order mismatch.
"""

using JLD
using Printf
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
pyplot()


include("..\\DataGen.jl")
include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")


t_start = 0
t_stop = 150000

sig_init = [1.5]
sigma = [.35]
sigma_v = sigma
d = 1

gap = 10^2

dt = 10^-3
Δt = gap*dt

Time = range(t_start,t_stop, step = dt)
N_grid = range(t_start,t_stop, step = Δt)

discard_N = 10^5
discard_T = gap*discard_N

T = length(Time)
N = length(N_grid)

Psi(x) = [x; x.^3]

## Mean filter info
nu = 2
M_out = 20
Nen = 100
scheme = "EM"

params = Dict(
    "t_start" => t_start,
    "t_stop" => t_stop,
    "sig_init" => t_stop,
    "sigma" => sigma,
    "sigma_v" => sigma_v,
    "d" => d,
    "gap" => gap,
    "dt" => dt,
    "Δt" => Δt,
    "discard_N" => discard_N,
    "discard_T" => discard_T,
    "Psi" => "Double Well potential",
    "nu" => 2,
    "M_out" => 20,
    "Nen" => 40,
    "scheme" => "EM"
)


## Single run
dW = randn(1,T + discard_T)
W = cumsum(dW,dims = 2)
ΔW = zeros(d,N + discard_N)
for i = 1:(N + discard_N - 1)
    ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,gap*i] - W[:,gap*(i-1)+1])
end

ANS = @timed DataGen_DWOL(T,
    scheme = "EM",
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    sigma = sigma,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = dW) # 120sec

signal_T = ANS[1]
clock_time = ANS[2]

τ_exp, τ_int, P = @time auto_times(signal_T[1,:], plt = true) # 7 sec

τ_exp*dt
τ_int*dt

N_disc = 20*τ_exp
N_eff = T/τ_int

@time h_wf = get_wf(signal_T[:,1:gap:end],Psi)

dat = Dict(
    "dat_signal_T" => signal_T,
    "dat_dW" =>dW,
    "dat_h_wf" =>h_wf)

save("Model Reduction Julia Scripts\\Scratch\\data\\Singal_Run_May15_2020.jld",merge(params,dat))


## Mean filter


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
    @printf("done with run %i\n",k)
    h_wf_ens[:,:,:,k] = get_wf(signal[:,1:gap:end],Psi)
end   #764 sec # Nen = 40 9946 sec
h_wf_ens = real.(h_wf_ens)

dat = Dict(
    "dat_h_wf_ens" => h_wf_ens)

save("Model Reduction Julia Scripts\\Scratch\\data\\Mean_Wiener_filteer_100_May15_2020.jld",
    merge(params,dat))

ANA = analyse_h_ens(h_wf_ens; plt = true)
h_m_EM = ANA[1]
# Here I extra the mean filter
ANA[3]
h = abs.(h_m_FE .- h_m_EM)


## Run the reduced model

plot(h_m[1,:,:]')

@time sig_hat_wf = redmodrun(
    h_wf, Psi;
    sigma_v = sigma_v,
    steps = N,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_N,
    sig_init = sig_init,
    obs_noise = ΔW
    )# 25 sec

t_end = 150000
stop_d = Int64(round(t_end/dt))
plot(Time[1:gap:stop_d],signal_T[1,1:gap:stop_d],
    color = :red,
    line = :dot,
    xlabel = "time",
    label = "Full Model",
    title = "Preformance of Low Resolution Reduced model")
stop_Δ = Int64(round(t_end/Δt))
plot!(N_grid[1:stop_Δ],sig_hat_wf[1,1:stop_Δ],
    label = "Reduced Model",
    color = :green,
    line = :dash)

plot(Time[1:gap*10:stop_d],abs.(signal_T[1,1:gap*10:stop_d] .- sig_hat_wf[1,1:10:stop_Δ]),
    color = :red,
    line = :solid,
    xlabel = "time",
    label = "Full Model",
    leg = :none,
    title = "Absolute difference of Low Resolution Reduced Model and Full model")



H_true = my_hist(signal_T[1,1:gap:end],50)
H_rm = my_hist(sig_hat_wf[1,:],50)
plot!(H_true,title = "Histogram of full full-model run")
plot!(H_rm,title = "Histogram of full reduced-model run")

plot(H_true,H_rm)



h_trunc = VarTrunc(h_wf)
indx = [1,2,3,5,10,20]
A = zeros(10^6+1,length(indx)+1)
sig_hat = zeros(N,length(indx)+1)
A[:,1] = autocov_con(signal_T[1,1:gap:end],0:10^6)
sig_hat[:,1] = signal_T[1,1:gap:end]
for i = 1:length(indx)
    @time sig_hat_m = redmodrun(
        h_trunc[:,:,:,indx[i]], Psi;
        sigma_v = sigma_v,
        steps = N,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard_N,
        sig_init = sig_init,
        obs_noise = ΔW
        )#
    sig_hat[:,i+1] = sig_hat_m[1,:]
    A[:,i+1] = autocov_con(sig_hat_m[1,:],0:10^6)
end

plot(A[1:10000,:],
    title = "Autocovariance functions of reduced models with varying truncations of Wiener filter",
    line = (2,[:solid :dash :dashdot :dashdot :dot :dot :dot :dot]),
    label = ["original" "h_wf-1" "h_wf-2" "h_wf-3" "h_wf-5" "h_wf-10" "h_wf-20"])

plot(Time[1:gap:stop_d],signal_T[1,1:gap:stop_d],
    color = :red,
    line = :dot,
    xlabel = "time",
    label = "Full Model",
    title = "Preformance of Low Resolution Reduced model")

plot!(N_grid[1:stop_Δ],sig_hat[1:stop_Δ,6],
    label = "Reduced Model",
    line = :dash)

plot(Time[1:gap*10:stop_d],abs.(signal_T[1,1:gap*10:stop_d] .- sig_hat[1:10:stop_Δ,6]),
    color = :red,
    line = :solid,
    xlabel = "time",
    label = "Full Model",
    leg = :none,
    title = "Absolute difference of Low Resolution Reduced Model and Full model")

indx = [1,2,3,5,10,20]
AA = zeros(10^6+1,length(indx)+1)
AA[:,1] = autocov_con(signal_T[1,1:gap:end],0:10^6)
@time for i = 1:length(indx)
    @time AA[:,1+i] = Autocov(h_trunc[:,:,:,indx[i]], Psi,
        N,
        lag = 0:10^6,
        bN = 10,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard_N,
        sig_init = sig_init,
        sigma = sigma,
        vari = false)
        @printf("Finished with iteration %i\n",i)
end

plot(AA[1:10000,:],
    title = "Autocovariance functions of reduced models with varying truncations of Wiener filter",
    line = (2,[:solid :dash :dashdot :dashdot :dot :dot :dot :dot]),
    label = ["original" "h_m-1" "h_m-2" "h_m-3" "h_m-5" "h_m-10" "h_m-20"])
