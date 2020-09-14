"""
Here we show more carefully that autocovariance functions
Improve with more of the Weiener filter coefficients.
"""


using JLD
using Printf
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
# pyplot()

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


## Load Mean filter

data_dict = load("Model Reduction Julia Scripts\\Scratch\\data\\Mean_Wiener_filteer_100_May15_2020.jld")

h_wf_ens = data_dict["dat_h_wf_ens"]
data_dict["scheme"]

ANA = analyse_h_ens(h_wf_ens; plt = true)
h_m_EM = ANA[1]
# Here I extract the mean filter
ANA[3]

h_trunc = VarTrunc(h_m_EM)

indx = [1,2,3,5,10,20]
bN = 100
AA = zeros(10^6+1,bN,length(indx)+1)
AA_mean = zeros(10^6+1,length(indx)+1)
AA_var = zeros(10^6+1,length(indx)+1)
@time for k = 1:bN
    signal = DataGen_DWOL(T,
        scheme = scheme,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard_T,
        sig_init = sig_init,
        sigma = sigma,
        V_prime = x -> -x.*(x.^2 .- 1),
        SM1 = true)
    AA[:,k,1] = autocov_con(signal[1,1:gap:end],0:10^6)
end

dat = Dict(
    "bN" => bN,
    "dat_AA" => AA)

save("Model Reduction Julia Scripts\\Scratch\\data\\Autocov_ensamble.jld",
    merge(params,dat))

data_AA = load("Model Reduction Julia Scripts\\Scratch\\data\\Autocov_ensamble.jld")
AA = data_AA["dat_AA"]

Labels = ["original" "h_m-1" "h_m-2" "h_m-3" "h_m-5" "h_m-10" "h_m-20"]
line_sty = [:solid :dash :dashdot :dashdot :dot :dot :dot :dot]


AA_mean[:,1] = mean(AA[:,:,1],dims = 2)
AA_var[:,1] = var(AA[:,:,1],dims = 2)
σ = sqrt.(AA_var[1:2*10^4,1]/bN)
plot(AA_mean[1:2*10^4,1],
    grid = false,
    ribbon = 3*σ,
    label = Labels[1],
    line =  line_sty[1],
    fillalpha = .5)

@time for i = 1:length(indx)
    @time AA[:,:,1+i] = Autocov(h_trunc[:,:,:,indx[i]], Psi,
        N,
        lag = 0:10^6,
        bN = bN,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard_N,
        sig_init = sig_init,
        sigma = sigma,
        ALL = true)
    AA_mean[:,1+i] = mean(AA[:,:,1+i],dims = 2)
    AA_var[:,1+i] = var(AA[:,:,1+i],dims = 2)
    @printf("Finished with iteration %i\n",i)
end

dat = Dict(
    "dat_AA" => AA)

save("Model Reduction Julia Scripts\\Scratch\\data\\Autocov_ensamble.jld",
    merge(params,dat))

for i = 1:length(indx)
    σ = sqrt.(AA_var[1:2*10^4,1+i]/bN)
    plot!(AA_mean[1:2*10^4,1+i],
        grid = false,
        ribbon = 3*σ,
        label = Labels[1+i],
        line =  line_sty[1+i],
        fillalpha = .5)
end

plot!(title = "Autocovariance functions of reduced models with varying truncations of Wiener filter")
gui()

i = 1
AA[:,:,1+i] = Autocov(h_trunc[:,:,:,indx[i]], Psi,
    N,
    lag = 0:10^6,
    bN = bN,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_N,
    sig_init = sig_init,
    sigma = sigma,
    ALL = true)
