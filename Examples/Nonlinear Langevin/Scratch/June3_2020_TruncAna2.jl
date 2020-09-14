"""
Here we show more carefully that autocovariance functions
Improve with more of the Weiener filter coefficients. Unlike
the similar analysis on May 19, 2020 we use rejection
sampling to speed up the batch simulations. This is done by
drawing samples from the equilibrium distribution and using
these as initial values. Thius we iliminate the need to
discard.
"""


using JLD
using Printf
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
pyplot()

cd("C:\\Users\\JaredMcBride"*
    "\\Desktop\\Git Stuff\\ModelReduction\\"*
    "Spectral Methods\\Model Reduction Julia Scripts\\Scratch")

include("..\\DataGen.jl")
include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")

## Run Parameters
t_start = 0
t_stop = 15000

sig_init = [1.5]
sigma = [.35]
sigma_v = sigma
d = 1

gap = 10^2

dt = 10^-3
Δt = gap*dt

Time = range(t_start,t_stop, step = dt)
N_grid = range(t_start,t_stop, step = Δt)

discard_N = 0
discard_T = 0 # This is to implement the sampling from the equil. dist.

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

data_dict = load("data\\Mean_Wiener_filteer_100_May15_2020.jld")

h_wf_ens = data_dict["dat_h_wf_ens"]
data_dict["scheme"]

ANA = analyse_h_ens(h_wf_ens; plt = true)
h_m_EM = ANA[1]
# Here I extract the mean filter
ANA[3]

## Truncate mean filter
h_trunc = VarTrunc(h_m_EM)

## Compute the autocovariances

N
maxlag = N-1
indx = [1,2,3,5,7,10,15,20]
bN = 100
AA      = zeros(maxlag+1,bN,length(indx)+1)
AA_mean = zeros(maxlag+1,length(indx)+1)
AA_var  = zeros(maxlag+1,length(indx)+1)

@time for k = 1:bN
    signal = DataGen_DWOL(T,
        scheme = scheme,
        t_start = t_start,
        t_stop = t_stop,
        discard = discard_T,
        sig_init = sig_init,
        sigma = sigma,
        V_prime = x -> -x.*(x.^2 .- 1),
        SM1 = false)
    AA[:,k,1] = autocov_con(signal[1,1:gap:end],0:10^6)
end

AA_mean[:,1] = mean(AA[:,:,1],dims = 2)
AA_var[:,1] = var(AA[:,:,1],dims = 2)

@time for i = 1:length(indx)
    @time AA[:,:,1+i] = Autocov(h_trunc[:,:,:,indx[i]], Psi,
        N,
        lag = 0:maxlag,
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

# Save the autocaovariance functions. We save them all to
#   be able to estimate error bars.

dat = Dict(
    "bN" => bN,
    "dat_AA" => AA)

save("data\\Autocov_ensamble.jld",
    merge(params,dat))

## Load Data
# Loading the data in the future
data_AA = load("data\\Autocov_ensamble.jld")
AA = data_AA["dat_AA"]

maxlags, bN, l_indx = size(AA)

AA_mean = mean(AA[:,:,:],dims = 2)
AA_var = var(AA[:,:,:],dims = 2)

## Ploting the mean Autocovariance functions together
#      with error bars.
Labels = ["original" "h_m-1" "h_m-2" "h_m-3" "h_m-5" "h_m-7" "h_m-10" "h_m-15" "h_m-20"]
line_sty = [:solid :dash :dash :dash :dashdot :dashdot :dot :dot :dot :dot]

σ = sqrt.(AA_var[1:2*10^4,1,1]/bN)
plot(AA_mean[1:2*10^4,1,1],
    grid = false,
    ribbon = 3*σ,
    label = Labels[1],
    line =  (4, line_sty[1]),
    fillalpha = .5,
    title = "Autocovariance functions of reduced models with varying truncations of Wiener filter")

for i = 2:l_indx
    σ = sqrt.(AA_var[1:2*10^4,1,i]/bN)
    plot!(AA_mean[1:2*10^4,1,i],
        grid = false,
        yerror = 3*σ,
        label = Labels[i],
        line =  (4, line_sty[i]),
        fillalpha = .5)
end

plot!(title = "Autocovariance functions of reduced models with varying truncations of Wiener filter")
gui()
