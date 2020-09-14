

using JLD
using Printf
using StatsBase
using Statistics
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
pyplot()

cd("C:\\Users\\JaredMcBride"*
    "\\Desktop\\Git Stuff\\ModelReduction\\"*
    "Spectral Methods\\Model Reduction Julia Scripts\\Scratch")

data_AA = load("data\\Autocov_ensamble.jld")
AA = data_AA["dat_AA"]

maxlags, bN, l_indx = size(AA)

AA_mean = mean(AA[:,:,:],dims = 2)
AA_var = var(AA[:,:,:],dims = 2)

## Ploting the mean Autocovariance functions together
#      with error bars.
Labels = ["original" "h_m-1" "h_m-2" "h_m-3" "h_m-5" "h_m-7" "h_m-10" "h_m-15" "h_m-20"]
line_sty = [:solid :dash :dashdot :dot :solid :dash :dashdot :dot :solid]

σ = sqrt.(AA_var[1:10:2*10^4,1,1]/bN)
P = plot(AA_mean[1:10:2*10^4,1,1],
    grid = false,
    ribbon = 3*σ,
    label = Labels[1],
    line =  (4, line_sty[1]),
    fillalpha = .5,
    title = "Autocovariance functions of reduced models with varying truncations of Wiener filter")
[1,2,3,4,5,6, 7, 8, 9]
[o,1,2,3,5,7,10,15,20]

for i in [7,8,9]
    σ = sqrt.(AA_var[1:10:2*10^4,1,i]/bN)
    plot!(P,AA_mean[1:10:2*10^4,1,i],
        grid = false,
        ribbon = 3*σ,
        label = Labels[i],
        line =  (4, line_sty[i]),
        fillalpha = .5)
end
gui()
