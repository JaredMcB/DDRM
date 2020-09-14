"""
Here we investigate sampling from the equilibrium distirbution
using rejection sampling.
"""


using JLD
using Printf
using Distributions
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
pyplot()

include("..\\DataGen.jl")
include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")

cd("C:\\Users\\JaredMcBride"*
    "\\Desktop\\Git Stuff\\ModelReduction\\"*
    "Spectral Methods\\Model Reduction Julia Scripts\\Scratch")

t_start = 0
t_stop  = 150000

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

parameters = Dict(
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
    "nu" => nu,
    "M_out" => M_out,
    "Nen" => Nen,
    "scheme" => scheme
)


## Load Mean filter

data_dict = load("data\\Mean_Wiener_filteer_100_May15_2020.jld")

h_wf_ens = data_dict["dat_h_wf_ens"]
data_dict["scheme"]

# dW = randn(1,T + discard_T)
# W = cumsum(dW,dims = 2)
# ΔW = zeros(d,N + discard_N)
# for i = 1:(N + discard_N - 1)
#     ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,gap*i] - W[:,gap*(i-1)+1])
# end
#
# @time signal_T = DataGen_DWOL(T,
#     scheme = scheme,
#     t_start = t_start,
#     t_stop = t_stop,
#     discard = discard_T,
#     sig_init = sig_init,
#     sigma = sigma,
#     V_prime = x -> -x.*(x.^2 .- 1),
#     SM1 = false,
#     Obs_noise = false,
#     e = dW)
#
# dat = Dict(
#     "dat_signal_N" => signal_T[1,1:100:end])
#
# save("data\\full_model_run.jld",
#     merge(parameters,dat))

data_dict = load("data\\full_model_run.jld")
signal_N = data_dict["dat_signal_N"]

@time τ_exp, τ_int, P = auto_times(signal_N,plt = true)
τ_exp
τ_int

pdf_N, b_mpts, P  = emp_pdf(signal_N, plt_only = false)
P
function p(z)
    i = argmin(abs.(b_mpts .- z))
    if (i == 1) | (i == length(b_mpts))
        x = [1; 2; 3]
        y = zeros(3)
    else
        x = b_mpts[i-1:i+1]
        y = pdf_N[i-1:i+1]
    end
    a = [x.^2 x ones(3)] \ y
    P = a[1]*z^2 + a[2]*z + a[3]
end

plot!(b_mpts,pdf_N,t=:scatter,
    marker = (3,:h))

μ_pos = 1
μ_neg = -1
m = maximum(pdf_N)*1.05
σ = 1/(sqrt(2π)*m)

D_neg = Normal(μ_pos,σ)
D_pos = Normal(μ_neg,σ)


x = -2:0.01:2
Y_pos = map(z -> pdf(D_pos,z),x)
Y_neg = map(z -> pdf(D_neg,z),x)
Y = (Y_neg .+ Y_pos)/2

plot!(x,2*Y, line = :dash)

q_neg(x) = 1/(sqrt(2π)*σ)*exp(-(x-μ_neg)^2/(2σ^2))
q_pos(x) = 1/(sqrt(2π)*σ)*exp(-(x-μ_pos)^2/(2σ^2))
q(x) = (q_neg(x) + q_pos(x))/2
cq(x) = 2*q(x)

plot!(x,cq.(x), line = (2,:dot))

N_samp = 4*10^3
Z = [(rand() < .5 ? rand(D_neg) : rand(D_pos)) for i = 1 : N_samp]

P_z = emp_pdf(Z)

ZZ = [[z,rand()*cq(z)] for z in Z]


plot!([ZZ[1] for ZZ in ZZ],[ZZ[2] for ZZ in ZZ],t = :scatter,
    color = :red,
    markeralpha = 0.3,
    markerstrokewidth=0)

plot(b_midpts,pdf)

z = filter(z -> z[2] <= p(z[1]) ,ZZ)
plot!([z[1] for z in z],[z[2] for z in z],t = :scatter,
    color = :green,
    markeralpha = 0.3,
    markerstrokewidth=0)

function samp()
    z = (rand() < .5 ? rand(D_neg) : rand(D_pos))
    zz = [z,rand()*cq(z)]
    while zz[2] > p(zz[1])
        z = (rand() < .5 ? rand(D_neg) : rand(D_pos))
        zz = [z,rand()*cq(z)]
    end
    z
end

o = samp()
O = [samp() for i = 1:100000]

P_o = emp_pdf(O)
