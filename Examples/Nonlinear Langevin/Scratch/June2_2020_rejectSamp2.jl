using JLD
using Printf
using Distributions
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
pyplot()

cd("C:\\Users\\JaredMcBride"*
    "\\Desktop\\Git Stuff\\ModelReduction\\"*
    "Spectral Methods\\Model Reduction Julia Scripts\\Scratch")


include("..\\DataGen.jl")
include("..\\DWOL_eqidist_sampler.jl")


## Run Parameters

t_start = 0
t_stop  = 1500000

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


## Generate looong time series


dW = randn(1,T + discard_T)
W = cumsum(dW,dims = 2)
ΔW = zeros(d,N + discard_N)
for i = 1:(N + discard_N - 1)
    ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,gap*i] - W[:,gap*(i-1)+1])
end


@time signal_T = DataGen_DWOL(T,
    scheme = scheme,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    sigma = sigma,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = dW)

# Save looong time series
dat = Dict(
    "dat_signal_T" => signal_T[1,:])

save("data\\full_model_run.jld",
    merge(parameters,dat))

# load looong time series
data_dict = load("data\\full_model_run.jld")
signal_T = data_dict["dat_signal_T"]

@time τ_exp, τ_int, P = auto_times(signal_N,plt = true)
τ_exp
τ_int
N_disc = 20*τ_exp
N_eff = N/τ_int

@time τ_exp, τ_int, P = auto_times(signal_T[1,1:10^8],plt = true)
τ_exp
τ_int
N_disc = 20*τ_exp
N_eff = T/τ_int



pdf_N, b_mpts, P  = emp_pdf(signal_T, plt_only = false)

# Save the pdf of Looong timeseires
dat = Dict(
    "dat_pdf_N" => pdf_N,
    "dat_b_mpts" => b_mpts)

save("data\\full_model_run_pdf.jld",
    merge(parameters,dat))

# Load Pdf of Looong timeseires


Z = DWOL_dist_samp(10^6)
P_o = emp_pdf(Z)

plot(b_mpts,pdf_N)

x = -2:0.01:2

μ_pos = 1
μ_neg = -1
m = maximum(pdf_N)*1.05
σ = 1/(sqrt(2π)*m)

D_neg = Normal(μ_pos,σ)
D_pos = Normal(μ_neg,σ)


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
