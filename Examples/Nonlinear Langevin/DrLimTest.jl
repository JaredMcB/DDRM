using PyPlot
using Random

include("DataGen.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")
include("../../Tools/KLPowerSpec.jl")

import .KLPowerSpec

steps = 10^6 + 1
scheme = "FE"
t_start = 0
t_stop = 10^4
discard = 100000
sig_init = [1.5]
sigma = [.3]
V_prime = x -> -x.*(x.^2 .- 1)
SM1 = false
Obs_noise = false
d = 1
#e = randn(d,steps + discard)

# Get full model run
Random.seed!(2014)
X = DataGen_DWOL(
    steps;
    scheme, t_start, t_stop, discard,
    sig_init , sigma, V_prime,
    SM1, Obs_noise, d
    )

T = range(t_start,stop = t_stop, length = steps)
X

data = Dict("steps" => steps,
            "t_stop" => t_start,
            "sigma" => sigma,
            "X" => X)
save("Examples/Nonlinear Langevin/data/data_10_23_2020.jld",data)

data = load("Examples/Nonlinear Langevin/data/data_10_23_2020.jld")
X = data["X"]

auto_times(X[1,:])

# Put in Psi functions
Psi(x) = [x; x.^3]

# Model reduction Parameters
M_out = 50
n = 3
p = 500
par = 55
ty = "bin"
xspec_est = "DM"
rl = true
Preds = false
PI = false
rtol = 1e-6


### Varing parameters
###            xspect_est, par    , nfft    , n    , p
#
Parms = [["DM"       , 5000  , 2^17    , 2    , 5],
         ["SP"       , 5000  , 2^17    , 2    , 5],
         ["DM"       , 100    , 2^10    , 3    , 500],
         ["DM"       , 500    , 2^16    , 3    , 500],
         ["DM"       , 1000   , 2^16    , 3    , 500],
         ["DM"       , 5000   , 2^16    , 3    , 500],
         ["DM"       , 10000  , 2^16    , 3    , 500],
         ["DM"       , 5000   , 2^20    , 3    , 500],
         ["DM"       , 10000  , 2^20    , 3    , 500],
         ["SP"       , 10000  , 2^10    , 2    , 5],
         ["SP"       , 10000  , 2^16    , 2    , 5],
         ["SP"       , 10000  , 2^17    , 2    , 5],
         ["SP"       , 10000  , 2^17    , 3    , 10]]

nfft      = Parms[1][3]

P = 2#length(Parms)

h_wf_packs  = []
times = zeros(P)
for i = 1:P
    Out = @timed get_wf(X, Psi;
        M_out, ty, rl, Preds, PI, rtol, info = true,
        xspec_est = Parms[i][1],
        par       = Parms[i][2],
        nfft      = Parms[i][3],
        n         = Parms[i][4],
        p         = Parms[i][5]);

    append!(h_wf_packs,Out.value)
    times[i]      = Out.time
end

h_wf_dm = h_wf_packs[1]
h_wf_sp = h_wf_packs[8]

h_wf_dm[1,:,1]
h_wf_sp[1,:,1]
