using Random
using Distributions

# Get software to generate model
gen = include("../../Examples/Nonlinear Langevin/DataGenDWOL.jl")

# Get model reduction software being tested
mr  = include("../Model_Reduction_Dev.jl")
at  = include("../AnalysisToolbox.jl")

# Model run Parameters
sigma    = [.4]
sig_init = [1.5]
# Numerical estimate parameters
scheme   = "FE"
steps    = 10^7 # Number of time steps (not including those discarded)
h        = .01
discard  = steps # Number of time steps discarded
gap      = 1

V_c_prime  = x -> -x.*(x.^2 .- 1)

h_wf = zeros(1,3,1)
h_wf[1,:,1] = [1.01 0 -0.01]

Psi(x)  = [x; x.^2; x.^3]

noise = true
noise_dist = MvNormal(zeros(1),sqrt(h)*sigma)
rand(noise_dist)

Random.seed!(2014)
X_rm_c = real(mr.redmodrun(reshape(sig_init,1,:), h_wf, Psi; steps,discard, noise, noise_dist))

Random.seed!(2014)
steps_tot = steps + discard
e = 1/sqrt(h)/sigma[1]*rand(noise_dist,steps_tot)

X_c = @time gen.DataGen_DWOL(; sigma, V_prime = V_c_prime, sig_init, scheme, steps, h, discard, gap, ObsNoise = true,e)

sum(abs.(X_rm_c[1,2:end] - X_c[1][1,1:end-1]))
