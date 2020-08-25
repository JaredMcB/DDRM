include("modgen_LSDE.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")

using JLD
using PyPlot
using DSP: nextfastfft


A       = reshape([-0.5],1,1)
σ       = reshape([1],1,1)
Xo      = [1]
t_disc  = 1000
gap     = 10
scheme  = "EM"
d       = size(A,1)
t_start = 0
t_stop  = 1e6
h       = 1e-2
Δt      = h*gap
M_out   = 100

@time X = modgen_LSDE(t_start,t_stop,h,
    A = A,
    σ = σ,
    Xo = Xo,
    t_disc = t_disc,
    gap = gap,
    scheme = scheme)

N       = size(X,2)
nfft    = nextfastfft(N)
X = [X zeros(d,nfft - N)]

τ_exp, τ_int    = auto_times(X[:])*Δt
N_eff           = N*Δt/τ_int

println("Time to get h_wf: ")
Psi(x) = x
@time h_wf_num = get_wf(X,Psi, M_out = M_out)

h_wf_ana = zeros(1,1,M_out)
h_wf_ana[1,1,1] = (1 .+ h*A)[1]

err     = abs.(h_wf_ana - h_wf_num)
err_sum = sum(err)
