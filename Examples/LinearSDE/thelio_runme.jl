using JLD
using PyPlot
using DSP: nextfastfft

include"../../Tools/modgen_LSDE.jl")
include("../../Tools/Model_Reduction_Dev.jl")

#Parameters
A       = reshape([-0.5],1,1)
σ       = reshape([1],1,1)
Xo      = [1]
t_disc  = 1000
gap     = 10
scheme  = "EM"
d       = size(A,1)
t_start = 0
t_stop  = 1e7
h       = 1e-2
Δtry    = h*gap

# Study the convergence of the wiener filter.

T_stop = [1e4, 1e5, 1e6, 1e7]

for i = 1:length(T_stop)
    println("Generating data for t_stop = $T_stop[i]")
    println("it took :")
    @time X = modgen_LSDE(t_start,t_stop,h,
        A = A,
        σ = σ,
        Xo = Xo,
        t_disc = t_disc,
        gap = gap,
        scheme = scheme)

    # Groom the data
    N = size(X,2)
    nfft = nextfastfft(N)
    X = [X zeros(d,nfft - N)]

    τ_exp, τ_int = auto_times(X1[:])*Δt
    N_ef = floor(N*Δt/τ_int)
    println("tau_exp = $τ_exp")
    println("tau_int = $τ_int")
    println("N_ef = $N_ef")

    Psi(x) = x

    @time h_wf = get_wf(X,Psi)

    
