include("modgen_LSDE.jl")
include("../../Tools/Model_Reduction_Dev.jl")

# include("modgen_LSDE.jl")
# include("..\\..\\Tools\\Model_Reduction_Dev.jl")

using JLD
using PyPlot
using DSP: nextfastfft

function runner(;
    A       = reshape([-0.5],1,1),
    σ       = reshape([1],1,1),
    Xo      = [1],
    t_disc  = 1000,
    gap     = 10,
    scheme  = "EM",
    d       = size(A,1),
    t_start = 0,
    t_stop  = 1e6,
    h       = 1e-2,
    Δt      = h*gap,
    M_out   = 100)
    println("===========t_stop = $t_stop===========")
    println("Time to get data: ")
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
    comp1_err = err[1,1,1]
    tail_err = maximum(err[:,:,2:end])
    err_sum = sum(err)

    println("======================")
    println("N_eff : $N_eff")
    println("comp1_err : $comp1_err")
    println("tail_err : $tail_err")
    println("======================")
    [N_eff; comp1_err; tail_err]
end


T_stop = map(x -> 10^x, 4:.5:7)
data = zeros(3,length(T_stop))
M = 20

for i in 1:length(T_stop)
    for j in 1:M
        dat = runner(t_stop = T_stop[i])
        data[:,i] .+=  dat/M
    end
end

save("/u5/jaredm/data/LSDE_Data/NoiseVNeff.jld", "data", data)
# save("c:\\Users\\JaredMcBride\\Desktop\\DDMR\\Examples\\LinearSDE\\LSDE_Data\\NoiseVNeff.jld", "data", data)
