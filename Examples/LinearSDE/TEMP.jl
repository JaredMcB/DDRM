# Test save seed

include("modgen_LSDE.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")

using Random

A       = reshape([-0.5],1,1)
σ       = reshape([1],1,1)
Xo      = [1]
t_disc  = 100
gap     = 10
d       = size(A,1)
t_start = 0
t_stop  = 1e3
h       = 1e-2
Δt      = h*gap
M_out   = 100

println("===========t_stop = $t_stop===========")
println("Time to get data: ")


seed = Random.seed!(seed).seed

@time X = modgen_LSDE(t_start,t_stop,h;
    A,
    σ,
    Xo,
    t_disc,
    gap)

X2 = copy(X)

X1 - X2
