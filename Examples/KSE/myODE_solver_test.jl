"""
Here I want to test my ODE solver.

The first test will be the ODE wll be the second order
d²x/dt² = -x , x(0) = 0, x'(0) = 1
This sound give the sine function.
"""

using PyPlot
using LinearAlgebra
mos = include("myODE_solver.jl")

scheme = mos.scheme_ETDRK4

L = [-1;-.5]
F_etd = [L, x->zeros(Complex{Float64},size(x))]

A = diagm(L)
F(x) = A*x
init = complex([1;1.0])

h       = .1
steps   = ceil(Int,2π/h)
discard = 0
gap     = 1

x = mos.my_ODE_solver(scheme,init;
    F = F_etd,
    steps,
    discard,
    h,
    gap)
x
t       = range(0,2π;length = steps)
plot(t,x[:,:]',"--")
plot(t,[exp.(-t) exp.(-.5*t)])

### True
F_true(t) = exp(A*t)*init
x_true  = zeros(Complex128,2,steps)

for i in 1:steps
    x_true[:,i] = F_true(t[i])
end




n = size(init,1) # Dimension of system
x = zeros(Complex128,n,(steps - 1) ÷ gap + 1)

step! = scheme(F, h)

# main stepping loop
temp = init
for n = 1:steps+discard
    global temp
    if (n > discard) & ((n - discard - 1) % gap == 0)   # save state
            x[:,(n-discard-1)÷gap+1] = temp
    end
    temp = step!(temp,temp)
end
x
init

## # Stiff example
scheme = mos.scheme_ETDRK4

δ = .0001
P = 2/δ

L = [ 1; 1.0]
F_etd = [L, x->[x[1]^2 - x[1]^3- x[1]; x[2]^2 - 2x[2]^3 - x[2]]]

A = diagm(L)
F = x -> [x[1]^2 - x[1]^3; x[2]^2 - 2x[2]^3]
init = complex([δ;δ])

h       = 1
steps   = ceil(Int,P/h)
discard = 0
gap     = 3


x = mos.my_ODE_solver(scheme,init;
    F = F_etd,
    steps,
    discard,
    h,
    gap)

t       = range(0,P;length = (steps - 1) ÷ gap + 1)
plot(t,x[:,:]',"--")
plot(t,[exp.(-t) exp.(-.5*t)])

F(init)

## # Lorenz63

scheme = mos.scheme_ETDRK4

σ = 10
β = 8/3
ρ = 28

L = [ -σ; -1.0; -β]
F_etd = [L, x->[σ*x[2]; ρ*x[1] - x[1]*x[3]; x[1]*x[2]]]

F =  x->[σ*(x[2]-x[1]); x[1]*(ρ - x[3]) - x[2]; x[1]*x[2] - β*x[3]]
init = complex([1;1;1.0])

h       = .14
steps   = ceil(Int,50/h)
discard = 0
gap     = 1


x = mos.my_ODE_solver(scheme,init;
    F = F_etd,
    steps,
    discard,
    h,
    gap)

plot3D(x[1,:],x[2,:],x[3,:])
