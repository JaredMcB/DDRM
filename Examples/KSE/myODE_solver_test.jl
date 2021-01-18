"""
Here I want to test my ODE solver.

The first test will be the ODE wll be the second order
d²x/dt² = -x , x(0) = 0, x'(0) = 1
This sound give the sine function.
"""

using PyPlot

mos = include("myODE_solver.jl")

scheme = mos.scheme_ETDRK4
F(x) = [0 1;-1 0]*x
init = complex([0;1.0])

h       = .01
steps   = ceil(Int,2π/h)
discard = 0
gap     = 1

x = mos.my_ODE_solver(scheme,init; #F,
    steps,
    discard,
    h,
    gap)
x
t       = range(0,2π;length = steps)
plot(t,x[1,:])
plot(t,sin.(t))

### True
F_true(t) = exp([0 1;-1 0]*t)*init
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
