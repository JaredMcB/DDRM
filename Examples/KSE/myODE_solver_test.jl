"""
Here I want to test my ODE solver.

The first test will be the ODE wll be the second order
d²x/dt² = -x , x(0) = 0, x'(0) = 1
This sound give the sine function.
"""

mos = include("myODE_solver.jl")

scheme = mos.scheme_FE
F(x) = [0 1;-1 0]*x
init = complex([0;1])

h       = .01
steps   = ceil(Int,2π/h)
discard = 0
gap     = 1



x = mos.my_ODE_solver(scheme,init,F;
    steps,
    discard,
    h,
    gap)
