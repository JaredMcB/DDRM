"""

file: modgen_linearSDE.jl
Author: Jared McBride (08/17/2020)

This file produces times series data for a linear SDE. The form of the model
is
    dX_t = AX_t dt + σdW_t
the parameters for the equation are
    A  = dxd-matrix
    σ  = dxnu-matrix
    Xo = d-vector the initial condition

the parameters for the approximated solution are
    t_start = the start time
    t_stop  = the end time
    t_disc  = the times that are discarded
    h       = the time step
    gap     = the observation frequency
    scheme  = the numerical solver scheme

"""

# function stepper(scheme, ...)

using LinearAlgebra

function modgen_LSDE(t_start,t_stop,h;
    A, σ, Xo, t_disc = 0, gap = 1)

    d, nu = size(σ)

    size(A,1) == size(A,2) || throw(DimensionMismatch("A must be square"))
    d == size(A,1) || throw(DimensionMismatch("A and σ are noncompatable"))

    steps_tot = floor(Int,(t_stop - t_start)/h) + 1
    steps_disc = floor(Int,(t_disc - t_start)/h) + 1
    steps = steps_tot - steps_disc

    num_obs = floor(Int, steps/gap) + 1

    X = zeros(d,num_obs)
    x = Xo
    for i = 2:steps_tot
        x = (I + h*A)*x + sqrt(h)*σ*randn(nu,1)
        if i >= steps_disc && (i - steps_disc) % gap == 0
            X[:,(i - steps_disc) ÷ gap + 1] = x
        end
    end
    X
end


## Testing it
