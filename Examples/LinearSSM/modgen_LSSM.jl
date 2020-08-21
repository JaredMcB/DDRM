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


function modgen_LSSM(t_start,t_stop,h;
    F, G = I, H = I, R = I, Q, Xo, t_disc = 0, gap = 1)

    dx = size(F,1)
    dx == size(F,2) || throw(DimensionMismatch("F must be square"))

    G == I ? nu = dx : nu = size(G,2)
    H == I ? dy = dx : dy = size(H,1)

    size(F,1) == size(F,2) || throw(DimensionMismatch("F must be square"))

    steps_tot = floor(Int,(t_stop - t_start)/h) + 1
    steps_disc = floor(Int,(t_disc - t_start)/h) + 1
    steps = steps_tot - steps_disc

    u = sqrt(R + zeros(nu,nu)) * randn(nu,steps_tot)
    v = sqrt(Q + zeros(dy,dy)) * randn(dy,steps_tot)
    X = zeros(dx,steps_tot); X[:,1] = Xo
    Y = zeros(dy,steps_tot)
    for i = 2:steps_tot
        X[:,i] = F*X[:,i-1] + G*u[:,i]
        Y[:,i-1] = H*X[:,i-1] +  v[:,i-1]
    end
    X = X[:,steps_disc+1:gap:end]
    Y = Y[:,steps_disc+1:gap:end]
end
