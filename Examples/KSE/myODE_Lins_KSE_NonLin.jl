using FFTW
using PyPlot
using DSP: conv

at    = include("../../Tools/AnalysisToolbox.jl")
kse   = include("myKSE_solver.jl")
ks    = include("../../klin/ks.jl")

## Kassam and Trefethen paramaters

T       = 150
P       = 32π
n       = 64
h       = .25
g       = x -> cos(x/16)*(1 + sin(x/16))
T_disc  = 0
n_gap   = 6
aliasing= false

using Statistics: mean, var
using LinearAlgebra: diagm
mos = include("myODE_solver.jl")

N = 2n+1

## Spatial grid and initial conditions:
x = P*(0:N-1)/N
u = g.(x)
v = fft(u)/N        # The division by N is to effect the DFT I want.

## Now we set up the equations
# dv_k/dt = (q^2_k - q^4_k)*v_k - i*q_k/2*(convolution) for k = -n:n
q = 2π/P*[0:n; -n:-1]

# (diagonal) Linear part as vector
L = q.^2 - q.^4

# Now we get NonLin

NonLin(v) = -0.5im*q .* ifftshift(conv(fftshift(v),fftshift(v))[N-n:N+n])
NonLino = v -> -0.5im*q .* (fft(real(ifft(v))).^2)*N

# NonLin! = ks.make_ks_field(n; alpha = 0, beta = 0, L = P)
#
# NonLin!(v[2:n+1],v[2:n+1])
# v
#
# NonLin = function (v)
#     NonLin!(v,v)
#     [0; v ; reverse(conj(v),dims = 1)
# end

## Now we Now we use the solver
scheme = mos.scheme_ETDRK4

F = x -> diagm(L)*x + NonLin(x)     #Solely for using this with FE or RK4 solvers

F(v)

F_etd = [L, NonLino]

steps   = ceil(Int,T/h)
discard = ceil(Int,T_disc/h)

vv = mos.my_ODE_solver(scheme,v;
    F = F_etd,
    steps,
    discard,
    h,
    gap = n_gap)

findfirst(isnan, sum(vv,dims = 1))

S = size(vv,2)
uu = real(N*ifft(vv,1))
figure()
H1 = imshow(reverse(reverse(uu[:,:]',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")

EE = abs2.(vv)

i=0
i += 1
semilogy(mean(EE[:,200:250],dims=2)*2N^2,label = "me (Julia) $i")
title("Energy spectrum with aliasing")
legend()
