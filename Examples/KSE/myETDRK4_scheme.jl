"""
Here I make my own ETDRK$ scheme


"""
using Statistics: mean
using LinearAlgebra: diagm
using FFTW
mos = include("myODE_solver.jl")

T       = 150
P       = 32π
n       = 64
h       = .25
g       = x -> cos(x/16)*(1 + sin(x/16))
T_disc  = 0
n_gap   = 6

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

# Nonlinear part
ℓ = -0.5im*q

K = 3N + 1
v_pad = zeros(ComplexF64,K)
Fp = plan_fft(v_pad)
iFp = plan_bfft(v_pad)

function NonLin(v)
    v_pad[1:n+1]        = v[1:n+1]
    v_pad[end-n+1:end]  = v[n+2:end]
    nv = Fp*(real(iFp*(v_pad)).^2)/K
    return ℓ .* [nv[1:n+1]; nv[end-n+1:end]]
end


function cont_quad(f::Function;
    M = 64,
    r = 2)
    gam = r*exp.(im*2pi*( (1:M) .- .5 )/M)
    function (hL)
        Gam = hL*ones(M)' + ones(N)*gam'
        mean(f.(Gam),dims = 2)[:]
    end
end

##############################################################

# module myETDRK4
#
# function scheme_ETDRK4(F,       # packet with diagonal linear part and nonlinear part seperated
#                        h)       # Assume autonomus RHS
#
# L       = F[1]           # (linear part) Assumed to be diagonal here, just Col vector
# NonLin  = F[2]           # Nonlinear part

N = size(L,1)
E = exp.(h*L)
E2 = exp.(h/2*L)

f_Q(z)     = (exp(z/2) - 1)/z
f_alpha(z) = (-4 - z + exp(z)*(4 - 3z + z^2))/z^3
f_beta(z)  = (2 + z + exp(z)*(-2 + z))/z^3
f_gamma(z) = (-4 - 3z - z^2 + exp(z)*(4 - z))/z^3

F_Q = cont_quad(f_Q)
F_alpha = cont_quad(f_alpha)
F_beta = cont_quad(f_beta)
F_gamma = cont_quad(f_gamma)

Q     = h*F_Q(h*L)
alpha = h*F_alpha(h*L)
beta  = h*F_beta(h*L)
gamma = h*F_gamma(h*L)

a = Complex.(zeros(N)) # required because of the @. macro
b = Complex.(zeros(N))
c = Complex.(zeros(N))

function step!(u,v)
       Nu = NonLin(u)
       @.  a  =  E2*u + Q*Nu
       Na = NonLin(a)
       @. b  =  E2*u + Q*Na
       Nb = NonLin(b)
       @. c  =  E2*a + Q*(2Nb-Nu)
       Nc = NonLin(c)
       @. v[:] =  E*u + h*alpha*Nu + 2h*beta*(Na+Nb) + h*gamma*Nc
       v
end
# end # module

step!(v,v)
v
u = copy(v)

Nu = NonLin(u)
@.  a  =  E2*u + Q*Nu
Na = NonLin(a)
@. b  =  E2*u + Q*Na
Nb = NonLin(b)
@. c  =  E2*a + Q*(2Nb-Nu)
Nc = NonLin(c)
@. v[:] =  E*u + alpha*Nu + 2*beta*(Na+Nb) + gamma*Nc





















x = P*(0:N-1)/N
u = g.(x)

init = fft(u)/N

step!(init,init)


N = size(init,1) # Dimension of system

steps   = ceil(Int,T/h)
discard = ceil(Int,T_disc/h)
gap = n_gap

x = zeros(ComplexF64,N,(steps - 1) ÷ gap + 1)

# if we use a more gerenal one step method
step! = scheme(F, h)

# main stepping loop
temp = copy(init)
for i = 1:steps+discard
    global temp
        if (i > discard) & ((i - discard - 1) % gap == 0)   # save state
                x[:,(i-discard-1)÷gap+1] = temp
        end
        temp = step!(temp,temp)                               # advance state
end
x
T = step!(v,v)
V = copy(v)
T - v


using PyPlot

uu = real((2n+1)*ifft(x,1))
ender = findfirst(isnan, sum(uu,dims = 1))[2]-10
figure()
H1 = imshow(reverse(reverse(uu[:,1:end]',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")
