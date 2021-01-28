"""
Here I make my own ETDRK4 scheme

This module allows for the output of a ETDRK4 stepper.

It is designed for compatability with myODE_solver.jl

"""

module myETDRK4

using Statistics: mean

function cont_quad(f::Function;
    N,
    M = 64,
    r = 2)
    gam = r*exp.(im*2pi*( (1:M) .- .5 )/M)
    function (hL)
        Gam = hL*ones(M)' + ones(N)*gam'
        mean(f.(Gam),dims = 2)[:]
    end
end

##############################################################


function scheme_ETDRK4(F,       # packet with diagonal linear part and nonlinear part seperated
                       h)       # Assume autonomus RHS

    L       = F[1]           # (linear part) Assumed to be diagonal here, just Col vector
    NonLin  = F[2]           # Nonlinear part

    N = size(L,1)
    E = exp.(h*L)
    E2 = exp.(h/2*L)

    f_Q(z)     = (exp(z/2) - 1)/z
    f_alpha(z) = (-4 - z + exp(z)*(4 - 3z + z^2))/z^3
    f_beta(z)  = (2 + z + exp(z)*(-2 + z))/z^3
    f_gamma(z) = (-4 - 3z - z^2 + exp(z)*(4 - z))/z^3

    F_Q = cont_quad(f_Q;N)
    F_alpha = cont_quad(f_alpha;N)
    F_beta = cont_quad(f_beta;N)
    F_gamma = cont_quad(f_gamma;N)

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
           @. v[:] =  E*u + alpha*Nu + 2beta*(Na+Nb) + gamma*Nc
           v
    end
end

end # module
