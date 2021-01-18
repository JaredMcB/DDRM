"""
file: New_KSE_solver.jl
author: Jared McBride (1-15-2021; American Fork, Utah)

Discription: This is a module containing a KSE solver that works better than the last one.

"""

module myODE_solver

using Statistics: mean

Complex128 = Complex{Float64}

function my_ODE_solver(scheme,
                       init;
                       F,
                       steps,           # after discard
                       discard,         # this is in steps
                       h,
                       gap,             # so output is (steps - 1) ÷ gap + 1 long
                       flags...)

        n = size(init,1) # Dimension of system
        x = zeros(Complex128,n,(steps - 1) ÷ gap + 1)

        # if we use a more gerenal one step method
        step! = scheme(F, h)

        # main stepping loop
        temp = copy(init)
        for n = 1:steps+discard
                if (n > discard) & ((n - discard - 1) % gap == 0)   # save state
                        x[:,(n-discard-1)÷gap+1] = temp
                end
                temp = step!(temp,temp)                               # advance state
        end
        x
end

function scheme_RK4(F, h) # Assume autonomus RHS
        step! = function (u,v)
                k1 = F(u)
                k2 = F(u+h*k1/2)
                k3 = F(u+h*k2/2)
                k4 = F(u+h*k3)

                v[:] = u + h/6*(k1 + 2k2 + 2k3 + k4)
        end
end


function scheme_ETDRK4(F,       # packet with diagonal linear part and nonlinear part seperated
                       h)       # Assume autonomus RHS

       L       = F[1]           # (linear part) Assumed to be diagonal here, just Col vector
       NonLin  = F[2]           # Nonlinear part

       N = size(L,1)
       E = exp.(h*L); E2 = exp.(h/2*L)

       M = 16                          # no. of pts use in contour integration
       r = exp.(im*π*((1:M) .-.5)/M)   # roots of unit suggested by Kassam and Trefethen
       LR = h*L*ones(M)' + ones(N)*r'  # the second dim varies r the first vaeries L

       Q = h*real(mean((exp.(LR/2) .- 1)./LR, dims=2))[:]
       f1 = h*real(mean((-4 .- LR+exp.(LR).*(4 .- 3*LR + LR.^2))./LR.^3,dims=2))[:]
       f2 = h*real(mean((2 .+ LR+exp.(LR).*(-2 .+ LR))./LR.^3,dims=2))[:]
       f3 = h*real(mean((-4 .- 3*LR-LR.^2+exp.(LR).*(4 .- LR))./LR.^3,dims=2))[:]

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
               @. v[:] =  E*u + Nu*f1 + 2*(Na+Nb)*f2 + Nc*f3
               v
       end
end

function scheme_FE(F,h)
        step! = function (u,v)
                v[:] = u + h*(F(u))
        end
end




end #module
