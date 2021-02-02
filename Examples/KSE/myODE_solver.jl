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
        for n = 1:steps
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


function scheme_FE(F,h)
        step! = function (u,v)
                v[:] = u + h*(F(u))
        end
end




end #module
