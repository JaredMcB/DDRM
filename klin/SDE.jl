################################
## SDE.jl

## Solving SDEs with additive noise.

module SDE

export EulerStepper, RK2Stepper, onepath

## Euler solver
function EulerStepper(f,dt; sigma=0.0)
    sigma_dt = sigma*sqrt(dt)
    function(x,w=randn())
        x + dt*f(x) + sigma_dt*w
    end
end

## RK2 solver (not more accurate, but is more stable)
function RK2Stepper(f,dt; sigma=0.0)
    sigma_dt = sigma*sqrt(dt)
    function(x,w=randn())
        k = x + 0.5*dt*f(x)
        return x + dt*f(k) + sigma_dt*w
    end
end

## Generate one sample path.
function onepath(drift, nsteps, dt;
                 init = 0.0,
                 sigma=0.0,
                 Stepper = EulerStepper,
                 noise=[],
                 ninit=0,
                 stride=1,
                 )

    nx = div(nsteps - ninit, stride)
    x  = zeros(nx)
    x0 = init

    onestep = Stepper(drift, dt; sigma=sigma)

    ## initialize
    for i=1:ninit
        x0 = onestep(x0)
    end
    x[1] = x0
    
    if length(noise)>0
        if stride > 1
            error("unsupported")
        end

        for n=2:nx
            x[n] = onestep(x[n-1],noise[n-1])
        end
    else
        for n=2:nx
            x[n] = onestep(x[n-1])
            for i=2:stride
                x[n] = onestep(x[n])
            end
        end
    end
    return x
end

end#module
