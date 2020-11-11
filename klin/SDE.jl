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
                 sigma=0.0,
                 Stepper = EulerStepper,
                 noise=[],
                 flags...)

    x       = zeros(nsteps)
    t       = (0:(nsteps-1))*dt
    onestep = Stepper(drift, dt; sigma=sigma)

    if length(noise)>0
        for n=2:nsteps
            x[n] = onestep(x[n-1],noise[n-1])
        end
    else
        for n=2:nsteps
            x[n] = onestep(x[n-1])
        end
    end
    return x
end

end#module
