```
This version was committed on Aug 20.

With the earlier neive from of dealiasing the code runs the full 2017
run with no Nan's. It also runs the Trefethen data (extended by ten times)
with no nan's
```
module Model_KSE

using FFTW, Statistics

function my_KSE_solver(
    T :: Real = 150; # Length (in seconds) of time of run
    P :: Real = 32π, # Period
    n :: Int64 = 64, # Number of fourier modes used
    h :: Real = 1/4, # Timestep
    g = x -> cos(x/16)*(1 + sin.(x/16)), # Initial condition function
    T_disc = T/2,
    n_gap = 100, # 1 +  No. of EDTRK4 steps between reported data
    aliasing = false
    )

    N = 2n

    ## Spatial grid and initial conditions:
    x = P*(1:N)/N
    u = g.(x)
    v = fft(u)/N            # The division by N is to effect the DFT I want.

    ## Precompute various ETDRK4 scalar quantities:
    q = 2π/P*[0:n-1;0; -n+1:-1]
    L = q.^2 - q.^4
    E = exp.(h*L); E2 = exp.(h/2*L)

    M = 16 # no. of pts use in contour integration
    r = exp.(im*π*((1:M) .-.5)/M) # roots of unit suggested by Kassam and Trefethen
    LR = h*L*ones(M)' + ones(N)*r' # the second dim varies r the first vaeries L

    Q = h*real(mean((exp.(LR/2) .- 1)./LR, dims=2))[:]
    f1 = h*real(mean((-4 .- LR+exp.(LR).*(4 .- 3*LR + LR.^2))./LR.^3,dims=2))[:]
    f2 = h*real(mean((2 .+ LR+exp.(LR).*(-2 .+ LR))./LR.^3,dims=2))[:]
    f3 = h*real(mean((-4 .- 3*LR-LR.^2+exp.(LR).*(4 .- LR))./LR.^3,dims=2))[:]

    ## Some declareations

    a = Complex.(zeros(N))
    b = Complex.(zeros(N))
    c = Complex.(zeros(N))
    Nv = Complex.(zeros(N))
    Na = Complex.(zeros(N))
    Nb = Complex.(zeros(N))
    Nc = Complex.(zeros(N))

    # Main time-stepping loop
    n_max = round(Int,T/h)
    n_obs = floor(Int,n_max/n_gap)
    n_disc = floor(Int,T_disc/h/n_gap)
    ℓ = -0.5im*q

    padding = aliasing ? 0 : N

    v_pad = [v[1:n]; zeros(padding);v[n+1:N]]
    F = plan_fft(v_pad)
    iF = plan_ifft(v_pad)

    if aliasing
        NonLin = v -> F*(real(iF*v)).^2*N
    else
        NonLin = function (v)
            v_pad = [v[1:n]; zeros(padding);v[n+1:N]]
            nv = F*(real(iF*v_pad)).^2*4N
            [nv[1:n]; nv[N+n+1:2N]]
        end
    end

    vv = complex(zeros(N, n_obs+1)); vv[:,1]= v
    uu = zeros(N, n_obs+1); uu[:,1]= u
    tt = zeros(n_obs+1); tt[1] = 0
    for n = 1:n_max
        t = n*h
        Nv .= ℓ.* NonLin(v)
        @.  a  =  E2*v + Q*Nv
        Na .= ℓ.* NonLin(a)
        @. b  =  E2*v + Q*Na
        Nb .= ℓ.* NonLin(b)
        @. c  =  E2*a + Q*(2Nb-Nv)
        Nc .= ℓ.* NonLin(c)
        @. v =  E*v + Nv*f1 + 2*(Na+Nb)*f2 + Nc*f3
        if n % n_gap == 0
            ni = Int64(n÷n_gap) + 1
            u = real.(ifft(v))*N
            uu[:,ni] = u
            vv[:,ni] = v
            tt[ni] = t
        end
    end

    start = n_disc+1
    uu[:,start:end], vv[:,start:end], tt[1:end-start+1]
end
end #module
