module Model_KSE

using FFTW, Statistics

function my_KSE_solver(
    T :: Real = 150; # Length (in seconds) of time of run
    P :: Real = 32π, # Period
    N :: Int64 = 128, # Number of fourier modes used
    h :: Real = 1/4, # Timestep
    g = x -> cos(π*x/16)*(1 + sin.(π*x/16)), # Initial condition function
    T_disc = T/2,
    n_gap = 100 # 1 +  No. of EDTRK4 steps between reported data
    )


    ## Spatial grid and initial conditions:
    x = P*(1:N)/N
    u = g.(x)
    v = ifft(u)                  # ifft to get division by N this definition is better here

    ## Precompute various ETDRK4 scalar quantities:
    q = 2π/P*[0:N÷2-1; 0; N÷2-N+1:-1]
    L = q.^2 - q.^4
    E = exp.(h*L); E2 = exp.(h/2*L)

    M = 32 # no. of pts use in contour integration
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

    v_pad = [v; zeros(N)]
    F = plan_fft(v_pad)        # julia's ifft is my fft for this problem.
    iF = plan_ifft(v_pad)        # julia's fft is my ifft for this problem.

    function NonLin(v)
        v_pad = [v; zeros(N)]
        nv = F*(real(iF*v_pad)).^2
        nv[1:N]
    end

    # ## Not correcting for aliasing
    # F = plan_ifft(v)          # julia's ifft is my fft for this problem.
    # iF = plan_fft(v)          # julia's fft is my ifft for this problem.
    # NonLin(v) = F*(real(iF*v)).^2

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
            u = real.(fft(v))   # julia's fft is my ifft for this problem.
            uu[:,ni] = u
            vv[:,ni] = v
            tt[ni] = t
        end
    end
    # Energy spectrum
    EE = log.(abs.(vv).^2)

    start = n_disc+1
    uu[:,start:end], vv[:,start:end], tt[1:end-start+1]
end

end # module
