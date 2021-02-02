"""
This version was committed on Aug 20.

With the earlier neive from of dealiasing the code runs the full 2017
run with no Nan's. It also runs the Trefethen data (extended by ten times)
with no nan's

I worked on it a while and could not get the code to dealias. So, I started over
    with a different a approach. The result is the stand-alone )DE solver
    "myOde_solver" with ETDrk4 capabilities. I use that model in the smodel
    myKSE_solver.jl. This was on Jan 18, 2021.
"""
module Model_KSE

using FFTW, Statistics

function my_KSE_solver(
    T :: Real = 150;    # Length (in seconds) of time of run
    P :: Real = 32π,    # Period
    N :: Int64 = 128,   # Number of fourier modes used
    h :: Real = 1/4,    # Timestep
    g = x -> cos(x/16)*(1 + sin(x/16)),     # Initial condition function
    T_disc = T/2,
    n_gap = 100,        # 1 +  No. of EDTRK4 steps between reported data
    aliasing = false
    )


    nf = floor(Int, N/2)        # Notice if nc+nf+1 = N
    nc = ceil(Int, N/2) - 1     # The idea is for these to give the
                                # Number of positve (nc) and negative
                                # (nf) modes used in approximation
                                # so, in keeping with the convention
                                # in FFTW if N is odd we simply have
                                # nf = nc = (N-1)/2
                                # However, if N even nf = nc + 1

    ## Spatial grid and initial conditions:
    x = P*(0:N-1)/N
    u = g.(x)
    v = fft(u)/N                # The division by N is to effect the DFT I want.
    isodd(N) || (v[nc+2] = 0)   # We want set this equal to zero since the signal
                                # is real

    ## Precompute various ETDRK4 scalar quantities:
    q = 2π/P*[0:nc; -nf:-1]
    isodd(N) || (q[nc+2] = 0)   # Here we definitly need the odd-man-out (that extra
                                # negitive mode when N is even) to be zeroed out.

    L = q.^2 - q.^4
    E = exp.(h*L); E2 = exp.(h/2*L)

    M = 16                          # no. of pts use in contour integration
    r = exp.(im*π*((1:M) .-.5)/M)   # roots of unit suggested by Kassam and Trefethen
    LR = h*L*ones(M)' + ones(N)*r'  # the second dim varies r the first vaeries L

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

    pad = aliasing ? 0 : ceil(Int,3N/2)

    v_pad = [v[1:nc+1]; zeros(pad); v[nc+2:N]]
    K = size(v_pad,1)
    F = plan_fft(v_pad)
    iF = plan_ifft(v_pad)

    if aliasing
        NonLin = v -> F*(real(iF*v)).^2*N
    else
        NonLin = function (v)
            v_pad = [v[1:nc+1]; zeros(pad);v[nc+2:N]]
            nv = F*(real(iF*(v_pad)).^2)*K
            Nv_dealiased = [nv[1:nc+1]; nv[end-nf+1:end]]
            # ifftshift(conv(fftshift(v),fftshift(v))[N-(n-1):N+n])/N
            # v_pad = [v[1:n]; zeros(pad);v[n+1:N]]
            # nv = F*(real(iF*v_pad)).^2*K/N
            # [nv[1:n]; nv[end-n+1:end]]
        end
    end

    vv = complex(zeros(N, n_obs+1));    vv[:,1]= v
    uu = zeros(N, n_obs+1);             uu[:,1]= u
    tt = zeros(n_obs+1);                tt[1] = 0

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
