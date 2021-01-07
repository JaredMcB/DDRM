############################################################################
## ks.jl

## This merges two modules.  The eventual goal is to provide

## 1. a self-contained solver for the Kuramoto-Sivashinsky
##    equation (KSE) on the circle; and

## 2. support for fitting model to data, making much the same
##    choices as Fei did for the Physica D paper (so as to
##    allow direct comparisons).

## First priority is #2.

module ks

Complex128 = Complex{Float64}

include("util.jl")      #include("WienerROM/Util/util.jl")
include("fftwutil.jl")  #include("WienerROM/Extra/fftwutil.jl")
include("etdrk.jl")     #include("WienerROM/Extra/etdrk.jl")

using FFTW,.fftwutil,.etdrk


######################################################
## Vector field

## Note: the dispersion relation is parametrized, so the
## same function can actually be used to solve e.g. the
## viscous Burgers equation.

export   make_ks_field
function make_ks_field(N; alpha=1.0, beta=1.0, L=2*pi,
                       dispersion = v->alpha*v^2-beta*v^4)

    ## parameters, temporary storage, etc
    N3 = 3*N+1  # for dealiasing
    C  = 2*pi/L

    Freal    = zeros(N3)
    do_hc2r! = make_r2r!(Freal, FFTW.HC2R)
    do_r2hc! = make_r2r!(Freal, FFTW.R2HC)

    ## spectrum of linear part
    spec = [ dispersion(C*k) for k=1:N ]

    ## function to evaluate nonlinear term

    function ks_field!(U,F)

        Freal[1] = 0.0

        let i  = 2,
            ii = N3

            for k=1:N
                Freal[i]  = U[k].re
                Freal[ii] = U[k].im
                i = i+1
                ii = ii-1
            end

            while i <= ii
                Freal[i] = Freal[ii] = 0.0
                i = i+1
                ii = ii-1
            end
        end

        do_hc2r!(Freal)
        for i=1:N3
            Freal[i] = Freal[i]^2/N3
        end
        do_r2hc!(Freal)

        let i  = 2,
            ii = N3
            for k=1:N
                F[k] = -(0.5*C*k) * Complex(-Freal[ii],Freal[i]) + spec[k]*U[k]
                i += 1
                ii -= 1
            end
        end
    end
end

######################################################
## Solver

export   make_ks_stepper
function make_ks_stepper(N, dt; alpha=1.0, beta=1.0, L=2*pi,
                         dispersion = v->alpha*v^2-beta*v^4,
                         delay=60.0)

    step! = make_etdrk4_stepper(make_ks_field(N; alpha=0., beta=0., L=L),
                                [ dispersion(2*pi*k/L) for k=1:N ],
                                dt)
end

export   run_ks
function run_ks(nsteps, dt; flags...)
        run_ks(kt_init(), nsteps, dt; flags...)
end

function run_ks(init, nsteps, dt;
                nsave=length(init),
                nsubsteps=1,
                verbose=false,
                delay=60.0,
                flags...)

    function echo(s)
        if verbose
            println(s)
            flush(stdout)
        end
    end

    dt /= nsubsteps

    N      = length(init)
    z      = zeros(Complex128, nsave, nsteps+1)
    z[:,1] = init[1:nsave]
    step!  = make_ks_stepper(N, dt; flags...)

    U = copy(init)
    V = zeros(Complex128,N)

    function onebigstep!(n)
        step!(U,nsubsteps,V)
        z[:,n+1] = V[1:nsave]
        let tmp = U
            U = V
            V = tmp
        end
    end

    echo("#modes == $N")

    if verbose
        foreach(onebigstep!, 1:nsteps, "run_ks"; delay=delay)
    else
        for n=1:nsteps
            onebigstep!(n)
        end
    end

    return z
end

export   run_ks0
function run_ks0(nsteps, dt; flags...)
    run_ks0(kt_init(), nsteps, dt; flags...)
end

function run_ks0(init, nsteps, dt; flags...)

    N = length(init)
    z = copy(init)

    step! = make_ks_stepper(N, dt; flags...)
    step!(z,nsteps,z)

    return z
end


################################
## Initial conditions

## From [Kassam-Trefethen]
function kt_init(;n=128, real=false)
    dx = 32*pi/n
    init = [cos(x/16)*(1+sin(x/16)) for x in (0:(n-1))*dx]
    if real
        init
    else
        fourier(init)
    end
end

## Randomly (gaussian) excite some modes
function random_init(N,M=N)
    @assert M <= N
    [[(randn() + im*randn())/sqrt(M) for k=1:M]...,
     zeros(N-M)...]
end

## Fei's init
function fei_init(;n=2*96, real=false)
    dx   = 2*pi/sqrt(0.085)/n
    init = [(1+sin(x))*cos(x) for x in (0:(n-1))*dx]
    if real
        init
    else
        fourier(init)
    end
end


######################################################
## Ansatz

function aim_vec!(target, u::Vector{Complex128}, k::Int)
    K = length(u)  # num modes, Fei's notation
    for j=1:K
        # target[j] = u__(u,j+K)*u__(u,j+K-k)
        target[j] = im*u__(u,j+K)*conj(u__(u,j+K-k))
    end
end

## This implements $\tilde{u}^n_j$ from the Physica D paper:
function u__(u::Vector{Complex128}, j::Int)

    K = length(u)  # num modes, Fei's notation

    if 1 <= j <= K
        return u[j]
    elseif K < j <= 2*K
        sum=0.0
        for l=(j-K):K
            sum += u[l]*u[j-l]
        end
        return im*sum
    else
        error("j==$j when K==$K")
    end
end


################################
## RK4 is part of the ansatz

function make_rk4_stepper(f!, dt, ndim)
    make_rk4_stepper(f!, dt, Float64, ndim)
end

function make_rk4_stepper(f!, dt, thetype, ndim)
    k0   = zeros(thetype, ndim)
    k1   = zeros(thetype, ndim)
    k2   = zeros(thetype, ndim)
    k3   = zeros(thetype, ndim)
    tmp  = zeros(thetype, ndim)
    xout = zeros(thetype, ndim)

    function step(x,xx)
        f!(x,k0)
        for k=1:ndim
            tmp[k] = x[k] + 0.5*dt*k0[k]
        end
        f!(tmp,k1)
        for k=1:ndim
            tmp[k] = x[k] + 0.5*dt*k1[k]
        end
        f!(tmp,k2)
        for k=1:ndim
            tmp[k] = x[k] + dt*k2[k]
        end
        f!(tmp, k3)
        for k=1:ndim
            xx[k] = x[k] + (dt/6)*k3[k] + (dt/3)*k2[k] +(dt/3)*k1[k] + (dt/6)*k0[k]
        end
        return xx
    end

    function step(x)
        return step(x,xout)
    end
end


## Lorenz 63 is just for testing RK4
function make_lorenz63(;sigma=10.0, rho=28.0, beta=8.0/3.0)

    Xdot = zeros(3)

    function f(X,Xdot)
        x = X[1]
        y = X[2]
        z = X[3]
        Xdot[1] = sigma*(y-x)
        Xdot[2] = x*(rho-z)-y
        Xdot[3] = x*y-beta*z
    end

    function f(X)
        f(X,Xdot)
        return Xdot
    end
    # function df(X)
    #     x = X[1]
    #     y = X[2]
    #     z = X[3]
    #     return [-sigma sigma 0.0; rho-z -1.0 -x; y x -beta]
    # end
    # (f,df)
end


################################
## Package it all up

#=

Columns of predictors consist of

1. State vector :: cols 1:K
2. RK4 update   :: cols (K+1):(2K)
3. AIM terms    :: cols (2K+1):(2K+K^2)

=#

## new version
export   make_observer
function make_observer(dt, num_modes; aim=true, flags...)
    if aim
        make_observer_with_aim(dt, num_modes; flags...)
    else
        make_observer_without_aim(dt, num_modes; flags...)
    end
end

function make_observer_with_aim(dt, num_modes; flags...)

    ## Temporary storage
    K          = num_modes
    thefield   = make_ks_field(K; flags...)
    thestepper = make_rk4_stepper(thefield, dt, Complex128, K)
    bigmat     = zeros(Complex128, K, 2*K+K^2)
    znext      = zeros(Complex128, K)

    function observe(z,n)

        thestepper(z[n], znext)

        for k=1:K
            bigmat[k,k]   = z[n][k]
            bigmat[k,K+k] = (znext[k]-z[n][k])/dt
            let kmin=2*K+1,
                kmax=2*K+K
                for k=1:K
                    aim_vec!(view(bigmat,k,kmin:kmax), z[n], k)
                    kmin,kmax += K
                end
            end
        end
        return bigmat
    end
end

function make_observer_without_aim(dt, num_modes; aim=true, flags...)

    ## Temporary storage
    K          = num_modes
    thefield   = make_ks_field(K; flags...)
    thestepper = make_rk4_stepper(thefield, dt, Complex128, K)
    bigmat     = zeros(Complex128, K, 2*K)
    znext      = zeros(Complex128, K)

    function observe(z,n)

        thestepper(z[n], znext)

        for k=1:K
            bigmat[k,k]   = z[n][k]
            bigmat[k,K+k] = (znext[k]-z[n][k])/dt
        end

        return bigmat
    end
end


################################
## Abusing the language
import Base.+
function +(xy::Tuple{T,T},z::T) where T
    (xy[1]+z,xy[2]+z)
end

end#module
