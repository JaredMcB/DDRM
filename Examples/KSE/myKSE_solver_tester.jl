using FFTW
using PyPlot

using Statistics: mean # for energy spectrum

at    = include("../../Tools/AnalysisToolbox.jl")
mykse = include("myKSE_solver.jl")
ks = include("../../klin/ks.jl")
kse   =  include("Model_kse.jl")

## Kassam and Trefethen paramaters
T       = 1500
P       = 32π
n       = 64
h       = 1/4
g       = x -> cos(x/16)*(1 + sin(x/16))
T_disc  = 0
n_gap   = 6

# My New Solver
vv = mykse.my_KSE_solver(T; P, n, h, g, T_disc, n_gap)

uu = real((2n+1)*ifft(vv,1))
figure()
H1 = imshow(reverse(reverse(uu[:,:]',dims=2),dims = 1), extent=[0,P,0,150], aspect="auto")

# Dr. Lin's
steps = ceil(Int,T/h/n_gap) - 1

vv_k = ks.run_ks(ks.kt_init(;n=129),steps,1.5;
             nsubsteps = 6,
             verbose = false, L = P)

uu_k = real((2n+1)*ifft([zeros(1, steps+1); vv_k; reverse(conj(vv_k),dims = 1)],1))
H1 = imshow(reverse(reverse(uu_k[:,1:100]',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")

## Plot trajectories
function test(vv,vv_k)
    figure()
    modes = 1:5
    for i in modes
        subplot(length(modes),1,i)
        plot(real(vv[i+1,:]),label ="mine")
        plot(real(vv_k[i,:]),":",label ="klin")
        title("mode $i")
        legend()
    end

    # Autocorellations
    lags = -100:100
    figure()
    for i in modes
        subplot(length(modes),1,i)
        A = at.my_autocor(real(vv[i+1,:]),lags)
        plot(1.5*lags,A,label ="mine")
        A_k = at.my_autocor(real(vv_k[i,:]),lags)
        plot(1.5*lags,A_k,":",label ="klin")
        title("mode $i")
        legend()
    end

    # Energy spectrum
    figure()
    E = mean(abs2.(vv[2:65,:]),dims = 2)[:]
    semilogy(E,label ="mine")
    E_k = mean(abs2.(vv_k),dims = 2)[:]
    semilogy(E_k,label ="klin")
    legend()
end

## 2017 parameters
T        = 10^4     # Length (in seconds) of time of run
T_disc   = T÷ 2     # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
n        = 96       # Number of fourier modes used
h        = .001     # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
n_gap    = 100


# My New Solver
@time vv = mykse.my_KSE_solver(T; P, n, h, g, T_disc, n_gap)

uu = real((2n+1)*ifft(vv,1))
H1 = imshow(reverse(reverse(uu[:,end-1000+1 : end]',dims=2),dims = 1), extent=[0,P,0,150], aspect="auto")

# Dr. Lin's
steps   = ceil(Int,T/h/n_gap) - 1
discard = ceil(Int,T_disc/h/n_gap)

steps += discard

@time vv_k = ks.run_ks(ks.kt_init(;n = 193),steps,.1;
             nsubsteps = 100,
             verbose = false, L = P)

vv_k = [zeros(1, steps+1); vv_k; reverse(conj(vv_k),dims = 1)]
uu_k = real((2n+1)*ifft(vv_k,1))
figure()
H1 = imshow(reverse(reverse(uu_k[:,end-1000+1 : end]',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")

## Plot trajectories
test(vv[:,end-50000+1 : end], vv_k[:,end-50000+1 : end])

## Energy Spectrum

V = [vv[:,end - 50000 +1 : end],vv_k[:,end - 50000 +1 : end]]

using Statistics: mean
for v in V
    E = mean(abs2.(v)[1:30],dims = 2)
    plot(E)
end














# My New Solver
vv = mykse.my_KSE_solver(T;
                        P,
                        n,
                        h,
                        g,
                        T_disc,
                        n_gap,
                        aliasing)``

uu = real((2n+1)*ifft(vv,1))
ender = findfirst(isnan, sum(uu,dims = 1))[2]-1000
H1 = imshow(reverse(uu[:,1:ender]',dims = 1), extent=[0,P,0,1000], aspect="auto")

EE = abs2.(vv)
figure()
semilogy(mean(EE[:,1:ender],dims=2),label = "me (Julia) $i")

lags = -50:50
A = at.my_autocor(vv[3,1:ender],lags)
plot(A)
# MyOld Solver
uu_o, vv_o, tt = kse.my_KSE_solver(T;
                        P,
                        N = 2n+1,
                        h,
                        g,
                        T_disc,
                        n_gap,
                        aliasing)

ender = findfirst(isnan, sum(uu_o,dims = 1))[2]-10
figure()
H1 = imshow(reverse(uu_o[:,1:ender]',dims = 1), extent=[0,P,0,150], aspect="auto")

EE_o = abs2.(vv_o)
semilogy(mean(EE_o[:,1:ender],dims=2),label = "me (Julia) $i")

A_o = at.my_autocor(vv_o[3,1:ender],lags)
plot(A_o)
