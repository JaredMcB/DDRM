using DSP,FFTW

ks = include("ks.jl")


function jaredm(P :: Real = 32π,    # Period
                n :: Int64 = 64,    # Number of linearly independent fourier modes used (beyond 0)
                )

    N = 2n+1

    ## Spatial grid and initial conditions:
    x = P*(0:N-1)/N

    ## Now we set up the equations
    # dv_k/dt = (q^2_k - q^4_k)*v_k - i*q_k/2*(convolution) for k = -n:n
    q = 2π/P*[0:n; -n:-1]

    # Nonlinear part
    ℓ = -0.5im*q
    pad = ceil(Int,n)

    function (v)
        ## if input is too short, pad it
        if length(v) < n
            v = vcat(v,zeros(n-length(v)))
        end
        ## input is a vector of Fourier coefficents v[k] for
        ## k>0; we assume v[-k] = conj(v[k]) and v[0]=0.
        v = [0,v...,conj(reverse(v))...]
        v_pad = [v[1:n+1]; zeros(pad);v[n+2:N]]
        nv = fft(real(bfft(v_pad)).^2)/length(v_pad)
        return ℓ .* [nv[1:n+1]; nv[end-n+1:end]]
    end
end

function klin(P :: Real = 32π,    # Period
              n :: Int64 = 64,    # Number of linearly independent fourier modes used (beyond 0)
              )

    ks! = ks.make_ks_field(n; L=P,alpha=0,beta=0)

    function(v)
        ## input is a vector of Fourier coefficents v[k] for
        ## k>0; we assume v[-k] = conj(v[k]) and v[0]=0.

        ## if input is too short, pad it
        if length(v) < n
            v = vcat(v,zeros(n-length(v)))
        end
        v = float(complex(v))
        F = zero(v)
        ks!(v,F)
        return [0,F...,conj(reverse(F))...]
    end
end

# julia> v=round.(8*randn(5))
# 5-element Array{Float64,1}:
#   2.0
#   5.0
#  -2.0
#   6.0
#   5.0

# julia> V=[reverse(v)...,0,v...]
# 11-element Array{Float64,1}:
#   5.0
#   6.0
#  -2.0
#   5.0
#   2.0
#   0.0
#   2.0
#   5.0
#  -2.0
#   6.0
#   5.0

# julia> truth=(-0.5*im*[0:10; -10:-1] .* ifftshift(conv(V,V)))[1:9]
# 9-element Array{Complex{Float64},1}:
#  -0.0 - 0.0im
#  -0.0 - 17.999999999999993im
#  -0.0 - 36.0im
#  -0.0 - 141.0im
#  -0.0 - 74.00000000000001im
#  -0.0 - 10.0im
#  -0.0 - 252.0im
#  -0.0 - 90.99999999999996im
#  -0.0 - 64.00000000000004im

# julia> jaredm(2pi,8)(v)[1:9]
# 9-element Array{Complex{Float64},1}:
#                     0.0 - 0.0im
#   8.526512829121202e-16 - 18.000000000000004im
#   5.684341886080802e-16 - 36.00000000000001im
#   7.673861546209083e-15 - 141.0im
#                     0.0 - 74.00000000000001im
#  2.1624521566426516e-14 - 9.999999999999984im
#  -1.663258776038883e-14 - 252.0im
#   8.849968619334368e-15 - 91.0im
#   4.043140489317773e-15 - 64.00000000000003im

# julia> klin(2pi,8)(v)[1:9]
# 9-element Array{Complex{Float64},1}:
#                      0.0 + 0.0im
#    5.551115123125783e-16 - 18.000000000000007im
#   1.7763568394002505e-15 - 35.99999999999999im
#   -4.923684827751382e-15 - 141.0im
#    6.366171895982347e-15 - 74.0im
#   -2.733532100385408e-14 - 10.0im
#  -1.2224603334267127e-14 - 251.99999999999994im
#   -7.646126771277965e-15 - 91.00000000000003im
#    8.329609062070828e-15 - 64.00000000000003im


## Set up

P = 32π
n = 64
g = x -> cos(x/16)*(1 + sin(x/16))
N = 2n+1

# Spatial grid and initial conditions:
x = P*(0:N-1)/N
u = g.(x)
v = fft(u)/N        # The division by N is to effect the DFT I want.

# Now we set up the equations
q = 2π/P*[0:n; -n:-1]
ℓ = -0.5im*q

Diff(x,y) = maximum(abs.(x - y))

## Dr. Lins function
fillout(v::Array{Complex{Float64},1}) = [0; v; reverse(conj(v))]

NonLin! = ks.make_ks_field(n; alpha = 0, beta = 0, L = P)

NonLin = function (v)
    u = zeros(Complex{Float64},n)
    NonLin!(v[2:n+1],u)
    fillout(u)
end

## My function
pad = n
K = 2n+1 + pad
NonLinNA = function (v)
    v_pad = [0;v[2:n+1]; zeros(pad);v[n+2:N]]
    nv = fft(bfft(v_pad).^2)/K
    Nv_dealiased = ℓ .* [nv[1:n+1]; nv[end-n+1:end]]
end

# N(v) truth
A       = zeros(ComplexF64,n)
a       = [-.25im, .25, -.25im,-.0625]
A[1:4]  = a
A       = fillout(A)
A       .*= ℓ

# julia> Diff(A,klin()(v[2:65]))
# 1.6553884238885772e-16
#
# julia> Diff(A,jaredm()(v[2:65]))
# 1.649639905177841e-16
#
# julia> Diff(jaredm()(v[2:65]),klin()(v[2:65]))
# 1.7542365736134082e-16
#
# julia> Diff(klin()(vv[2:65]),NonLin(vv))
# 0.0
#
# julia> Diff(jaredm()(vv[2:65]),NonLinNA(vv))
# 2.282727361246348e-16
#
# julia> Diff(jaredm()(vv[2:65]),klin()(vv[2:65]))
# 3.2191158424353993e-16
