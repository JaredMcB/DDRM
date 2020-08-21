using DSP
using FFTW
using LinearAlgebra
using Plots
pyplot()

include("..\\DataGen.jl")

x = rand(100)

y = rand(100)

Fx = fft(x)
Fy =fft(y)

z = conv(x,y)

z_fft = ifft(Fx .* Fy)

plot(z)
plot!(real.(z_fft))

yy = reverse(y)

z_dot = dot(x,yy)
z[100]
z_fft = real.(z_fft)

z_fft[100]

x = rand(1000)
y = rand(1000)

z = dot(x,y)

zz =dot_fast(x,y)
N = 10000
z_dif = zeros(N)
for i = 1:N
    x = 100*randn(10^6)
    y = 100*randn(10^6)
    z_dif[i] = dot(x,y) - dot_fast(x,y)
end

d = rand()

function dot_fast(x::Vector{Number},y::Vector{Number})
    d = size(x,1)
    d == size(y,1) || throw(DimensionMismatch("x and y must have same length"))

    yy = reverse(y)
    real(ifft(fft(x) .* fft(yy))[d])
end

_autodot(x::Vector{Number}, lx::Int, l::Int) = dot_fast(view(x, 1:(lx-l)), view(x, (1+l):lx))

function autocov_fft(x::Vector{Number}, lags::UnitRange{Int})
    lx = length(x)
    m = length(lags)

    z = x .- mean(x)
    r = zeros(m)
    for k = 1 : m  # foreach lag value
        r[k] = _autodot(z, lx, lags[k]) / lx
    end
    return r
end
