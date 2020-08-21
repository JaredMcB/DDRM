using LinearAlgebra
using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()

Xo = 1.5
σ = 1

t_start = 0
t_stop = 100

dt = 2^-10
Δt = 2^-5

scale = Int( Δt/dt )

Time = range(t_start,t_stop, step = dt)
N_grid = range(t_start,t_stop, step = Δt)
T = length(Time)
N = length(N_grid)

dW = sqrt(dt)*randn(T)
W = cumsum(dW)

## 'Truth'
X = zeros(T)
X[1] = Xo
for i = 1:T-1
    X[i+1] = X[i] - X[i]*(X[i]^2 - 1)*dt + σ*dW[i]
end


## 'Approximation'

Y = zeros(N)
Y[1] = Xo
for i = 1:N-1
    ΔW = W[scale*i] - W[scale*(i - 1) + 1]
    Y[i+1] = Y[i] - Y[i]*(Y[i]^2 - 1)*Δt + σ*ΔW
end


## Explicit Midpoint
F(x) = -x*(x^2 - 1)

Y_em = zeros(N)
Y_em[1] = Xo
for i = 1:N-1
    ΔW = W[scale*i] - W[scale*(i - 1) + 1]
    k1 = Y_em[i] + Δt/2*F(Y_em[i])
    Y_em[i+1] = Y_em[i] + F(k1)*Δt + σ*ΔW
end


## Explicit Tranaziod
Y_et = zeros(N)
Y_et[1] = Xo
for i = 1:N-1
    ΔW = W[scale*i] - W[scale*(i - 1) + 1]
    k0 = F(Y_et[i])
    k1 = Y_et[i] + Δt*k0
    Y_et[i+1] = Y_et[i] + 1/2*F(k1)*Δt + 1/2*k0*Δt + σ*ΔW
end


## RK4
Y_RK4 = zeros(N)
Y_RK4[1] = Xo
for i = 1:N-1
    ΔW = W[scale*i] - W[scale*(i - 1) + 1]
    k1 = Y_RK4[i]
    k2 = Y_RK4[i] + Δt/2*F(k1)
    k3 = Y_RK4[i] + Δt/2*F(k2)
    k4 = Y_RK4[i] + Δt*F(k3)
    Y_RK4[i+1] = Y_RK4[i] +
                 Δt/6*( F(k1) + 2*F(k2) + 2*F(k3) + F(k4) ) +
                 σ*ΔW
end



plot(Time,X, label = "truth")
plot!(N_grid,[Y Y_em Y_et Y_RK4],
    line = (2,[:dash :dash :dot :dot]),
    label = ["FE" "EM" "ET" "KR4"])


x = -3:.01:3

plot(x,(x.^2 .- 1).^2)
