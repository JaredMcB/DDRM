using Plots
pyplot()
using StatsBase

Xo = 2.0; sigma = .2; T = 10

dt = 2^(-10)
deltat = 2^(-4)

S = Int(1/dt)*T
N = Int(1/deltat)*T

Time = T*(0:S-1)/(S-1)
Ngrid = T*(0:N-1)/(N-1)

# Wiener process
dW = [0; randn(S-1)*sqrt(dt)]
W = cumsum(dW)

# Exact
X = zeros(S);
scale = Int(deltat/dt)

X[1] = Xo;
for n = 1:S-1
    deltaW = (W[n+1] - W[n])
    X[n+1] = X[n] + dt*(-X[n]*(X[n]^2-1)) + sigma*deltaW
end

# Forward Euler
Y1 = zeros(N);

Y1[1] = Xo;
for n = 1:N-1
    deltaW = (W[scale*(n)+1] - W[scale*(n-1)+1])
    Y1[n+1] = Y1[n] + deltat*(-Y1[n]*(Y1[n]^2-1)) + sigma*deltaW
end

# 3rd-order Taylor in Deterministic part.
Y2 = zeros(N);
scale = Int(deltat/dt)

Y2[1] = Xo;
for n = 1:N-1
    deltaW = (W[scale*(n)+1] - W[scale*(n-1)+1])
    Y2[n+1] = Y2[n] + deltat*(-Y2[n]*(Y2[n]^2-1)) +
            1/2*deltat^2*(-3*Y2[n]^2 + 1) +
            1/6*deltat^3*(-6*Y2[n]) +
            sigma*deltaW
end



# Runge-Kutta 4 in Deterministic part.

f = x -> -x.*(x.^2 .- 1)

Y3 = zeros(N);
scale = Int(deltat/dt)

Y3[1] = Xo;
for n = 1:N-1
    deltaW = (W[scale*(n)+1] - W[scale*(n-1)+1])

    k1 =  f(Y3[n])
    k2 =  f(Y3[n] + deltat*k1/2)
    k3 =  f(Y3[n] + deltat*k2/2)
    k4 =  f(Y3[n] + deltat*k3)

    Y3[n+1] = Y3[n] + 1/6*deltat*(k1 + 2*k2 + 2*k3 + k4) + sigma*deltaW
end

cut = 2^0

view_t = 1 : Int(S/cut)
view_N = 1 : Int(N/cut)

plot(Time[view_t],X[view_t],
    color = :red,
    line = (1, :solid),
    label = "exact")
plot!(Ngrid[view_N],[Y1[view_N] Y2[view_N] Y3[view_N]],
    line = (2,[:dash :dashdot :dot]),
    label = ["Foward Euler" "Taylor" "RK4"])
