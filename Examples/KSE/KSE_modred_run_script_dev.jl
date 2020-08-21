
using JLD
using DSP # For conv function in Psi
using Dates

include("Model_KSE.jl")
include("Model_Reduction.jl")

## Parameters for KSE model

T = 400 # Length (in seconds) of time of run
T_disc = 100 # Length (in seconds) of time discarded
P = 32π  # Period
N = 128  # Number of fourier modes used
h = 1e-2 # Timestep
g = x -> cos(π*x/16)*(1 + sin.(π*x/16))
q = 2π/P*(0:N-1)

## Parameters for model reduction
d = 5 # No. of lowest modes taken in reduced model
M_out = 20 # No. of coeficinets in Wiener filter output

## Get full model data ######################################################

uu, vv, tt =  my_KSE_solver(T,
         T_disc  = T_disc,
         P = P,
         N = N,
         h = h,
         g = g)

# Save this data alog with parameters

paramaters = Dict(
   "T" => T,
   "P" => P,
   "N" => N,
   "h" => h,
   "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
   "q" => q,
   "d" => d,
   "tm" => now(),
   "M_out" => M_out)
dat = Dict(
   "dat_uu" => uu,
   "dat_vv" => vv,
   "dat_tt" => tt)
Data = merge(paramaters,dat)
# save("../data/KSE_Data/KSE_sol.jld",Data)
save("Data\\KSE_sol_wienerfilter.jld",Data)
println("data saved")

# # Load Old Data
# @time Data = load("../data/KSE_Data/KSE_sol.jld")
# uu = Data["dat_uu"]
# vv = Data["dat_vv"]

## Get Reduced model #########################################################

# collect observations
obs_gap = 1
V_obs = vv[2:d+1,1:obs_gap:end]

# Build PSI
function InvBurgRK4_1step(x)
   lx = length(x)
   function F(x)
       𝑥 = [conj(reverse(x, dims = 1));0; x]
       conv(𝑥,𝑥)[2*lx+2:3*lx+1]
   end
   k1 = F(x)
   k2 = F(x .+ h*k1/2)
   k3 = F(x .+ h*k2/2)
   k4 = F(x .+ h*k3)
   A = @. x + h/6*(k1 + 2k2 + 2k3 + k4)
end

function Inertialman_part(x)
   lx = length(x)
   𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

   L = complex(zeros(lx^2))
   for j = 1:lx
      for k = 1:lx
         L[ (j-1)*lx+k] = 𝑥(j+lx)*𝑥(j+lx-k)
      end
   end
   L
end

Psi(x) = [x; InvBurgRK4_1step(x); Inertialman_part(x)]

# Get Wiener filter
#@time h_wf = get_wf(V_obs,Psi, M_out = M_out,PI = true)
signal = V_obs
rl = true
PI = false
rtol = 1e-6

sig = signal[:,2:end]
d, steps = size(sig)
nu = size(Psi(zeros(d,1)),1)

pred = complex(zeros(nu, steps))
for n = 1:steps
     pred[:,n] = Psi(signal[:,n])
end
pred

# h_wf = rl ? real(vector_wiener_filter_fft(pred, sig, M_out,
#                     PI = PI, rtol = rtol)) :
#                  vector_wiener_filter_fft(pred, sig, M_out,
#                     PI = PI, rtol = rtol)

par = 55
Nex = 2^10
win = "Par"

d, stepsy = size(sig)
nu, stepsx = size(pred)

stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])

Nexh = Int(floor(Nex/2))

L = par
lags = 0:L;

# Smoothed viewing window
if win == "Bar"
     lam = 1 .- (0:L)/L
elseif win == "Tuk"
     lam = .5*(1 .+ cos.(pi/L*(0:L)))
elseif win == "Par"
     LL = Int(floor(L/2))
     lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
     lam2 = 2*(1 .- (LL+1:L)/L).^3
     lam = [lam1; lam2]
else
     lam = ones(L+1)
end

R_pred_smoothed = complex(zeros(nu,nu,length(lags)))
for i = 1 : nu
     for j = 1 : nu
         R_pred_smoothed[i,j,:] = lam .* my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
     end
end

RR = R_pred_smoothed[:,:,1]
ishermitian(RR)

sqrt(RR)



























# Save Wienerfilter
dat = Dict("dat_h_wf" => h_wf)
Data = merge(paramaters,dat)
# save("../data/KSE_Data/KSE_sol_wienerfilter.jld",Data)
save("Data\\KSE_sol_wienerfilter.jld",Data)
println("Wiener filter saved")

# # Load old Wiener Filter
# LD = load("Data\\KSE_sol_wienerfilter.jld")
# h_wf = LD["dat_h_wf"]
# println("Wiener filter load








n_gap = 1


## Spatial grid and initial conditions:
x = P*(0:N-1)/N
u = g.(x)
v = fft(u)

## Precompute various ETDRK4 scalar quantities:
q = 2π/P*(0:N-1) # Wavenumbers
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
F = plan_fft(v_pad)
iF = plan_ifft(v_pad)

function NonLin(v)
    v_pad = [v; zeros(N)]
    nv = F*(real(iF*v_pad)).^2
    nv[1:N]
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
        u = real.(ifft(v))
        uu[:,ni] = u
        vv[:,ni] = v
        tt[ni] = t
    end
end
# Energy spectrum
EE = log.(abs.(vv).^2)


size(vv)

start = n_disc+1
uu[:,start:end], vv[:,start:end], tt[end-start+1]
