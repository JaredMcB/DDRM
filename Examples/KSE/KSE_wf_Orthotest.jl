
using JLD
using DSP # For conv function in Psi
using Dates

mr = include("../../Tools/WFMR.jl")
at = include("../../Tools/AnalysisToolbox.jl")

# Load Old Data

gen = "lin1e3"     # this is just a reference designation it shows up in the
                   # output file. I think of generatrion.

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")

sol_file = server ? "../../../data/KSE_Data/ks_sol_$gen.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"
println("Sol load location: " * sol_file)

@time vv = load(sol_file,"dat_vv")

## Get Reduced model #########################################################
# Model reductrion parameters

d = 5
h = 0.1
# collect observations
obs_gap = 1
V_obs = vv[2:d+1,1:obs_gap:end]
vv = []

# Build PSI
function InvBurgRK4_1step(x)
   lx = length(x)
   function F(x)
       𝑥 = [conj(@view x[lx:-1:1]) ;0; x]
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

signal = V_obs
V_obs = []
M_out = 20
n = 3; p = 1500; par = 1500
ty = "bin"
xspec_est = "old"
nfft = 0
rl = true
Preds = false
N_ckms = 3000
PI = false
rtol = 1e-6
Verb = false
tm = now()

paramaters = Dict(
    "M_out" => M_out,
    "n" => n,
    "p" => p,
    "par" => par,
    "ty" => ty,
    "xspec_est" => xspec_est,
    "nfft" => nfft,
    "rl" => rl,
    "Preds" => Preds,
    "N_ckms" => N_ckms,
    "rtol" => rtol,
    "tm" => tm
)

d, steps = size(signal)
nu = size(Psi(zeros(d,1)),1)

sig = @view signal[:,2:steps]   # sig is now one a head of signal
steps -= 1                      # this makes steps the length of sig
pred = mr.get_pred(signal,Psi)

wf_file = server ? "../../../data/KSE_Data/ks_wf_$gen-Len2000.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_wf_$gen-Len2000.jld"
println("WF load location: " * wf_file)

@time h_wf = load(wf_file,"dat_h_wf")

# Pred is ment to estimate sig
sig_hat = zeros(ComplexF64, d, steps)
sig_hat[:,1:M_out-1] = @view sig[:,1:M_out-1]
for i = M_out : 1000
    @views sig_hat[:,i] = sum(h_wf[:,:,k]*pred[:,i-k+1] for k = 1: M_out)
end

Err = sig - sig_hat

C = zeros()





using PyPlot
plot(rand(10))
for dd = 1:1
    plot(real(sig_hat[dd,1:1000]),label = "real hat $dd")
    plot(real(sig[dd,1:1000]),label = "real $dd")
    legend()
end
sig[:,1:100]
