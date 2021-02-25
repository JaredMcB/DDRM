using JLD
using DSP # For conv function in Psi
using Dates


mrf = include("../Tools/WFMR_fast.jl")
at = include("../Tools/AnalysisToolbox.jl")


gen = "lin1e3"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")

sol_file = server ? "../../../data/KSE_Data/ks_sol_$gen.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"
println("Sol load location: " * sol_file)

@time vv = load(sol_file,"dat_vv")


d = 5
h = 0.1
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
M_out = 20
n = 3; p = 1500; par = 1500
ty = "bin"
xspec_est = "old"
nfft = 0
rl = true
Preds = false
PI = false
rtol = 1e-6
info = false
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
    "rtol" => rtol,
    "tm" => tm
)

Len = 500

@time h_wf = mrf.get_wf(signal[:,1:Len], Psi; M_out,verb = true)


wf_file = "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_wf_$gen-Len$Len.jld"
h_wf_slow = load(wf_file,"dat_h_wf")

maximum(abs.(h_wf - h_wf_slow))

sig = @view signal[:,2:end]
pred = mrf.get_pred(signal, Psi)

d, stepsy = size(sig)
nu, stepsx = size(pred)

win = "par"
stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])
nfft = nfft == 0 ? nextfastfft(steps) : nfft
nffth = nfft ÷ 2
L = min(par,steps-1)

R_pred_smoothed = mrf.matrix_autocov_seq(pred; L, steps, nu, win)

P = R_pred_smoothed
N_ckms = 10^5

using SparseArrays
using LinearAlgebra

d = size(P)[1];
m = size(P)[3] - 1

NN = reverse((@view P[:,:,2:m+1]),dims = 3)
Re = Rr = p0 = @view P[:,:,1]

F = sparse([[spzeros(d,d*(m-1)); sparse(I,d*(m-1),d*(m-1))] spzeros(d*m,d)])
h = sparse([spzeros(d,d*(m-1)) sparse(I,d,d)])

K = complex(zeros(d*m,d))
for i = 0 : m-1
    K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
end
FL = K
i = 0
errK = errR = 1
rtol = 1e-6
@time hL = h*FL
@time FL = F*FL

@time Rr_pinv = pinv(Rr,rtol = rtol)
@time Rr_inv  = inv(Rr)
@time Re_pinv = pinv(Re, rtol = rtol)
rtol = 1e-32
@time FL_RrhLt = FL * Rr_inv

@time FL_RrhLt_o = FL/Rr

maximum(abs.(FL_RrhLt - FL_RrhLt_o))

maximum(abs.(Rr_pinv - Rr_inv))

@time S = svd(Rr)
sort(S.S)

# Stopping criteria stuff
i += 1

@time FL_RrhLt = FL * Rr_pinv * hL'
hL_RrhLt = hL * Rr_pinv * hL'
errK = norm(FL_RrhLt)
errR = norm(hL_RrhLt)

FL -= K * Re_pinv * hL
K  -= FL_RrhLt
Rr -= hL' * Re_pinv * hL
Re -= hL_RrhLt
    end

maximum(abs.(FL_RrhLt - FL_RrhLt_o))
