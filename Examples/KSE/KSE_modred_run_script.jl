
using JLD
using DSP # For conv function in Psi
using Dates

mr = include("../../Tools/WFMR.jl")
at = include("../../Tools/AnalysisToolbox.jl")

# Load Old Data

gen = "lin1e5"     # this is just a reference designation it shows up in the
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

Len = 5000

h_wf = @time mr.get_wf(signal[:,1:Len], Psi; M_out, verb = true)

# Save Wienerfilter
dat = Dict("dat_h_wf" => h_wf)
Data = merge(paramaters,dat)
# save("../data/KSE_Data/KSE_sol_wienerfilter.jld",Data)

wf_file = server ? "../../../data/KSE_Data/ks_wf_$gen-Len$Len.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_wf_$gen-Len$Len.jld"
save(wf_file,Data)
println("Wiener filter saved")

# # Load old Wiener Filter
# LD = load("Data\\KSE_sol_wienerfilter.jld")
# h_wf = LD["dat_h_wf"]
# println("Wiener filter load
