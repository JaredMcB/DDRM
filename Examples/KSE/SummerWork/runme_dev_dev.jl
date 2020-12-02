
using JLD
using DSP: conv, nextfastfft # For conv function in Psi
using Dates

include("Model_KSE.jl")
include("Model_Reduction_Dev.jl")

gen = "_Lin"
T = 10^5 # Length (in seconds) of time of run
T_disc = Int(T/2) # Length (in seconds) of time discarded
P = 21.55  # Period
N = 96  # Number of fourier modes used
h = 1e-3 # Timestep
g = x -> randn()*cos(2π*x/P)*(randn() + sin.(2π*x/P))
q = 2π/P*(0:N-1)
obs_gap = 100
n = 3
p = 1500
par = 1500
d = 5 # No. of lowest modes taken in reduced model
M_out = 2000
short = false
loadsol = true
loadwf = false

# Find out if I'm on the server
server = pwd() == "/u5/jaredm/Server_scripts" ? true : false
println("on server = $server")
# set save destinations
sol_file = server ? "../data/KSE_Data/KSE_sol$gen.jld" :
                   "Data\\KSE_sol$gen.jld"
println("Sol save location: "*sol_file)
wf_file = server ? "../data/KSE_Data/KSE_wf$gen-Mo$M_out.jld" :
                   "Data\\KSE_wf$gen-Mo$M_out.jld"
println("WF save location: "*wf_file)
rmrun_file = server ? "../data/KSE_Data/KSE_rmrun$gen-Mo$M_out.jld" :
                   "Data\\KSE_rmrun$gen-Mo$M_out.jld"
println("redmodrun save location: "*rmrun_file)

paramaters = Dict(
   "gen" => gen,
   "T" => T,
   "T_disc" => T_disc,
   "P" => P,
   "N" => N,
   "h" => h,
   "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
   "q" => q,
   "d" => d,
   "tm" => now(),
   "M_out" => M_out,
   "n" => n,
   "p" => p,
   "par" => par,
   "obs_gap" => obs_gap,
   "d" => d, # No. of lowest modes taken in reduced model
   "M_out" => M_out,
   "short" => short,
   "loadsol" => loadsol,
   "loadwf" => loadwf)


println("the Parameters ===================")
for x in keys(paramaters)
   println(x*" : ", paramaters[x])
end
println("==================================")

if loadsol
   # Load Old Data
   @time Data = load(sol_file)
   print("Data Loaded")
   uu = Data["dat_uu"]
   vv = Data["dat_vv"]
   tt = Data["dat_tt"]
else
   uu, vv, tt =  my_KSE_solver(T,
            T_disc  = T_disc,
            P = P,
            N = N,
            h = h,
            g = g,
            n_gap = obs_gap)

   dat = Dict(
      "dat_uu" => uu,
      "dat_vv" => vv,
      "dat_tt" => tt)
   Data = merge(paramaters,dat)
   save(sol_file,Data)
   println("data saved")
end

## Get Reduced model #########################################################

# collect observations
V_obs = vv[2:d+1,1:end]
nfft = nextfastfft(size(V_obs,2))
M_out = min(M_out, nfft)

# Build PSI
function InvBurgRK4_1step(x)
  lx = length(x)
  function F(x)
      𝑥 = [conj(reverse(x, dims = 1));0; x]
      -im/2*(2π/P*(1:lx)/N) .* conv(𝑥,𝑥)[2*lx+2:3*lx+1]
  end

  Δt = h*obs_gap

  k1 = F(x)
  k2 = F(x .+ Δt*k1/2)
  k3 = F(x .+ Δt*k2/2)
  k4 = F(x .+ Δt*k3)
  A =  @. x + Δt/6*(k1 + 2k2 + 2k3 + k4)
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

function Inertialman_part_short(x)
   lx = length(x)
   𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

   L = complex(zeros(binomial(lx+1,2)))
   i = 1
   for j = 1:lx
      for k = j:lx # k should normaly go from 1 to lx but i changed it to go from j to lx.
         L[i] = 𝑥(j+lx)*𝑥(j+lx-k)
         i += 1
      end
   end
   L
end

Psi(x) = short ? [x; InvBurgRK4_1step(x); Inertialman_part_short(x)] :
                 [x; InvBurgRK4_1step(x); Inertialman_part(x)]


# Get Wiener filter
if loadwf
   # Load old Wiener Filter
   LD = load(wf_file)
   short == LD["short"] || throw(DimensionMismatch("wf and psi don't agree"))
   h_wf = LD["dat_h_wf"]
   println("Wiener filter load")
else
   print("Get_wf computation time: ")
   @time h_wf = get_wf(V_obs,Psi,
      M_out = M_out,
      n = n, p = p,
      par = par,
      rl = false,
      PI = false)

   # Save Wienerfilter
   dat = Dict("dat_h_wf" => h_wf)
   Data = merge(paramaters,dat)
   save(wf_file,Data)
   println("Wiener filter saved")
end

## Run reduced model with no noise
steps = size(V_obs,2)
V_rm = [V_obs[:,1:M_out] complex(zeros(d,steps-M_out))]
nu = size(Psi(V_obs[:,1]),1)

# load presamples
PSI_past = complex(zeros(nu,steps))
for i=1:M_out
    PSI_past[:,i] = Psi(V_obs[:,i])
end

# Move forward without original data
print("Reduced Model Run Time: ")
@time for i = M_out+1:steps
    V_rm[:,i] = sum(h_wf[:,:,k]*PSI_past[:,i-k] for k = 1:M_out)
    isnan(V_rm[1,i]) && break
    PSI_past[:,i] = Psi(V_rm[:,i])
end

# Process reduced model run
steps_rm = size(V_rm,2)
V_rm_end = conj(reverse(V_rm,dims = 1))
VV_rm = [zeros(1,steps_rm); V_rm; zeros(N-2d-1,steps_rm);V_rm_end]
UU_rm = real(ifft(VV_rm,1))
tt_rm = tt[1:end]

dat = Dict(
   "dat_UU_rm" => UU_rm,
   "dat_VV_rm" => VV_rm,
   "dat_tt_rm" => tt_rm)
Data = merge(paramaters,dat)
save(rmrun_file,Data)
println("Reduced Model Run saved")
