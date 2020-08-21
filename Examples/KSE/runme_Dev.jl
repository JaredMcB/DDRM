include("KSE_modredrun_Dev.jl")

## Parameters for KSE model
gen = "_Lin"
T = 10^5 # Length (in seconds) of time of run
T_disc = Int(T/2) # Length (in seconds) of time discarded
P = 21.55  # Period
N = 96  # Number of fourier modes used
h = 1e-3 # Timestep
g = x -> randn()*cos(2π*x/P)*(randn() + sin.(2π*x/P))
obs_gap = 100
n = 3
p = 1500
par = 1500
d = 5 # No. of lowest modes taken in reduced model
M_out = 10000 # No. of coeficinets in Wiener filter output

KSE_modredrun(
   gen = gen,
   T = T, # Length (in seconds) of time of run
   T_disc = T_disc, # Length (in seconds) of time discarded
   P = P,  # Period
   N = N,  # Number of fourier modes used
   h = h, # Timestep
   g = g,
   obs_gap = obs_gap,
   d = d, # No. of lowest modes taken in reduced model
   M_out = M_out,
   n = n, p = p,
   par = par,
   short = false,
   loadsol = false,
   loadwf = false)
