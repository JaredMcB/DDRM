include("KSE_modredrun.jl")

## Parameters for KSE model
gen = 1
T = 400 # Length (in seconds) of time of run
T_disc = 100 # Length (in seconds) of time discarded
P = 32π  # Period
N = 128  # Number of fourier modes used
h = 1e-2 # Timestep
g = x -> cos(π*x/16)*(1 + sin.(π*x/16))
q = 2π/P*(0:N-1)
obs_gap = 10
d = 5 # No. of lowest modes taken in reduced model
M_out = 10 # No. of coeficinets in Wiener filter output

KSE_modredrun(
   gen = gen,
   T = T, # Length (in seconds) of time of run
   T_disc = T_disc, # Length (in seconds) of time discarded
   P = P,  # Period
   N = N,  # Number of fourier modes used
   h = h, # Timestep
   g = g,
   q = q,
   obs_gap = obs_gap,
   d = d, # No. of lowest modes taken in reduced model
   M_out = M_out,
   short = false,
   loadsol = true,
   loadwf = false)
