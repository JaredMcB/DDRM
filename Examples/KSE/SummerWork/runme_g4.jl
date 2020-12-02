include("KSE_modredrun.jl")

## Parameters for KSE model
gen = 4
T = 400 # Length (in seconds) of time of run
T_disc = 100 # Length (in seconds) of time discarded
P = 2π  # Period
N = 16  # Number of fourier modes used
h = 1e-2 # Timestep
g = x -> cos(x)*(1 + sin.(x))
q = 2π/P*(0:N-1)
obs_gap = 1
d = 5 # No. of lowest modes taken in reduced model
M_out = 10 # No. of coeficinets in Wiener filter output

@time KSE_modredrun(
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
   short = true,
   loadsol = false,
   loadwf = false)
