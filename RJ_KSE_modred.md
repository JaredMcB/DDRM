# Monday, November 23, 2020

11:31 AM - Here I start again, (though much wiser now) my study of model reduction by the Weiner projection of the KSE model. So, to get started I will need some data. I will us the data I produced last summer from parameters used by Lu, Lin and Chorin in

Lu, Fei, Kevin K. Lin, and Alexandre J. Chorin. "Data-based stochastic model reduction for the Kuramoto–Sivashinsky equation." *Physica D: Nonlinear Phenomena* 340 (2017): 46-57.

These values can be found in Sec. 5.1
- L = 2π/sqrt(0.085) ≈ 21.55 (this is the period)
- N = 32*floor(Int,L/(2\pi)) = 96 (number of modes used in the full model run)
- dt = 0.001 (time step of full model run)
- v0(x) = (1 + sin(x))*cos(x) (initial value)
- gap = 100
- In the paper the model was run for 6×10^4 time units and the first 10^4 were discarded. Here we will run it for 10^5 time units (10^8 steps) and discard the first half.   

The data was run and saved in the summer. I just today loaded the old data and plot it and it looks good. The data is in "Examples/KSE/Data/KSE_sol_lin.jld" on my laptop. Since, `.jld` files are in `.gitignore` it would be nice to have a set of this data on the other machines I work with. The data was generated using `my_KSE_solver` from `Model_KSE.jl` in "Examples/KSE". Particularly with the following call:
```julia
T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2, # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085),  # Period
N        = 96,  # Number of fourier modes used
h        = 1e-3, # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16)),
obs_gap  = 100

uu, vv, tt =  my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)
```

I suppose setting the seed and verifying that the data is the same may be good too. So,
we will want
```julia
using Random
seed = 2020
Random.seed!(seed)
```
running the last command immediately before the above call.

This is how I obtained the data that I will be using for the initial exploration of the Wiener filtering algorithm on KSE data.

12:48 PM - Now to get the Wiener filter for this data. I will proceed as I did in the summer and us the 5 lowest (nonconstant) modes. That is
```julia
X = vv[2:6,:]
```
This give a 5 by 500,001 data timeseries.

#### Experiment Nov 23, 2020 1 (Wiener filter for Lin data)
The script for this experiment was written fresh for the occasion, but used basically the top half of the function 'KSE_mdredrun' in "KSE/KSE_mdredrun.jl".  It is saved as "KSE/KSE_wf.jl".

The data is as described above.

```julia
Exp = 11_23_20_1

## Parameters
# run parameters
T        = 10^5              # Length (in seconds) of time of run
T_disc   = T÷2               # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)    # Period
N        = 96                # Number of fourier modes used
h        = 1e-3              # Timestep
g        = x -> cos(x)*(1 + sin.(x))

q        = 2π/P*(0:N-1)      # Wve numbers (derivative)

# Observation parameters
obs_gap  = 100
d        = 5 # No. of lowest modes taken in reduced model

# Wiener filtering parameters
M_out    = 20
nfft     = 2^14
par      = 5000
xspec_est= "old" # Default
short    = true
loadsol  = true
```

I ran this experiment on my laptop and got a out-of-memory error. I will now conduction a memory requirement analysis.

2:49 PM - I created the file "KSE/KSE_data_gen.jl." It just generates and saves KSE runs. I ran it just now according to the above parameters. Displayed below for conveniences.
```julia
T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2 # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100
seed     = 2020
```
The run was saved in "Examples/KSE/Data/KSE_sol_lin.jld" on my desktop. 
