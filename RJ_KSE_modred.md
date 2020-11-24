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

11:38 PM - After working on getting the experiment to go through, off and on all evening, I located the problem. I realized that I had done something sort of silly. This was the code in the function `spectfact_matrix_CKMS` in "Tools/Model_Reduction_Dev.jl"

```julia
F = sparse([[zeros(d,d*(m-1)); I] zeros(d*m,d)])
h = sparse([zeros(d,d*(m-1)) I])
```

So that before I gained any benefit from the sparse functionality it needed to store a traditional matrix with `dm × dm = (25)(5001) × (25)(5001) = 15,631,250,625`
entries, each entry being `float64` (over 120 GB). So, I change those lines to:

```julia
F = sparse([[spzeros(d,d*(m-1)); sparse(I,d*(m-1),d*(m-1))] spzeros(d*m,d)])
h = sparse([spzeros(d,d*(m-1)) sparse(I,d,d)])
```

And this fixed it for me, I can get the filter without the memory error. And here is the wiener filter I get: (it is saved in "Examples/KSE/Data/KSE_wf11_23_20_1-Mo20.jld")

```
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  2920.78+86.0131im  -217.362-6.38956im  -173.838+47.5813im  -111.006+2.8904im    372.945-32.5095im  -2924.89-56.8678im  …    -0.0015088+0.000684852im   0.000923013+0.000153075im   0.00142293+0.00194967im
  -356.12+557.155im   2902.06+122.798im  -391.346+61.9222im   784.721+112.349im  -150.986+18.2402im   751.572-541.415im     -0.000434052-0.00210347im   -0.000490083+0.00367567im    0.00426563-0.00315144im
 -341.911-74.2766im  -94.8209+71.2338im   3307.71+179.298im  -754.658+83.3726im   1316.66+11.0441im   1222.65+202.034im      -0.00377062-0.000221395im    0.00159182+0.00189452im     0.0130932-0.00136477im
   209.17-143.634im   197.418+37.1193im  -176.397+11.8324im   1196.92+116.858im   136.348-94.5001im    41.066+147.243im      0.000525555-0.00125459im    0.000378801+0.000326265im  -0.00071374-9.65751e-5im
 -8.61591-33.4598im  -19.7647+34.7131im   11.5672-16.5271im  -33.1281+40.0582im  -38.4866+22.9555im   3.17366+51.8029im      -8.57798e-5-0.000443054im   -0.00012088-0.000119477im   0.00121908-0.00104132im

[:, :, 2] =
 -7759.88-236.393im   627.664+7.43828im   454.368-130.546im    268.91+3.53185im  -973.635+75.8821im   7799.94+161.334im  …   0.00523096-0.00223101im  -0.00296975-0.000640396im  -0.00491493-0.00684972im
   1193.5-1602.2im   -7578.86-332.4im     1077.58-238.184im  -2185.59-364.206im   479.451-7.19992im  -2241.97+1586.67im      0.00157828+0.00724294im   0.00164119-0.0118443im     -0.0152026+0.0106635im
  1084.48+202.078im   346.952-305.459im  -8937.85-558.427im    2044.7-217.277im  -3458.33+7.27727im  -3330.49-556.14im        0.0134938+0.00119259im  -0.00522827-0.00619827im    -0.0463804+0.00392117im
 -510.371+388.588im  -582.244-117.914im   523.068-45.546im   -2904.42-345.602im  -350.641+260.528im  -58.6537-367.444im     -0.00183893+0.00419061im  -0.00124291-0.00103709im    0.00255997+0.000264235im
   37.392+93.4376im   53.0811-94.2821im  -30.6937+47.3465im   88.5567-99.8301im   147.114-53.0441im  -19.9734-142.918im     0.000264682+0.00159728im  0.000391158+0.000383516im  -0.00433801+0.00347736im

[:, :, 3] =
  5350.89+170.584im  -483.071+6.41652im  -283.844+95.8814im  -138.705-16.1145im   624.176-47.7117im  -5416.91-123.801im  …   -0.00561517+0.00224266im    0.00277731+0.00082384im   0.00524701+0.0075443im
 -1106.74+1264.89im   5101.91+231.883im  -773.819+282.258im   1566.78+335.031im  -421.983-63.836im    1812.06-1292.95im      -0.00186083-0.00773265im   -0.00167567+0.0113642im     0.0170725-0.011323im
 -875.357-153.443im  -338.745+368.925im   6403.77+477.529im  -1391.52+146.515im   2215.17-56.0657im   2279.09+410.525im       -0.0150817-0.00187566im    0.00506962+0.00594193im    0.0515195-0.0034035im
  315.193-256.752im   445.416+103.85im   -413.394+56.8152im   1756.98+277.483im   228.016-189.623im  -30.4105+209.335im       0.00197981-0.00429985im    0.00120888+0.000959219im  -0.0028879-0.000243004im
 -39.3075-69.9776im  -35.7928+66.5243im   22.3769-36.5236im  -59.7104+59.4334im  -131.261+25.5174im   26.6048+105.603im     -0.000239227-0.00181095im  -0.000382883-0.00036542im    0.0048411-0.00361176im

...

[:, :, 18] =
 -82.0288-3.92259im    9.93911+0.2223im      9.29949+0.449136im    7.4362+0.911555im  -14.8361-2.27226im    80.9818+3.67635im    …   -8.04113e-5+8.10665e-6im    4.07315e-5+2.91395e-5im    6.26269e-5+0.000106657im      
  21.2101-10.8716im   -73.5278+2.0558im      11.4004+5.2175im    -43.6966-1.11043im    9.62675-6.81308im   -34.9782+8.75682im         1.66236e-5-0.000104063im  -2.63923e-5+0.000177502im   9.21832e-5-8.00233e-5im       
  37.5627+9.76838im    1.85439+0.753894im   -84.9788-9.51861im     35.078+1.48597im   -60.5251-4.55755im   -66.5973-13.9621im       -0.000174641-1.87823e-5im    5.38448e-5+0.000106241im  0.000463127-1.41534e-5im       
  1.80034+11.9169im   -10.3489-4.41569im     11.0303+4.11451im   -10.5075-1.59132im   -2.22282+4.35251im   -2.05436-11.5135im         2.50605e-5-5.79293e-5im    2.18811e-5+1.12321e-5im     1.4162e-6+9.45823e-7im       
  1.78383-0.505166im  -1.91347-1.27441im   0.0430545-0.863538im  0.300517-1.42996im    1.82886-0.963437im  -1.62982+0.0791231im       1.61193e-6-1.44552e-5im   -8.79652e-6-5.34033e-6im    4.05752e-5-3.88802e-5im       

[:, :, 19] =
 -73.9919-2.65793im   5.40608+0.562001im    6.99861-1.00879im     5.47576-0.630584im   …  -1.27101e-6+1.13302e-5im   6.77226e-6-6.62535e-6im  5.09724e-5+1.62828e-5im     1.98787e-5-3.05382e-5im
 0.169747-2.68968im  -70.0764-1.45003im     8.13096+7.10543im    -33.0754-0.582066im       4.53443e-6-7.18842e-6im   2.04146e-5-1.60579e-5im  -1.3893e-5+0.000166537im  -0.000102771+0.000108316im
   22.985-5.60736im   3.55963+13.4174im    -63.6217+1.01673im      24.055+1.36557im       -7.26848e-6+1.40743e-6im  -2.05499e-5+4.44906e-6im  8.32945e-5+0.000139682im  -0.000243081+4.71682e-5im
 0.311492+6.68581im  -11.4999+0.257542im     7.3424+0.525371im   -17.2546+0.0196872im      5.07448e-7-1.27567e-5im   7.71129e-8-6.11894e-6im  7.51987e-6+2.02777e-5im     2.37272e-5+9.74479e-6im
 0.126487+0.70794im   1.81071-1.07538im   0.0282389+0.0543418im   1.56894-1.7382im         5.87888e-6+1.61594e-6im  -1.06219e-5+8.49857e-6im  1.47784e-6-2.42374e-6im    -1.99084e-5+1.41607e-5im

[:, :, 20] =
 -24.7332+0.10752im    0.946861+0.231071im  0.00273346-3.14355im    2.44789-0.72574im  -3.79639+1.26521im   …  -4.36156e-6+9.91228e-6im   5.60006e-5-5.98667e-6im   4.84142e-5-9.3454e-6im    -1.98256e-5-0.000146238im   
 -26.1412-0.411731im   -34.6857-8.15459im       1.6021+2.8524im    -7.95984-2.29227im  -6.14417+1.43367im       8.75043e-6-9.37365e-6im   3.23351e-5+9.90862e-5im   2.18319e-5+7.04962e-5im  -0.000306368+0.000113951im   
 -50.4341+0.0650817im   -1.2381+7.98817im     -13.0444+0.336579im  0.366277-1.69105im  -14.1211-1.99876im       9.19965e-6-6.95808e-6im   0.00017311+9.11037e-5im   4.01433e-5+9.42063e-5im  -0.000674968-0.000119863im
  3.85341+0.84301im    -2.40369-0.725644im    -2.93544+2.57677im   -13.7484-2.01821im  -4.01097+3.17342im      -1.36109e-6-1.39479e-5im  -1.47302e-5+3.15551e-5im  -1.68483e-6+1.4567e-5im    -2.00679e-5+4.08112e-5im    
 -1.05924-1.41698im    0.537742-0.882067im   -0.225708-0.503753im  0.568337-1.09505im  0.864916-0.911282im      5.36164e-6+1.65332e-6im  -3.32195e-6+2.82276e-5im   1.25603e-6-2.07771e-6im   -7.82065e-5+3.80681e-5im
 ```
