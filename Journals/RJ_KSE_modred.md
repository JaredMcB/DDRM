Nonlinear# Monday, November 23, 2020

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


# Tuesday, November 24, 2020

9:55 AM - Yesterday, I got a result for the Wiener filter using KSE data and short psi predictors. The Wiener filter was suspect, though, as it had coefficients in the thousands. Today, I want to vary the parameters a bit and see if I can get a more reasonable output for the Wiener filter. I will run these now on thelio because it takes my computer a very long time.

So, first I will rerun *Experiment Nov 23, 2020 1* on thelio.

10:31 AM - Ran "Examples/KSE/KSE_data_gen.jl" (job 157) to generate KSE data and save it as "data/KSE_Data/KSE_sol_lin.jld". Here, `gen = "lin"`. So, now we have a copy of the data of thelio, with the standard `gen = "lin"` parameters:

```julia
gen = "lin"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2 # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100
seed     = 2020
```


#### Experiment Nov 24, 2020 1 (job 158) (WF from Lin data on thelio)
Here I am just verifying the Wiener filter I got on my desktop, on thelio. To do this I ran the script "Examples/KSE/KSE_wf.jl". the only modification of this script from yesterday is the experiment ID `Exp = "11_24_20_1"`. As a reminded the WF parameters are:

```julia
# Wiener filtering parameters
M_out    = 20
nfft     = 2^14
par      = 5000
xspec_est= "old" # Default
```
 Started running at: `12:00:00`
 received email at: `1:08 PM`

```
 julia> h_wf = load("KSE_wf11_24_20_1-Mo20.jld","dat_h_wf")
 5×25×20 Array{Complex{Float64},3}:
 [:, :, 1] =
    3109.9+84.7565im  -39.7975-77.5235im   7.11213+49.9121im  -118.967-49.9446im  -501.476+0.966074im  …  0.000549076+3.44095e-5im   -6.97895e-5+0.000344297im  0.000566329-0.000193866im
  -1148.84-30.6693im   3232.96+128.743im   -309.17-200.155im  -345.561+112.045im  -558.048-22.5339im       -0.0013568-0.00172819im   0.000896153-0.00181537im   -0.00935467+0.00513726im
  -252.427-124.57im   -566.325+44.6072im   3168.12-31.7975im  -795.521-115.311im   1148.17+202.584im       0.00283151-0.000930736im   0.00167334+0.00157873im   -0.00511643-0.00830197im
  -313.978+142.862im  -198.487+134.461im   -224.11-155.008im   1126.97+84.6742im  -542.366-21.3349im      0.000363192+0.000308685im  0.000418992-0.000837426im  -0.00275977+0.00427073im
  -63.7546+72.5455im   18.5082+2.35742im  -24.4098-10.0775im   4.88718+33.977im    -7.9182-33.0331im       -9.8937e-5+8.42246e-5im    0.00033948-0.000124782im   5.69483e-5-0.000957531im

 [:, :, 2] =
  -8241.82-242.109im   135.053+207.524im  -12.9576-119.311im   301.812+121.994im    1387.8-19.1814im  …  -0.00195107-0.000174871im  0.000243494-0.00115168im    -0.00212351+0.000514835im
   3176.33+57.9795im  -8546.96-375.495im   907.993+582.686im   863.419-305.348im   1514.91+63.7025im      0.00487494+0.00629004im   -0.00286146+0.00560815im      0.0335742-0.0167639im
   712.358+422.22im    1725.28-176.239im  -8548.54+113.918im    2320.7+306.155im  -3242.37-613.841im      -0.0100169+0.00305788im    -0.0055161-0.00517093im      0.0180225+0.0296957im
   755.089-408.552im   502.806-364.732im   586.962+433.273im  -2765.08-251.093im   1363.19+61.1068im     -0.00129827-0.00115494im   -0.00127364+0.00284297im     0.00966918-0.0146489im
   162.213-193.49im   -51.2965-12.1497im   69.0221+31.1499im  -19.8986-87.379im    65.3114+87.3529im     0.000358549-0.000279009im  -0.00112634+0.000366455im  -0.000158865+0.00337956im

 [:, :, 3] =
    5639.9+176.303im  -126.541-148.986im   2.71145+65.3569im  -186.128-76.8613im  -999.969+29.4081im  …    0.00217686+0.000274976im  -0.000271429+0.00116757im    0.002519-0.000310408im
  -2297.08-19.0417im   5867.16+288.83im   -717.542-468.913im   -506.85+213.04im   -1064.88-45.3805im      -0.00551137-0.00724598im     0.00261658-0.00502219im  -0.0379232+0.0169118im
  -533.702-398.15im   -1435.49+230.13im    6062.35-123.22im   -1807.65-207.598im   2385.44+508.092im        0.0111752-0.00312159im     0.00535357+0.00497043im  -0.0201502-0.0335603im
  -427.501+304.684im  -299.637+253.99im   -388.203-317.879im   1692.25+194.269im   -839.64-40.0804im       0.00148482+0.00138537im     0.00107456-0.00297393im  -0.0105005+0.0157904im
  -101.092+132.772im   36.2929+17.125im   -52.5339-26.3798im    19.838+56.4945im   -71.346-59.6032im     -0.000410404+0.000282763im    0.00112048-0.00031188im  8.36847e-5-0.00376227im

 ...

 [:, :, 18] =
  -104.542-5.45384im     2.62601-0.208246im  0.0518849-1.35528im     3.82022-0.154201im   21.2783-0.988257im  …   2.35915e-5+1.1028e-5im   -5.20795e-6+1.7476e-5im      2.85348e-5+2.08824e-5im
   44.6575-0.521279im   -82.9852-10.8645im     14.3028+1.90606im     13.2923-3.58912im      16.02+1.43101im      -0.00010039-5.72663e-5im   5.25571e-5-6.6691e-5im    -0.000546477+0.000186319im
   6.24042+1.34578im      36.388+5.20176im    -97.2304+0.123756im    40.3401+4.26719im   -62.1176-7.31098im       4.47861e-5-2.45649e-5im   8.41217e-5+0.000120517im   -0.00022193-0.000314449im
   5.19614-4.157im       11.1009-3.26051im      10.516+6.05245im    -22.2856-3.73709im    11.0738+1.72661im       5.55222e-7+4.29189e-6im   2.13319e-5-4.37111e-5im   -0.000192747+8.24675e-5im
   1.72253-2.77361im   -0.704536+0.47509im     1.09941+0.240471im  -0.625231-0.597646im   1.48926+1.03673im      -6.03275e-6+2.6377e-6im    1.86673e-5-1.30561e-6im     9.19437e-6-3.00916e-5im

 [:, :, 19] =
  -72.5065-4.56631im   0.595983+5.4208im     -0.616405-3.37867im     3.61536+0.0379607im   12.1755-1.0053im     …  -1.00548e-5-2.06446e-5im  1.94206e-6+1.29777e-5im   -1.77571e-5-6.3334e-5im
   24.9339+2.13021im   -87.4924-2.37502im      3.58082+1.82149im      15.659-4.02461im     16.2158+0.0504848im      4.36741e-5+5.04838e-5im    6.475e-5-8.13762e-5im   0.000381882-3.98933e-5im
  0.405442-6.28824im   -3.73396+0.0867596im   -56.3007-6.56403im      21.719+2.65926im     -29.514-3.78666im       -3.99717e-5-3.74547e-6im  9.26336e-5+0.000127987im  0.000247422+0.000278405im
   7.61729-3.68829im    6.81353-4.44026im       4.7274+3.24823im    -23.0566-3.96673im     12.1187+1.65507im       -6.31901e-6-4.25119e-6im  2.97968e-5-9.54946e-6im    8.40925e-5-6.26588e-5im
   2.52379-1.15261im  -0.649033-0.0867026im   0.345576-0.131901im  0.0063477-0.604838im   0.874488+1.28978im        3.37441e-6+2.17833e-6im  1.26601e-5-1.2747e-6im     1.18745e-5+3.30193e-5im

 [:, :, 20] =
  -30.9002-1.29305im  -0.798362-3.37331im   -1.73166-1.36627im    2.92803+0.11047im    1.81021+3.19735im   …  -3.27863e-5+3.14469e-6im  6.53143e-6+1.31502e-6im  -5.00288e-5+4.67293e-5im
   5.40012-2.39241im   -32.8446-3.03259im   -8.18909-1.96476im    11.4529-1.39188im    3.20283-2.07963im       1.95392e-5+0.0001478im   6.27334e-5-6.51242e-5im  0.000480937-2.58763e-5im
  -4.45359-7.63609im    -9.4026-2.08173im   -8.07437-4.1751im    -5.73963-0.884132im   2.88312+1.33092im      -0.00017352+6.12357e-6im  6.56691e-5+8.80677e-5im   0.00013205+0.000515736im
   7.15529-3.7861im     6.95083+1.78727im    2.44135-1.90392im   -19.4209-2.29633im    9.81436-1.56204im      -2.82034e-5-3.21568e-5im  3.09948e-5+2.25751e-5im  0.000141767-0.000220439im
   1.61407-1.19554im   0.273568+0.191752im  0.325172-0.333677im  0.189027-0.576582im  0.291242+0.930842im      6.82463e-6-2.39109e-6im  3.95178e-6-4.47809e-7im  -1.26377e-5+4.83518e-5im
 ```

 Here are more output results:

 ```
Get_wf computation time: 4137.379482 seconds (841.26 M allocations: 8.396 TiB, 7.69% gc time)

CKMS Computation time: 3898.175304 seconds (556.87 M allocations: 8.307 TiB, 8.03% gc time)
Number of CKMS iterations: 20236
errK errR : 9.974873558782727e-11 3.1617462324742043e-15
```


I want to rerun the above experiment, now with more principled choices of parameters.  So, it will be very simmilar to *Experiment Nov 24, 2020 1*. The difference being in the WF parameters, namely `par` and `nfft`. These will be determined based on characteristics of the data. Such as the time it takes the autocorrolation to decay to a 2 margins of error.

So, once the data set was completed on thelio, I created the notebook "KSE_data_analyzer" in "Examples/KSE". In it I load the data, get the `sig` and `pred` arrays and then look at there autocovariances and cross covariances.

Here is what I looked at specifically:
1. The autocovariance sequence for each of the 25 predictors. I just did:
   ```
lags = -4000:4000
let plt = plot
    for m = 1:25
        figure()
        A = at.my_crosscov(pred[m,:],pred[m,:],lags)
        plt(Δt*lags,A)
        legend()
        title("modes $m and $m")
    end
end
   ```
   I noticed that the shapes of the first 5 and the next 5 were very similar respectively (I mean 1 looked like 6, but scaled, 2 looked like 7, but scaled, and so on). The next 15 had extremely high

2. I began to think about the fact that these are not mean zero processes and how that effects model reduction in this way.

5:21 PM - Done for the day.



# Wednesday, November 25, 2020
2:47 PM - Today I have been investigating the question of implementing model reduction with noncentered data. Since this uses Wiener filtering theory and Weiner filters assumes the processes are mean zero, there is something to be done here. Is it so simple as centering `X` and `Psi(X)` taking the wiener filter and then adding means appropriately to the reduced model?

4;09 PM - It seems to be exactly that simple. The code as it is only contemplates estimated second order information from the time series, the spectral density of  `pred` and the cross spectral density of `sig` and `pred`. De-meaning in these functions is built-it. So, the Wiener filter is spatial shift invariant. The thing to do while running the reduced model is to add the mean of `sig` and subtract `h_wf` acting on the mean of `pred`.
this means running the reduce model should look like this:

```julia
d, steps = size(X)
nu = size(Psi(zeros(d,1)),1)

sig  = X[:,2:end]

pred = complex(zeros(nu, steps-1))
for n = 1:steps-1
    pred[:,n] = Psi(X[:,n])
end

PSI  = complex(zeros(nu,steps))
X_rm = complex(zeros(d,steps))
X_rm[:,1:M_out] = X[:,1:M_out]

pred_m = mean(pred,dims = 2)
sig_m  = mean(sig,dims = 2)

C = sig_m - sum(h_wf[:,:,k]*pred_m for k = 1:M_out,dims = 2)

for i = 1:M_out
    PSI[:,i] = Psi(X_rm[:,i])
end

for i = M_out+1:steps
   X_rm[:,i] = C + sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out,dims = 2) + Noise[:,i]
   PSI[:,i] = Psi(X_rm[:,i])
end
return X_rm
```

I verified this experimentally just to be sure. I used the following potentials with the following derivatives:
```julia
V_c_prime  = x -> -x.*(x.^2 .- 1)
V_s_prime  = x -> -x.^3 + 6x.^2 - 11x .+ 6
```
the first is centered (at zero) and the other has an axis of symmetry at `x=2`. I used this with the following parameters to generate two time series:
```julia
# Model run Parameters
sigma    = [.4]
sig_init = [1.5]
# Numerical estimate parameters
scheme   = "FE"
steps    = 10^7 # Nomber of time steps (not including those discarded)
h        = .1
discard  = steps # Number of time steps discarded
gap      = 1

V_c_prime  = x -> -x.*(x.^2 .- 1)
V_s_prime  = x -> -x.^3 + 6x.^2 - 11x .+ 6


# Get full model run
Random.seed!(2014)
X_c = @time gen.DataGen_DWOL(; sigma, V_prime = V_c_prime, sig_init, scheme, steps, h, discard, gap)

Random.seed!(2014)
X_s = @time gen.DataGen_DWOL(; sigma, V_prime = V_s_prime, sig_init, scheme, steps, h, discard, gap)
```

The series were identical, but for a shift of 2 in the x-direction. I then computed the Wiener filters using

```julia
Psi(x)  = [x; x.^2; x.^3]

par   = 5000
nfft  = 2^14
Preds = true
M_out = 20

h_wf_c, pred_c = mr.get_wf(X_c,Psi; M_out, par, nfft, Preds);
h_wf_s, pred_s = mr.get_wf(X_s,Psi; M_out, par, nfft, Preds);
```
and got the following results:
```
julia> h_wf_c
1×3×20 Array{Float64,3}:
[:, :, 1] =
 1.09548  -0.00017354  -0.0983164

[:, :, 2] =
 0.00221948  0.000275644  -0.000893722

[:, :, 3] =
 0.00053129  -0.000162947  -0.000390774

...


julia> h_wf_s
1×3×20 Array{Float64,3}:
[:, :, 1] =
 -0.0836275  0.589727  -0.098317

[:, :, 2] =
 -0.00960635  0.00563605  -0.00089298

[:, :, 3] =
 -0.0035093  0.00218446  -0.000391526

...
```
Observe that these both agree with the what is expected since the model for each (using forward Euler) is
```
X_c[:,n+1] = -h*X_c[:,n].^3 + (1 + h)*X_c[:,n] + Noise
X_s[:,n+1] = -h*X_s[:,n].^3 + 6h*X_s[:,n].^2 + (1 - 11h)*X_s[:,n] .+ 6h + Noise
```
this leads an analytic of
```
h_wf_c_ana[:,:,1] = [1.1   0    -0.1]
h_wf_c_ana[:,:,1] = [-0.1  0.6  -0.1]
```
Which matches pretty well.

To verify the additive adjustment above I looked at the one step transition using the function `one_step_pred` from "Examples/Testers/testertool.jl"
```julia
m_pred_s = mean(pred_s,dims = 2)
S = 2 .- sum(h_wf_s[:,:,j]*m_pred_s for j = 1 : M_out)

X_hat_c = one_step_pred(sig_c, h_wf_c, pred_c)
X_hat_s = one_step_pred(sig_s, h_wf_s, pred_s) .+ S

wind = (1:100) .+ 24000
plot([sig_c[1,wind] X_hat_c[1,wind.+1]])
plot([sig_s[1,wind] X_hat_s[1,wind.+1]])
```
The predictors and the original match in both cases.


# Monday, November 30, 2020

12:52 PM - The centered and noncentered data was good to investigating. Now, I will get back to computing the Wiener filters with better parameters. I return to thelio and the notebook from Tuesday, Nov 24, 2020.

Looking at the autocorrelations of the predictor series, it looks like they decay to roughly a margin of error after about 150 units of time, so this would be 1500 time steps, which means `par = 1500` should be enough to capture the correlations in the data. I will also try a little less than this and see how it goes.

Now we want to get a picture of the necessary resolution in frequency space. For this I look at the cross covariances. It seemed to me that there was very weak (if any) correlation between any pair of distinct predictors excepts the identity terms (1-5) and the inviscid Burger terms (6-10). That is to say I might suspect the cross covariance (`sig` and `pred`) could look like this:
```
⋆ 0 0 0 0 ⋆ 0 0 0 0
0 ⋆ 0 0 0 0 ⋆ 0 0 0
0 0 ⋆ 0 0 0 0 ⋆ 0 0           zeros(5,15)
0 0 0 ⋆ 0 0 0 0 ⋆ 0
0 0 0 0 ⋆ 0 0 0 0 ⋆
```
where the stars an nonzero, the covariance of the predictors seems to have this structure
```
⋆ 0 0 0 0 ⋆ 0 0 0 0
0 ⋆ 0 0 0 0 ⋆ 0 0 0
0 0 ⋆ 0 0 0 0 ⋆ 0 0           
0 0 0 ⋆ 0 0 0 0 ⋆ 0
0 0 0 0 ⋆ 0 0 0 0 ⋆           zeros(10,15)
⋆ 0 0 0 0 ⋆ 0 0 0 0
0 ⋆ 0 0 0 0 ⋆ 0 0 0
0 0 ⋆ 0 0 0 0 ⋆ 0 0           
0 0 0 ⋆ 0 0 0 0 ⋆ 0
0 0 0 0 ⋆ 0 0 0 0 ⋆

zeros(15,10)                  ⋆-diagonal(15,15)
```

I don't quite know what to do with this structure yet. But I might as well experiment on my purposed values of `par = 1500` I will first use


#### Experiment Nov 30, 2020 1 (Langevin data with smallest possible `nfft`)

This experiment is a run of the double-welled overdamped Langevin process:

| parameter | value |
|--- |--- |
|`sigma` | 0.3 |
|`V_prime` | x -> -x.*(x.^2 .- 1) |
|`sig_init`| 1.5 |
|`scheme` | "FE" |
|`steps` | 10^7 |
|`h` (time step)| 0.1 |
|`discard` | 'steps'|
|`gap` | 1|

The Wiener filtering parameters are as follows:
```julia
# Put in Psi functions
Psi(x) = [x; x.^3]

# Model reduction Parameters
M_out = 20
ty = "bin"

### Varing parameters
###       xspect_est , par    , nfft    , n    , p
Parms = [["DM"       , 50     , 2^17    , 2    , 5],
         ["DM"       , 50     , 120     , 2    , 5],
         ["DM"       , 5000   , 2^14    , 2    , 5]]
```

The Weiner filtering code is run three times. Each with it's own associated row of the above `Parms` variables.

here are the results:
```
julia> times
3-element Array{Float64,1}:
 29.8823148
 30.949746601
 43.4594022

julia> h_wf_packs[1][1,:,1:10]'
10×2 LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},2}}:
       1.09824+2.64328e-14im     -0.099226-7.93992e-15im
    0.00180432-3.31373e-14im  -0.000795544+1.03531e-14im
  -0.000809017-8.12924e-15im    6.95014e-5+2.60707e-15im
    0.00118388-4.62527e-16im   -0.00021893-1.16555e-16im
  -0.000236858+9.45978e-15im    4.61908e-5-2.87852e-15im
  -0.000561057+6.48942e-15im   0.000165708-2.3882e-15im
   0.000879032+3.28968e-14im  -0.000350527-9.49214e-15im
    0.00110459-5.01987e-14im  -0.000240615+1.52815e-14im
   -0.00213206+3.43862e-14im   0.000517865-1.08814e-14im
    0.00210348-4.83778e-14im  -0.000620188+1.51543e-14im

julia> h_wf_packs[8][1,:,1:10]'
10×2 LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},2}}:
      1.09748+5.69274e-14im    -0.0989922-1.6957e-14im
   0.00197576-1.32301e-13im  -0.000842647+4.03524e-14im
 -0.000734024+1.0566e-13im     4.83054e-5-3.17053e-14im
    0.0011944-4.93628e-14im  -0.000222027+1.50126e-14im
  -0.00025543-8.21179e-15im    5.11942e-5+1.61224e-15im
 -0.000577133+7.15856e-14im   0.000169312-2.20492e-14im
  0.000890196-1.77782e-14im  -0.000355921+6.84521e-15im
   0.00116187-9.08725e-14im  -0.000260733+2.67209e-14im
  -0.00201997+9.09829e-14im   0.000480553-2.75675e-14im
   0.00226817-3.75216e-14im  -0.000673646+1.1115e-14im

julia> h_wf_packs[15][1,:,1:10]'
10×2 LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},2}}:
     1.09698+3.44339e-14im    -0.0988348-1.03353e-14im
  0.00183088-3.28098e-14im  -0.000800667+9.90803e-15im
-0.000678346-1.24045e-14im    2.57955e-5+4.46203e-15im
  0.00154754-2.70635e-14im  -0.000333177+8.09354e-15im
 0.000106247+4.18996e-14im   -6.75428e-5-1.28916e-14im
-0.000196658-4.25912e-14im    4.28827e-5+1.26922e-14im
  0.00143237+1.07269e-13im  -0.000537704-3.24095e-14im
  0.00171027-1.09667e-13im  -0.000422434+3.3529e-14im
 -0.00220578+5.65421e-14im   0.000508484-1.76376e-14im
  0.00277172-7.97258e-14im  -0.000828119+2.48918e-14im
```

# Tuesday, December 1, 2020

12:57 PM - I just most of a paper from 1987 on constrained wiener filtering.

Picinbono, Bernard, and Michel Bouvet. "Constrained Wiener filtering (Corresp.)." *IEEE transactions on information theory* 33.1 (1987): 160-166.

according to Google scholar the article has been cited 14 times. So, maybe that is not so good.

I recently became interested in constrained Wiener filtering after having observed that many of the variables seem that they might be uncorrelated. However this does not amount, thinking about it now to a constraint on the Wiener filter itself.

1:32 PM - Today I want to run the Wiener filter with my new purposed values. This will be the experiment below.

#### Experiment Dec 1, 2020 1 (thelio job 161)

This is a variation of *Experiment Nov 24, 2020 1*. I will use the same script ("Examples/KSE/KSE_wf.jl") but with the following Wiener filtering parameters.
```julia
# Wiener filtering parameters
M_out = 20
nfft = 2^14
par = 1500
xspec_est = "old" # Default
short = true
loadsol = true
```
The experiment was run of thelio. Here is the resulting WF:
```
julia> h_wf = load("KSE_wf12_01_20_1-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  755.809+16.8891im  -7.99654-19.2875im    1.74191+9.5422im   -23.5897-12.8582im  -119.178+1.76331im  …   0.000137812+3.1384e-5im    -2.88503e-5+7.90048e-5im    0.000229118+3.43086e-5im
 -227.292-16.5817im   795.051+26.4794im   -74.0062-52.6848im  -70.0933+19.5675im  -132.655-4.29833im     -0.000407065-0.000309832im  0.000186247-0.000391773im   -0.00263106+0.00133284im
 -37.5647-24.1439im  -139.802+15.4142im    805.113-10.1535im  -199.727-24.9106im   268.982+47.0386im       0.00076189-0.000231712im  0.000459718+0.000373766im   -0.00137725-0.00182692im
 -77.6234+7.29485im  -50.4557+26.1065im   -57.7761-35.3984im   306.858+14.3475im  -133.931+2.2222im         5.9025e-5+7.73847e-5im    7.53998e-5-0.000218475im  -0.000797737+0.00113539im
 -18.0482+17.2948im   2.66778+0.449046im  -3.84975-1.16959im  -2.38131+6.4214im    19.0767-7.70804im      -3.30878e-5+1.52914e-5im    7.32968e-5-3.03003e-5im     2.64164e-5-0.0002416im

[:, :, 2] =
 -1569.82-38.3816im   23.3827+41.5954im  -3.65188-17.6442im   48.4977+25.2517im   265.503-5.60642im  …  -0.000381726-9.851e-5im       7.74235e-5-0.000209311im  -0.000655757-0.0001421im
  506.821+30.7912im  -1648.16-61.4394im   170.721+116.493im   140.971-43.905im    294.836+14.3094im       0.00117041+0.000912478im  -0.000457616+0.000948904im    0.00745586-0.00334468im
  86.1736+74.0488im   334.109-42.6377im  -1701.39+30.1481im   458.821+55.1135im   -612.09-114.012im      -0.00214064+0.000603504im   -0.00120283-0.000982147im    0.00385627+0.00521701im
  152.451-23.3338im   103.221-56.9231im   120.205+79.0413im  -578.024-35.6993im   278.444-5.32853im     -0.000166412-0.000228459im  -0.000181711+0.000579489im    0.00218827-0.00303165im
  37.5742-37.2256im  -5.77793-2.23601im   8.58938+2.92126im   3.60399-13.397im   -19.8948+16.2505im       9.48261e-5-3.91243e-5im    -0.00019237+6.85788e-5im    -5.88643e-5+0.000668121im

[:, :, 3] =
  719.315+19.3151im  -16.0441-20.9664im   1.87975+6.02222im   -20.7038-10.97im    -131.568+4.16373im  …   0.000286198+8.71628e-5im   -5.56506e-5+0.000146257im  0.000510602+0.000159712im
  -253.09-10.9498im   758.744+32.5812im  -89.7797-59.7251im   -55.6395+22.5034im  -145.417-10.1288im     -0.000919663-0.000743511im  0.000266214-0.000571507im  -0.00576746+0.00221613im
 -43.3236-53.5248im  -185.296+34.4566im   807.019-23.649im    -238.677-27.7674im   310.003+66.2367im       0.00165286-0.000424311im   0.00079691+0.000655438im  -0.00296839-0.00411608im
 -61.0435+16.4077im  -42.2946+27.6957im  -53.7297-39.6852im    233.151+20.123im   -122.747+3.82663im      0.000132639+0.00018845im   0.000103775-0.000408917im  -0.00160325+0.00219919im
  -16.631+17.7166im   2.67324+2.28866im  -4.53776-1.7623im   -0.547991+6.03799im   1.84314-7.66763im      -7.44482e-5+2.58022e-5im   0.000130867-3.77832e-5im    2.38248e-5-0.000505687im

...

[:, :, 18] =
   10.2358+0.531073im  -0.0928337-0.203508im   0.0841993+0.307561im   -0.494304-0.0132043im     -1.9992-0.0375501im  …  -1.18453e-6+6.2722e-8im    6.07048e-7-1.87745e-6im  -1.64603e-6+1.6099e-6im
  -3.24454-0.265899im     10.4284+0.614137im   -0.635679-0.206811im    -2.14139+0.295488im     -2.34184-0.196591im       3.54675e-6-3.37443e-8im  -7.33874e-6+8.81498e-6im   1.30385e-5-1.15063e-5im
 0.0939445+0.755843im    -1.21206-0.405799im     8.70961+0.742676im    -3.18528-0.288548im      5.06877+0.348294im      -1.48653e-6+1.31741e-6im  -1.41961e-5-1.70999e-5im   1.03116e-6+1.62778e-6im
  -1.11486+0.206345im    -1.38571+0.320466im   -0.748577-0.366684im     2.97784+0.436729im     -2.36271-0.0786643im      2.86315e-7-1.95628e-7im  -2.76108e-6+4.00493e-6im   1.03062e-5-3.29003e-6im
 -0.420108+0.259633im   0.0374631-0.0255983im  -0.025495+0.0282555im  -0.013918+0.0824061im  -0.0558794-0.149825im       2.41567e-7-2.5201e-7im   -1.84912e-6+2.09734e-7im  -1.82017e-6+3.37703e-7im

[:, :, 19] =
   5.62406+0.266855im    0.102505-0.155868im       0.156626+0.309237im    -0.416107-0.0624584im  -0.821485-0.241043im    …   1.79579e-6+1.40979e-6im  -6.84511e-8-8.77884e-7im    3.7089e-6+3.99403e-6im
  -1.18973-0.201809im     6.78185+0.107483im       0.557635+0.142722im     -1.81839+0.0359051im   -1.36928+0.0470964im      -4.35543e-6-9.92192e-6im  -8.35897e-6+8.25016e-6im  -4.83225e-5-1.56527e-6im
  0.461106+1.3594im        1.0254+0.192592im        3.32212+0.844007im    -0.471182+0.0867182im    1.35972-0.328128im        1.47953e-5-1.14852e-6im  -1.23547e-5-1.44486e-5im  -1.76507e-5-4.44814e-5im
  -1.12549+0.177664im   -0.854812+0.0518135im     -0.354863+0.0286334im     2.53934+0.280216im    -2.15298+0.00644204im      8.16351e-7+1.65752e-6im  -2.98549e-6+1.56807e-7im  -1.30641e-5+1.41783e-5im
 -0.309673+0.127419im  -0.0185458-0.00746015im  -0.00633013+0.0525899im  -0.0692969+0.0651832im  0.0182528-0.115419im       -5.97111e-7-1.24233e-7im  -1.01583e-6+2.3266e-7im   -5.79772e-8-4.35809e-6im

[:, :, 20] =
  0.748594-0.033087im      0.259389-0.0927239im   0.184712+0.244339im    -0.28087-0.092235im    0.302562-0.357833im   …   3.39889e-6+2.04469e-6im  -6.09092e-7+1.47444e-7im   6.20729e-6+4.70983e-6im
  0.800133-0.0942719im      2.45087-0.306475im     1.43901+0.359878im    -1.21679-0.174894im   -0.317788+0.204138im      -8.70227e-6-1.42877e-5im  -6.61357e-6+5.23074e-6im  -7.80137e-5+4.98976e-6im
  0.764098+1.61388im        2.72878+0.526509im     -1.6516+0.745632im     1.85887+0.373472im    -1.98783-0.809004im       2.18631e-5-2.29975e-6im   -7.2806e-6-7.71774e-6im  -2.41591e-5-6.28942e-5im
 -0.884206+0.0921793im    -0.352349-0.135436im   0.0589821+0.338333im      1.7194+0.098785im    -1.60485+0.0654073im      1.16684e-6+2.34122e-6im  -2.20493e-6-2.90078e-6im  -2.56124e-5+2.1326e-5im
 -0.148989+0.00313676im  -0.0685945+0.0028647im  0.0137244+0.0588467im  -0.102036+0.0421728im  0.0679917-0.0643922im     -1.00699e-6+2.09907e-8im  -5.28237e-8+1.83839e-7im   1.55833e-6-6.24226e-6im
```
Here are more output results:
```
Get_wf computation time: 494.230526 seconds (435.98 M allocations: 954.198 GiB, 2.22% gc time)

CKMS Computation time: 243.893564 seconds (169.32 M allocations: 864.214 GiB, 2.96% gc time)
Number of CKMS iterations: 6820
errK errR : 9.848219888213858e-11 3.2746787123590665e-15
```

#### Experiment Dec 1, 2020 2 (thelio job 162)

This experiment is a variation of the one before it (*Experiment Dec 1, 2020 1*), but with
```julia
nfft = 2^12
par = 500
```
The experiment was run of thelio. Here is the resulting WF:
```
h_wf = load("KSE_wf12_01_20_2-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
 203.584+3.09534im  -2.73014-5.95842im    0.0616854+2.07933im   -4.25374-3.2134im   -27.7598+0.916238im  …    3.54291e-5+8.48023e-6im  -7.62406e-6+1.64215e-5im   0.000112775+2.32236e-5im
-42.8952-7.91818im   201.806+5.01733im     -21.5038-14.4517im   -12.7938+4.05026im  -23.8118+2.69453im      -0.000160175-6.78447e-5im   4.13392e-5-8.30701e-5im  -0.000855814+0.000338669im
-4.46939-1.53189im  -46.1351+4.31942im      230.519-0.894398im  -56.4492-3.75571im    51.892+11.1621im       0.000205335-9.20794e-5im  0.000153782+8.43695e-5im  -0.000503379-0.000467105im
-18.8168-4.67129im   -12.523+5.21867im     -18.6121-8.85055im     93.564+2.37447im  -29.5831+0.571702im       2.46067e-5+1.66097e-5im    2.6309e-5-5.43579e-5im  -0.000178782+0.000304895im
-7.54098+3.15272im  0.167248-0.0399667im  -0.422169-0.138953im  -1.37431+1.25712im   15.4743-1.82322im       -7.47703e-6+9.26569e-6im   1.48772e-5-9.62472e-6im    2.90185e-5-5.24079e-5im

[:, :, 2] =
-335.014-5.54885im    6.25747+10.47im     -1.12689-2.98101im    7.85653+5.08438im   54.3551-1.97964im  …   -7.78269e-5-2.28751e-5im     1.70091e-5-3.55928e-5im   -0.000263477-7.29107e-5im
 81.7936+14.1785im   -322.121-7.29737im    38.0402+23.4197im    21.1328-7.9609im    47.4065-2.59343im      0.000407157+0.000171595im   -7.27576e-5+0.000160077im    0.00204946-0.000606909im
 5.50337+7.4361im     85.6178-9.33099im   -376.445+2.6297im     105.401+6.75082im  -111.104-21.1647im     -0.000484024+0.000195714im  -0.000336198-0.000181606im    0.00111332+0.00107764im
 33.4921+6.20733im    21.5971-9.64479im     34.336+16.3788im   -133.735-4.83519im   57.0795-1.4555im       -5.82902e-5-4.01614e-5im    -5.28118e-5+0.000112765im   0.000371606-0.000638284im
 13.6716-5.7004im   0.0437112-0.165886im  0.239548+0.268183im   2.56193-2.29049im  -15.5816+3.05984im       1.83619e-5-1.97876e-5im    -3.06454e-5+1.74626e-5im    -5.98895e-5+0.000109933im

[:, :, 3] =
 109.178+2.16627im    -3.02891-3.68523im    1.38635+0.587794im   -2.99669-1.60396im  -22.5612+1.01748im   …    3.96623e-5+1.59496e-5im   -9.26017e-6+1.79888e-5im   0.000148154+5.81522e-5im
-32.4015-4.80239im     94.8419+1.50666im   -11.3306-6.56923im    -5.28268+3.6158im   -19.9337-0.700065im     -0.000267965-0.000109666im   1.85941e-5-6.87889e-5im   -0.00122774+0.000200828im
 2.33927-5.95766im    -29.7785+5.44549im    116.625-1.90154im     -39.677-2.7592im    53.7026+8.60292im       0.000287575-9.57117e-5im   0.000168958+8.78634e-5im  -0.000564711-0.000592969im
-12.1831-0.999153im   -6.49443+3.88604im   -12.8824-6.23348im     36.4343+2.15856im  -23.0858+1.07689im        3.62753e-5+2.40499e-5im    2.31514e-5-5.38497e-5im  -0.000163579+0.000315496im
-5.09631+2.12097im   -0.315438+0.205453im  0.298656-0.105687im  -0.953936+0.83574im   2.34589-1.03835im       -1.12143e-5+9.74732e-6im    1.38752e-5-6.80201e-6im    2.65997e-5-5.27751e-5im

...

[:, :, 18] =
 -0.676249-0.0177268im  -0.0390984-0.0115825im   0.00113469-0.014869im     0.0116395+0.0242489im   -0.0996344+0.0240919im    …  -1.56936e-7-2.31869e-7im  2.76606e-8+7.84498e-8im  -9.70673e-7-1.09546e-6im
 -0.115476-0.0265847im   -0.777687-0.0926653im    -0.192549-0.0314538im    0.0882896-0.00558733im   -0.061598-0.0824794im        1.91187e-6+1.58777e-6im  1.32606e-6-6.10967e-7im   9.40108e-6+3.17627e-6im
0.00497367-0.173377im    -0.497138-0.0706967im    -0.432329-0.139627im    -0.0354698-0.0390218im     0.508442+0.113901im        -2.33989e-6+3.14138e-7im  1.87559e-6+1.87429e-6im   2.75776e-6+6.09233e-6im
 0.0790026+0.0665533im   0.0188128+0.0257761im    0.0160278-0.0266972im    -0.534553-0.00398996im   0.0140033+0.00775984im      -2.18616e-7-2.21576e-7im  4.32165e-7+1.98055e-8im   6.35932e-7-1.11434e-6im
0.00291016+0.0015479im  0.00208272+0.00101381im   0.0166983-0.00112902im   0.0116093-0.00226967im  -0.0939914-0.000645237im      8.30002e-8+1.38662e-8im  1.24766e-7-5.79864e-8im  -6.51902e-8+2.0183e-7im

[:, :, 19] =
 -0.107326-0.00354825im  -0.0529434-0.0275677im    -0.0147906-0.00291988im   -0.001733+0.0218804im    …  -3.79504e-7-2.83202e-7im  9.47979e-8-4.54791e-8im  -1.50251e-6-1.10328e-6im
  -0.26379-0.0617491im   -0.0754421-0.0267948im     -0.314637-0.0713229im   -0.0181322+0.00547635im       2.03747e-6+1.84545e-6im  1.09925e-6-4.02788e-7im   1.23223e-5+8.02536e-7im
 -0.130149-0.142068im     -0.716405-0.106423im       0.334922-0.0960096im    -0.278958-0.0593884im       -2.98152e-6+7.29125e-7im  9.98104e-7+1.09905e-6im   5.24631e-6+9.15017e-6im
 0.0306405+0.0388469im   -0.0643457+0.0274887im    -0.0445321-0.0583908im    -0.458963-0.00843077im      -2.66355e-7-3.0624e-7im    3.2179e-7+3.78091e-7im   2.57665e-6-2.36557e-6im
-0.0231621+0.0102394im   0.00476338-0.000219589im   0.0125938-0.00110315im  0.00742414+0.000189929im      1.26903e-7-4.26637e-8im  2.72867e-8-4.24455e-8im  -3.78451e-7+4.57819e-7im

[:, :, 20] =
  0.353687+0.00877631im  -0.0584751-0.0371183im    -0.025772+0.0072875im    -0.0124737+0.0171149im    -0.22694+0.0164599im  …  -4.47823e-7-2.51148e-7im  1.29967e-7-1.2875e-7im   -1.49052e-6-8.19747e-7im
 -0.351059-0.0805261im     0.494024+0.0322071im    -0.366634-0.0931863im      -0.10239+0.0137008im   -0.209124-0.0586322im      1.56954e-6+1.54308e-6im  6.76826e-7-1.57354e-7im   1.11851e-5-1.2214e-6im
 -0.222554-0.0957266im    -0.796024-0.116031im      0.893118-0.0416644im     -0.442046-0.0672578im    0.525374+0.0962727im      -2.6615e-6+8.40181e-7im  1.06708e-7+2.08138e-7im   5.65056e-6+9.04069e-6im
 -0.011304+0.0123095im    -0.123663+0.0244327im   -0.0858703-0.0754044im     -0.355745-0.0121001im   -0.151357-0.0224742im     -2.32593e-7-2.911e-7im    1.60499e-7+5.62768e-7im   3.39388e-6-2.60752e-6im
-0.0417339+0.0159985im   0.00568693-0.00112026im  0.00736822-0.000760594im  0.00302293+0.00200042im  -0.090206-0.0106053im      1.28267e-7-7.36953e-8im  -5.3469e-8-2.35215e-8im  -5.24846e-7+5.16608e-7im
```
Here are more output results:
```
Get_wf computation time: 276.235305 seconds (321.11 M allocations: 214.306 GiB, 1.42% gc time)

CKMS Computation time: 34.419597 seconds (58.60 M allocations: 125.967 GiB, 1.72% gc time)
Number of CKMS iterations: 2950
errK errR :  6.883824319059169e-11 2.153244140484029e-14
```

#### Experiment Dec 1, 2020 3 (thelio job 163)

This experiment is a variation of the one before it (*Experiment Dec 1, 2020 1*), but with
```julia
nfft = 2^14
par = 500
```
The experiment was run of thelio. Here is the resulting WF:
```
julia> h_wf = load("KSE_wf12_01_20_3-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  203.584+3.09534im  -2.73014-5.95842im    0.0616854+2.07933im   -4.25374-3.2134im   -27.7598+0.916238im  …    3.54291e-5+8.48023e-6im  -7.62406e-6+1.64215e-5im   0.000112775+2.32236e-5im
 -42.8952-7.91818im   201.806+5.01733im     -21.5038-14.4517im   -12.7938+4.05026im  -23.8118+2.69453im      -0.000160175-6.78447e-5im   4.13392e-5-8.30701e-5im  -0.000855814+0.000338669im
 -4.46939-1.53189im  -46.1351+4.31942im      230.519-0.894398im  -56.4492-3.75571im    51.892+11.1621im       0.000205335-9.20794e-5im  0.000153782+8.43695e-5im  -0.000503379-0.000467105im
 -18.8168-4.67129im   -12.523+5.21867im     -18.6121-8.85055im     93.564+2.37447im  -29.5831+0.571702im       2.46067e-5+1.66097e-5im    2.6309e-5-5.43579e-5im  -0.000178782+0.000304895im
 -7.54098+3.15272im  0.167248-0.0399667im  -0.422169-0.138953im  -1.37431+1.25712im   15.4743-1.82322im       -7.47703e-6+9.26569e-6im   1.48772e-5-9.62472e-6im    2.90185e-5-5.24079e-5im

[:, :, 2] =
 -335.014-5.54885im    6.25747+10.47im     -1.12689-2.98101im    7.85653+5.08438im   54.3551-1.97964im  …   -7.78269e-5-2.28751e-5im     1.70091e-5-3.55928e-5im   -0.000263477-7.29107e-5im
  81.7936+14.1785im   -322.121-7.29737im    38.0402+23.4197im    21.1328-7.9609im    47.4065-2.59343im      0.000407157+0.000171595im   -7.27576e-5+0.000160077im    0.00204946-0.000606909im
  5.50337+7.4361im     85.6178-9.33099im   -376.445+2.6297im     105.401+6.75082im  -111.104-21.1647im     -0.000484024+0.000195714im  -0.000336198-0.000181606im    0.00111332+0.00107764im
  33.4921+6.20733im    21.5971-9.64479im     34.336+16.3788im   -133.735-4.83519im   57.0795-1.4555im       -5.82902e-5-4.01614e-5im    -5.28118e-5+0.000112765im   0.000371606-0.000638284im
  13.6716-5.7004im   0.0437112-0.165886im  0.239548+0.268183im   2.56193-2.29049im  -15.5816+3.05984im       1.83619e-5-1.97876e-5im    -3.06454e-5+1.74626e-5im    -5.98895e-5+0.000109933im

[:, :, 3] =
  109.178+2.16627im    -3.02891-3.68523im    1.38635+0.587794im   -2.99669-1.60396im  -22.5612+1.01748im   …    3.96623e-5+1.59496e-5im   -9.26017e-6+1.79888e-5im   0.000148154+5.81522e-5im
 -32.4015-4.80239im     94.8419+1.50666im   -11.3306-6.56923im    -5.28268+3.6158im   -19.9337-0.700065im     -0.000267965-0.000109666im   1.85941e-5-6.87889e-5im   -0.00122774+0.000200828im
  2.33927-5.95766im    -29.7785+5.44549im    116.625-1.90154im     -39.677-2.7592im    53.7026+8.60292im       0.000287575-9.57117e-5im   0.000168958+8.78634e-5im  -0.000564711-0.000592969im
 -12.1831-0.999153im   -6.49443+3.88604im   -12.8824-6.23348im     36.4343+2.15856im  -23.0858+1.07689im        3.62753e-5+2.40499e-5im    2.31514e-5-5.38497e-5im  -0.000163579+0.000315496im
 -5.09631+2.12097im   -0.315438+0.205453im  0.298656-0.105687im  -0.953936+0.83574im   2.34589-1.03835im       -1.12143e-5+9.74732e-6im    1.38752e-5-6.80201e-6im    2.65997e-5-5.27751e-5im

...

[:, :, 18] =
  -0.676249-0.0177268im  -0.0390984-0.0115825im   0.00113469-0.014869im     0.0116395+0.0242489im   -0.0996344+0.0240919im    …  -1.56936e-7-2.31869e-7im  2.76606e-8+7.84498e-8im  -9.70673e-7-1.09546e-6im
  -0.115476-0.0265847im   -0.777687-0.0926653im    -0.192549-0.0314538im    0.0882896-0.00558733im   -0.061598-0.0824794im        1.91187e-6+1.58777e-6im  1.32606e-6-6.10967e-7im   9.40108e-6+3.17627e-6im
 0.00497367-0.173377im    -0.497138-0.0706967im    -0.432329-0.139627im    -0.0354698-0.0390218im     0.508442+0.113901im        -2.33989e-6+3.14138e-7im  1.87559e-6+1.87429e-6im   2.75776e-6+6.09233e-6im
  0.0790026+0.0665533im   0.0188128+0.0257761im    0.0160278-0.0266972im    -0.534553-0.00398996im   0.0140033+0.00775984im      -2.18616e-7-2.21576e-7im  4.32165e-7+1.98055e-8im   6.35932e-7-1.11434e-6im
 0.00291016+0.0015479im  0.00208272+0.00101381im   0.0166983-0.00112902im   0.0116093-0.00226967im  -0.0939914-0.000645239im      8.30002e-8+1.38662e-8im  1.24766e-7-5.79864e-8im  -6.51902e-8+2.0183e-7im

[:, :, 19] =
  -0.107326-0.00354825im  -0.0529434-0.0275677im    -0.0147906-0.00291988im   -0.001733+0.0218804im    …  -3.79504e-7-2.83202e-7im  9.47979e-8-4.54791e-8im  -1.50251e-6-1.10328e-6im
   -0.26379-0.0617491im   -0.0754421-0.0267948im     -0.314637-0.0713229im   -0.0181322+0.00547635im       2.03747e-6+1.84545e-6im  1.09925e-6-4.02788e-7im   1.23223e-5+8.02536e-7im
  -0.130149-0.142068im     -0.716405-0.106423im       0.334922-0.0960096im    -0.278958-0.0593884im       -2.98152e-6+7.29125e-7im  9.98104e-7+1.09905e-6im   5.24631e-6+9.15017e-6im
  0.0306405+0.0388469im   -0.0643457+0.0274887im    -0.0445321-0.0583908im    -0.458963-0.00843077im      -2.66355e-7-3.0624e-7im    3.2179e-7+3.78091e-7im   2.57665e-6-2.36557e-6im
 -0.0231621+0.0102394im   0.00476338-0.000219589im   0.0125938-0.00110315im  0.00742414+0.000189929im      1.26903e-7-4.26637e-8im  2.72867e-8-4.24455e-8im  -3.78451e-7+4.57819e-7im

[:, :, 20] =
   0.353687+0.00877631im  -0.0584751-0.0371183im    -0.025772+0.0072875im    -0.0124737+0.0171149im    -0.22694+0.0164599im  …  -4.47823e-7-2.51148e-7im  1.29967e-7-1.2875e-7im   -1.49052e-6-8.19747e-7im
  -0.351059-0.0805261im     0.494024+0.0322071im    -0.366634-0.0931863im      -0.10239+0.0137008im   -0.209124-0.0586322im      1.56954e-6+1.54308e-6im  6.76826e-7-1.57354e-7im   1.11851e-5-1.2214e-6im
  -0.222554-0.0957266im    -0.796024-0.116031im      0.893118-0.0416644im     -0.442046-0.0672578im    0.525374+0.0962727im      -2.6615e-6+8.40181e-7im  1.06708e-7+2.08138e-7im   5.65056e-6+9.04069e-6im
  -0.011304+0.0123095im    -0.123663+0.0244327im   -0.0858703-0.0754044im     -0.355745-0.0121001im   -0.151357-0.0224742im     -2.32593e-7-2.911e-7im    1.60499e-7+5.62768e-7im   3.39388e-6-2.60752e-6im
 -0.0417339+0.0159985im   0.00568693-0.00112026im  0.00736822-0.000760594im  0.00302293+0.00200042im  -0.090206-0.0106053im      1.28267e-7-7.36953e-8im  -5.3469e-8-2.35215e-8im  -5.24846e-7+5.16608e-7im
```
Here are more output results:
```
Get_wf computation time: 281.941092 seconds (321.62 M allocations: 215.771 GiB, 1.46% gc time)

CKMS Computation time: 34.760786 seconds (58.60 M allocations: 125.967 GiB, 1.89% gc time)
Number of CKMS iterations: 2950
errK errR :  6.883824319059169e-11 2.153244140484029e-14 (same as Exper 12/01/20 2)
```

3:30 PM - I would like to have another run of the data to test, just to check that this is pretty insensitive to the run itself. So I will run "Examples/KSE/KSE_data_gen.jl" again this time with `seed = 2021` as opposed to `seed = 2020` as before (see Nov 24, 2020, 10:31 AM). The data was save as "data/KSE_Data/KSE_sol_lin1.jld". Here, `gen = "lin1"`. This was Thelio job 165.


#### Experiment Dec 1, 2020 4 (thelio job 166) (Exper 12/01/20 1 on different data)

This experiment is exactly the same as *Experiment Dec 1, 2020 1* except that is runs on the data from "data/KSE_Data/KSE_sol_lin1.jld" (`seed = 2021`).

The experiment was run of thelio. Here is the resulting WF:
```
h_wf = load("KSE_wf12_01_20_4-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  755.809+16.8891im  -7.99654-19.2875im    1.74191+9.5422im   -23.5897-12.8582im  -119.178+1.76331im  …   0.000137812+3.1384e-5im    -2.88503e-5+7.90048e-5im    0.000229118+3.43086e-5im
 -227.292-16.5817im   795.051+26.4794im   -74.0062-52.6848im  -70.0933+19.5675im  -132.655-4.29833im     -0.000407065-0.000309832im  0.000186247-0.000391773im   -0.00263106+0.00133284im
 -37.5647-24.1439im  -139.802+15.4142im    805.113-10.1535im  -199.727-24.9106im   268.982+47.0386im       0.00076189-0.000231712im  0.000459718+0.000373766im   -0.00137725-0.00182692im
 -77.6234+7.29485im  -50.4557+26.1065im   -57.7761-35.3984im   306.858+14.3475im  -133.931+2.2222im         5.9025e-5+7.73847e-5im    7.53998e-5-0.000218475im  -0.000797737+0.00113539im
 -18.0482+17.2948im   2.66778+0.449046im  -3.84975-1.16959im  -2.38131+6.4214im    19.0767-7.70804im      -3.30878e-5+1.52914e-5im    7.32968e-5-3.03003e-5im     2.64164e-5-0.0002416im

[:, :, 2] =
 -1569.82-38.3816im   23.3827+41.5954im  -3.65188-17.6442im   48.4977+25.2517im   265.503-5.60642im  …  -0.000381726-9.851e-5im       7.74235e-5-0.000209311im  -0.000655757-0.0001421im
  506.821+30.7912im  -1648.16-61.4394im   170.721+116.493im   140.971-43.905im    294.836+14.3094im       0.00117041+0.000912478im  -0.000457616+0.000948904im    0.00745586-0.00334468im
  86.1736+74.0488im   334.109-42.6377im  -1701.39+30.1481im   458.821+55.1135im   -612.09-114.012im      -0.00214064+0.000603504im   -0.00120283-0.000982147im    0.00385627+0.00521701im
  152.451-23.3338im   103.221-56.9231im   120.205+79.0413im  -578.024-35.6993im   278.444-5.32853im     -0.000166412-0.000228459im  -0.000181711+0.000579489im    0.00218827-0.00303165im
  37.5742-37.2256im  -5.77793-2.23601im   8.58938+2.92126im   3.60399-13.397im   -19.8948+16.2505im       9.48261e-5-3.91243e-5im    -0.00019237+6.85788e-5im    -5.88643e-5+0.000668121im

[:, :, 3] =
  719.315+19.3151im  -16.0441-20.9664im   1.87975+6.02222im   -20.7038-10.97im    -131.568+4.16373im  …   0.000286198+8.71628e-5im   -5.56506e-5+0.000146257im  0.000510602+0.000159712im
  -253.09-10.9498im   758.744+32.5812im  -89.7797-59.7251im   -55.6395+22.5034im  -145.417-10.1288im     -0.000919663-0.000743511im  0.000266214-0.000571507im  -0.00576746+0.00221613im
 -43.3236-53.5248im  -185.296+34.4566im   807.019-23.649im    -238.677-27.7674im   310.003+66.2367im       0.00165286-0.000424311im   0.00079691+0.000655438im  -0.00296839-0.00411608im
 -61.0435+16.4077im  -42.2946+27.6957im  -53.7297-39.6852im    233.151+20.123im   -122.747+3.82663im      0.000132639+0.00018845im   0.000103775-0.000408917im  -0.00160325+0.00219919im
  -16.631+17.7166im   2.67324+2.28866im  -4.53776-1.7623im   -0.547991+6.03799im   1.84314-7.66763im      -7.44482e-5+2.58022e-5im   0.000130867-3.77832e-5im    2.38248e-5-0.000505687im

...

[:, :, 18] =
   10.2358+0.531073im  -0.0928337-0.203508im   0.0841993+0.307561im   -0.494304-0.0132043im     -1.9992-0.0375501im  …  -1.18453e-6+6.2722e-8im    6.07048e-7-1.87745e-6im  -1.64603e-6+1.6099e-6im
  -3.24454-0.265899im     10.4284+0.614137im   -0.635679-0.206811im    -2.14139+0.295488im     -2.34184-0.196591im       3.54675e-6-3.37443e-8im  -7.33874e-6+8.81498e-6im   1.30385e-5-1.15063e-5im
 0.0939445+0.755843im    -1.21206-0.405799im     8.70961+0.742676im    -3.18528-0.288548im      5.06877+0.348294im      -1.48653e-6+1.31741e-6im  -1.41961e-5-1.70999e-5im   1.03116e-6+1.62778e-6im
  -1.11486+0.206345im    -1.38571+0.320466im   -0.748577-0.366684im     2.97784+0.436729im     -2.36271-0.0786643im      2.86315e-7-1.95628e-7im  -2.76108e-6+4.00493e-6im   1.03062e-5-3.29003e-6im
 -0.420108+0.259633im   0.0374631-0.0255983im  -0.025495+0.0282555im  -0.013918+0.0824061im  -0.0558794-0.149825im       2.41567e-7-2.5201e-7im   -1.84912e-6+2.09734e-7im  -1.82017e-6+3.37703e-7im

[:, :, 19] =
   5.62406+0.266855im    0.102505-0.155868im       0.156626+0.309237im    -0.416107-0.0624584im  -0.821485-0.241043im    …   1.79579e-6+1.40979e-6im  -6.84511e-8-8.77884e-7im    3.7089e-6+3.99403e-6im
  -1.18973-0.201809im     6.78185+0.107483im       0.557635+0.142722im     -1.81839+0.0359051im   -1.36928+0.0470964im      -4.35543e-6-9.92192e-6im  -8.35897e-6+8.25016e-6im  -4.83225e-5-1.56527e-6im
  0.461106+1.3594im        1.0254+0.192592im        3.32212+0.844007im    -0.471182+0.0867182im    1.35972-0.328128im        1.47953e-5-1.14852e-6im  -1.23547e-5-1.44486e-5im  -1.76507e-5-4.44814e-5im
  -1.12549+0.177664im   -0.854812+0.0518135im     -0.354863+0.0286334im     2.53934+0.280216im    -2.15298+0.00644204im      8.16351e-7+1.65752e-6im  -2.98549e-6+1.56807e-7im  -1.30641e-5+1.41783e-5im
 -0.309673+0.127419im  -0.0185458-0.00746015im  -0.00633013+0.0525899im  -0.0692969+0.0651832im  0.0182528-0.115419im       -5.97111e-7-1.24233e-7im  -1.01583e-6+2.3266e-7im   -5.79772e-8-4.35809e-6im

[:, :, 20] =
  0.748594-0.033087im      0.259389-0.0927239im   0.184712+0.244339im    -0.28087-0.092235im    0.302562-0.357833im   …   3.39889e-6+2.04469e-6im  -6.09092e-7+1.47444e-7im   6.20729e-6+4.70983e-6im
  0.800133-0.0942719im      2.45087-0.306475im     1.43901+0.359878im    -1.21679-0.174894im   -0.317788+0.204138im      -8.70227e-6-1.42877e-5im  -6.61357e-6+5.23074e-6im  -7.80137e-5+4.98976e-6im
  0.764098+1.61388im        2.72878+0.526509im     -1.6516+0.745632im     1.85887+0.373472im    -1.98783-0.809004im       2.18631e-5-2.29975e-6im   -7.2806e-6-7.71774e-6im  -2.41591e-5-6.28942e-5im
 -0.884206+0.0921793im    -0.352349-0.135436im   0.0589821+0.338333im      1.7194+0.098785im    -1.60485+0.0654073im      1.16684e-6+2.34122e-6im  -2.20493e-6-2.90078e-6im  -2.56124e-5+2.1326e-5im
 -0.148989+0.00313676im  -0.0685945+0.0028647im  0.0137244+0.0588467im  -0.102036+0.0421728im  0.0679917-0.0643922im     -1.00699e-6+2.09907e-8im  -5.28237e-8+1.83839e-7im   1.55833e-6-6.24226e-6im
```
Here are more output results:
```
Get_wf computation time: 477.814471 seconds (435.98 M allocations: 954.198 GiB, 1.66% gc time)

CKMS Computation time: 225.586398 seconds (169.32 M allocations: 864.214 GiB, 1.77% gc time)
Number of CKMS iterations: 6820
errK errR :  9.848219888213858e-11 3.2746787123590665e-15
```


3:57 PM - I want to see what the WF function (`vector_wiener_filter_fft`) does with these uncorrelated modes. I will go back to the "KSE_data_analyzer" notebook and see what it does with these. We don't really get a zero cross spectrum of the modes that should be uncorrelated. I looked at a lot of the cross spectra and compared the direct estimator with the periodogram. They were very different for large values of `nfft` but pretty close for `nfft = 2^10`. This made me wonder if there would be an advantage to having `L` > `nfft` and I though no because that would no longer be smooth. there would be a truncation. I will look at this tomorrow it is now  little after 5 PM.

I think I would also like to rerun experiment 12/01/2020 2 and 3 tomorrow on different the other data.


# Wednesday, December 2, 2020

11:53 PM - To start off with I will go a head and run those other two experiments I wanted to run from yesterday. These are on the new KSE data this values of
```julia
nfft = 2^12
par = 500

nfft = 2^14
par = 500
```
respectively. I would also Like to try
```julia
nfft = 2^10
par = 500

nfft = 2^10
par = 100
```

Once I have these going I would like to look at the data again in the "KSE_data_analyzer" notebook.

#### Experiment Dec 2, 2020 1 (Exper 2 from yesterday on new data) (job 169)
This experiment is exactly the same as *Experiment Dec 1, 2020 2* except that is runs on the data from "data/KSE_Data/KSE_sol_lin1.jld" (`seed = 2021`). the key parameters are
```julia
nfft = 2^12
par = 500
```

The experiment was run of thelio. Here is the resulting WF:
```
julia> h_wf = load("KSE_wf12_02_20_1-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  203.584+3.09534im  -2.73014-5.95842im    0.0616854+2.07933im   -4.25374-3.2134im   -27.7598+0.916238im  …    3.54291e-5+8.48023e-6im  -7.62406e-6+1.64215e-5im   0.000112775+2.32236e-5im
 -42.8952-7.91818im   201.806+5.01733im     -21.5038-14.4517im   -12.7938+4.05026im  -23.8118+2.69453im      -0.000160175-6.78447e-5im   4.13392e-5-8.30701e-5im  -0.000855814+0.000338669im
 -4.46939-1.53189im  -46.1351+4.31942im      230.519-0.894398im  -56.4492-3.75571im    51.892+11.1621im       0.000205335-9.20794e-5im  0.000153782+8.43695e-5im  -0.000503379-0.000467105im
 -18.8168-4.67129im   -12.523+5.21867im     -18.6121-8.85055im     93.564+2.37447im  -29.5831+0.571702im       2.46067e-5+1.66097e-5im    2.6309e-5-5.43579e-5im  -0.000178782+0.000304895im
 -7.54098+3.15272im  0.167248-0.0399667im  -0.422169-0.138953im  -1.37431+1.25712im   15.4743-1.82322im       -7.47703e-6+9.26569e-6im   1.48772e-5-9.62472e-6im    2.90185e-5-5.24079e-5im

[:, :, 2] =
 -335.014-5.54885im    6.25747+10.47im     -1.12689-2.98101im    7.85653+5.08438im   54.3551-1.97964im  …   -7.78269e-5-2.28751e-5im     1.70091e-5-3.55928e-5im   -0.000263477-7.29107e-5im
  81.7936+14.1785im   -322.121-7.29737im    38.0402+23.4197im    21.1328-7.9609im    47.4065-2.59343im      0.000407157+0.000171595im   -7.27576e-5+0.000160077im    0.00204946-0.000606909im
  5.50337+7.4361im     85.6178-9.33099im   -376.445+2.6297im     105.401+6.75082im  -111.104-21.1647im     -0.000484024+0.000195714im  -0.000336198-0.000181606im    0.00111332+0.00107764im
  33.4921+6.20733im    21.5971-9.64479im     34.336+16.3788im   -133.735-4.83519im   57.0795-1.4555im       -5.82902e-5-4.01614e-5im    -5.28118e-5+0.000112765im   0.000371606-0.000638284im
  13.6716-5.7004im   0.0437112-0.165886im  0.239548+0.268183im   2.56193-2.29049im  -15.5816+3.05984im       1.83619e-5-1.97876e-5im    -3.06454e-5+1.74626e-5im    -5.98895e-5+0.000109933im

[:, :, 3] =
  109.178+2.16627im    -3.02891-3.68523im    1.38635+0.587794im   -2.99669-1.60396im  -22.5612+1.01748im   …    3.96623e-5+1.59496e-5im   -9.26017e-6+1.79888e-5im   0.000148154+5.81522e-5im
 -32.4015-4.80239im     94.8419+1.50666im   -11.3306-6.56923im    -5.28268+3.6158im   -19.9337-0.700065im     -0.000267965-0.000109666im   1.85941e-5-6.87889e-5im   -0.00122774+0.000200828im
  2.33927-5.95766im    -29.7785+5.44549im    116.625-1.90154im     -39.677-2.7592im    53.7026+8.60292im       0.000287575-9.57117e-5im   0.000168958+8.78634e-5im  -0.000564711-0.000592969im
 -12.1831-0.999153im   -6.49443+3.88604im   -12.8824-6.23348im     36.4343+2.15856im  -23.0858+1.07689im        3.62753e-5+2.40499e-5im    2.31514e-5-5.38497e-5im  -0.000163579+0.000315496im
 -5.09631+2.12097im   -0.315438+0.205453im  0.298656-0.105687im  -0.953936+0.83574im   2.34589-1.03835im       -1.12143e-5+9.74732e-6im    1.38752e-5-6.80201e-6im    2.65997e-5-5.27751e-5im

...

[:, :, 18] =
  -0.676249-0.0177268im  -0.0390984-0.0115825im   0.00113469-0.014869im     0.0116395+0.0242489im   -0.0996344+0.0240919im    …  -1.56936e-7-2.31869e-7im  2.76606e-8+7.84498e-8im  -9.70673e-7-1.09546e-6im
  -0.115476-0.0265847im   -0.777687-0.0926653im    -0.192549-0.0314538im    0.0882896-0.00558733im   -0.061598-0.0824794im        1.91187e-6+1.58777e-6im  1.32606e-6-6.10967e-7im   9.40108e-6+3.17627e-6im
 0.00497367-0.173377im    -0.497138-0.0706967im    -0.432329-0.139627im    -0.0354698-0.0390218im     0.508442+0.113901im        -2.33989e-6+3.14138e-7im  1.87559e-6+1.87429e-6im   2.75776e-6+6.09233e-6im
  0.0790026+0.0665533im   0.0188128+0.0257761im    0.0160278-0.0266972im    -0.534553-0.00398996im   0.0140033+0.00775984im      -2.18616e-7-2.21576e-7im  4.32165e-7+1.98055e-8im   6.35932e-7-1.11434e-6im
 0.00291016+0.0015479im  0.00208272+0.00101381im   0.0166983-0.00112902im   0.0116093-0.00226967im  -0.0939914-0.000645237im      8.30002e-8+1.38662e-8im  1.24766e-7-5.79864e-8im  -6.51902e-8+2.0183e-7im

[:, :, 19] =
  -0.107326-0.00354825im  -0.0529434-0.0275677im    -0.0147906-0.00291988im   -0.001733+0.0218804im    …  -3.79504e-7-2.83202e-7im  9.47979e-8-4.54791e-8im  -1.50251e-6-1.10328e-6im
   -0.26379-0.0617491im   -0.0754421-0.0267948im     -0.314637-0.0713229im   -0.0181322+0.00547635im       2.03747e-6+1.84545e-6im  1.09925e-6-4.02788e-7im   1.23223e-5+8.02536e-7im
  -0.130149-0.142068im     -0.716405-0.106423im       0.334922-0.0960096im    -0.278958-0.0593884im       -2.98152e-6+7.29125e-7im  9.98104e-7+1.09905e-6im   5.24631e-6+9.15017e-6im
  0.0306405+0.0388469im   -0.0643457+0.0274887im    -0.0445321-0.0583908im    -0.458963-0.00843077im      -2.66355e-7-3.0624e-7im    3.2179e-7+3.78091e-7im   2.57665e-6-2.36557e-6im
 -0.0231621+0.0102394im   0.00476338-0.000219589im   0.0125938-0.00110315im  0.00742414+0.000189929im      1.26903e-7-4.26637e-8im  2.72867e-8-4.24455e-8im  -3.78451e-7+4.57819e-7im

[:, :, 20] =
   0.353687+0.00877631im  -0.0584751-0.0371183im    -0.025772+0.0072875im    -0.0124737+0.0171149im    -0.22694+0.0164599im  …  -4.47823e-7-2.51148e-7im  1.29967e-7-1.2875e-7im   -1.49052e-6-8.19747e-7im
  -0.351059-0.0805261im     0.494024+0.0322071im    -0.366634-0.0931863im      -0.10239+0.0137008im   -0.209124-0.0586322im      1.56954e-6+1.54308e-6im  6.76826e-7-1.57354e-7im   1.11851e-5-1.2214e-6im
  -0.222554-0.0957266im    -0.796024-0.116031im      0.893118-0.0416644im     -0.442046-0.0672578im    0.525374+0.0962727im      -2.6615e-6+8.40181e-7im  1.06708e-7+2.08138e-7im   5.65056e-6+9.04069e-6im
  -0.011304+0.0123095im    -0.123663+0.0244327im   -0.0858703-0.0754044im     -0.355745-0.0121001im   -0.151357-0.0224742im     -2.32593e-7-2.911e-7im    1.60499e-7+5.62768e-7im   3.39388e-6-2.60752e-6im
 -0.0417339+0.0159985im   0.00568693-0.00112026im  0.00736822-0.000760594im  0.00302293+0.00200042im  -0.090206-0.0106053im      1.28267e-7-7.36953e-8im  -5.3469e-8-2.35215e-8im  -5.24846e-7+5.16608e-7im
```
Here are more output results:
```
Get_wf computation time: 274.970210 seconds (321.11 M allocations: 214.306 GiB, 1.46% gc time)

CKMS Computation time: 34.916413 seconds (58.60 M allocations: 125.967 GiB, 2.02% gc time)
Number of CKMS iterations: 2950
errK errR :  6.883824319059169e-11 2.153244140484029e-14
```
#### Experiment Dec 2, 2020 2 (Exper 3 from yesterday on new data) (job 170)
This experiment is exactly the same as *Experiment Dec 1, 2020 3* except that is runs on the data from "data/KSE_Data/KSE_sol_lin1.jld" (`seed = 2021`). the key parameters are
```julia
nfft = 2^14
par = 500
```

The experiment was run of thelio. Here is the resulting WF:
```
julia> h_wf = load("KSE_wf12_02_20_2-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  203.584+3.09534im  -2.73014-5.95842im    0.0616854+2.07933im   -4.25374-3.2134im   -27.7598+0.916238im  …    3.54291e-5+8.48023e-6im  -7.62406e-6+1.64215e-5im   0.000112775+2.32236e-5im
 -42.8952-7.91818im   201.806+5.01733im     -21.5038-14.4517im   -12.7938+4.05026im  -23.8118+2.69453im      -0.000160175-6.78447e-5im   4.13392e-5-8.30701e-5im  -0.000855814+0.000338669im
 -4.46939-1.53189im  -46.1351+4.31942im      230.519-0.894398im  -56.4492-3.75571im    51.892+11.1621im       0.000205335-9.20794e-5im  0.000153782+8.43695e-5im  -0.000503379-0.000467105im
 -18.8168-4.67129im   -12.523+5.21867im     -18.6121-8.85055im     93.564+2.37447im  -29.5831+0.571702im       2.46067e-5+1.66097e-5im    2.6309e-5-5.43579e-5im  -0.000178782+0.000304895im
 -7.54098+3.15272im  0.167248-0.0399667im  -0.422169-0.138953im  -1.37431+1.25712im   15.4743-1.82322im       -7.47703e-6+9.26569e-6im   1.48772e-5-9.62472e-6im    2.90185e-5-5.24079e-5im

[:, :, 2] =
 -335.014-5.54885im    6.25747+10.47im     -1.12689-2.98101im    7.85653+5.08438im   54.3551-1.97964im  …   -7.78269e-5-2.28751e-5im     1.70091e-5-3.55928e-5im   -0.000263477-7.29107e-5im
  81.7936+14.1785im   -322.121-7.29737im    38.0402+23.4197im    21.1328-7.9609im    47.4065-2.59343im      0.000407157+0.000171595im   -7.27576e-5+0.000160077im    0.00204946-0.000606909im
  5.50337+7.4361im     85.6178-9.33099im   -376.445+2.6297im     105.401+6.75082im  -111.104-21.1647im     -0.000484024+0.000195714im  -0.000336198-0.000181606im    0.00111332+0.00107764im
  33.4921+6.20733im    21.5971-9.64479im     34.336+16.3788im   -133.735-4.83519im   57.0795-1.4555im       -5.82902e-5-4.01614e-5im    -5.28118e-5+0.000112765im   0.000371606-0.000638284im
  13.6716-5.7004im   0.0437112-0.165886im  0.239548+0.268183im   2.56193-2.29049im  -15.5816+3.05984im       1.83619e-5-1.97876e-5im    -3.06454e-5+1.74626e-5im    -5.98895e-5+0.000109933im

[:, :, 3] =
  109.178+2.16627im    -3.02891-3.68523im    1.38635+0.587794im   -2.99669-1.60396im  -22.5612+1.01748im   …    3.96623e-5+1.59496e-5im   -9.26017e-6+1.79888e-5im   0.000148154+5.81522e-5im
 -32.4015-4.80239im     94.8419+1.50666im   -11.3306-6.56923im    -5.28268+3.6158im   -19.9337-0.700065im     -0.000267965-0.000109666im   1.85941e-5-6.87889e-5im   -0.00122774+0.000200828im
  2.33927-5.95766im    -29.7785+5.44549im    116.625-1.90154im     -39.677-2.7592im    53.7026+8.60292im       0.000287575-9.57117e-5im   0.000168958+8.78634e-5im  -0.000564711-0.000592969im
 -12.1831-0.999153im   -6.49443+3.88604im   -12.8824-6.23348im     36.4343+2.15856im  -23.0858+1.07689im        3.62753e-5+2.40499e-5im    2.31514e-5-5.38497e-5im  -0.000163579+0.000315496im
 -5.09631+2.12097im   -0.315438+0.205453im  0.298656-0.105687im  -0.953936+0.83574im   2.34589-1.03835im       -1.12143e-5+9.74732e-6im    1.38752e-5-6.80201e-6im    2.65997e-5-5.27751e-5im

...

[:, :, 18] =
  -0.676249-0.0177268im  -0.0390984-0.0115825im   0.00113469-0.014869im     0.0116395+0.0242489im   -0.0996344+0.0240919im    …  -1.56936e-7-2.31869e-7im  2.76606e-8+7.84498e-8im  -9.70673e-7-1.09546e-6im
  -0.115476-0.0265847im   -0.777687-0.0926653im    -0.192549-0.0314538im    0.0882896-0.00558733im   -0.061598-0.0824794im        1.91187e-6+1.58777e-6im  1.32606e-6-6.10967e-7im   9.40108e-6+3.17627e-6im
 0.00497367-0.173377im    -0.497138-0.0706967im    -0.432329-0.139627im    -0.0354698-0.0390218im     0.508442+0.113901im        -2.33989e-6+3.14138e-7im  1.87559e-6+1.87429e-6im   2.75776e-6+6.09233e-6im
  0.0790026+0.0665533im   0.0188128+0.0257761im    0.0160278-0.0266972im    -0.534553-0.00398996im   0.0140033+0.00775984im      -2.18616e-7-2.21576e-7im  4.32165e-7+1.98055e-8im   6.35932e-7-1.11434e-6im
 0.00291016+0.0015479im  0.00208272+0.00101381im   0.0166983-0.00112902im   0.0116093-0.00226967im  -0.0939914-0.000645239im      8.30002e-8+1.38662e-8im  1.24766e-7-5.79864e-8im  -6.51902e-8+2.0183e-7im

[:, :, 19] =
  -0.107326-0.00354825im  -0.0529434-0.0275677im    -0.0147906-0.00291988im   -0.001733+0.0218804im    …  -3.79504e-7-2.83202e-7im  9.47979e-8-4.54791e-8im  -1.50251e-6-1.10328e-6im
   -0.26379-0.0617491im   -0.0754421-0.0267948im     -0.314637-0.0713229im   -0.0181322+0.00547635im       2.03747e-6+1.84545e-6im  1.09925e-6-4.02788e-7im   1.23223e-5+8.02536e-7im
  -0.130149-0.142068im     -0.716405-0.106423im       0.334922-0.0960096im    -0.278958-0.0593884im       -2.98152e-6+7.29125e-7im  9.98104e-7+1.09905e-6im   5.24631e-6+9.15017e-6im
  0.0306405+0.0388469im   -0.0643457+0.0274887im    -0.0445321-0.0583908im    -0.458963-0.00843077im      -2.66355e-7-3.0624e-7im    3.2179e-7+3.78091e-7im   2.57665e-6-2.36557e-6im
 -0.0231621+0.0102394im   0.00476338-0.000219589im   0.0125938-0.00110315im  0.00742414+0.000189929im      1.26903e-7-4.26637e-8im  2.72867e-8-4.24455e-8im  -3.78451e-7+4.57819e-7im

[:, :, 20] =
   0.353687+0.00877631im  -0.0584751-0.0371183im    -0.025772+0.0072875im    -0.0124737+0.0171149im    -0.22694+0.0164599im  …  -4.47823e-7-2.51148e-7im  1.29967e-7-1.2875e-7im   -1.49052e-6-8.19747e-7im
  -0.351059-0.0805261im     0.494024+0.0322071im    -0.366634-0.0931863im      -0.10239+0.0137008im   -0.209124-0.0586322im      1.56954e-6+1.54308e-6im  6.76826e-7-1.57354e-7im   1.11851e-5-1.2214e-6im
  -0.222554-0.0957266im    -0.796024-0.116031im      0.893118-0.0416644im     -0.442046-0.0672578im    0.525374+0.0962727im      -2.6615e-6+8.40181e-7im  1.06708e-7+2.08138e-7im   5.65056e-6+9.04069e-6im
  -0.011304+0.0123095im    -0.123663+0.0244327im   -0.0858703-0.0754044im     -0.355745-0.0121001im   -0.151357-0.0224742im     -2.32593e-7-2.911e-7im    1.60499e-7+5.62768e-7im   3.39388e-6-2.60752e-6im
 -0.0417339+0.0159985im   0.00568693-0.00112026im  0.00736822-0.000760594im  0.00302293+0.00200042im  -0.090206-0.0106053im      1.28267e-7-7.36953e-8im  -5.3469e-8-2.35215e-8im  -5.24846e-7+5.16608e-7im
```
Here are more output results:
```
Get_wf computation time: 286.394543 seconds (321.62 M allocations: 215.771 GiB, 1.37% gc time)

CKMS Computation time: 34.425177 seconds (58.60 M allocations: 125.967 GiB, 1.75% gc time)
Number of CKMS iterations: 2950
errK errR :  6.883824319059169e-11 2.153244140484029e-14
```
#### Experiment Dec 2, 2020 3 (job 171)
This experiment is exactly the same as *Experiment Dec 1, 2020 2* except that is runs on the data from "data/KSE_Data/KSE_sol_lin1.jld" (`seed = 2021`). the key parameters are
```julia
nfft = 2^10
par = 500
```

The experiment was run of thelio. Here is the resulting WF:
```
julia> h_wf = load("KSE_wf12_02_20_3-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
  203.584+3.0948im   -2.73015-5.95829im    0.0617399+2.07921im   -4.25393-3.21317im  -27.7599+0.916234im  …    3.54295e-5+8.48031e-6im  -7.62393e-6+1.64217e-5im   0.000112777+2.322e-5im
 -42.8943-7.91897im   201.806+5.01742im     -21.5042-14.4515im    -12.794+4.04995im  -23.8119+2.69468im      -0.000160175-6.78434e-5im   4.13401e-5-8.30705e-5im  -0.000855817+0.000338665im
 -4.47069-1.53199im  -46.1352+4.31918im      230.517-0.895246im  -56.4484-3.75584im   51.8916+11.1618im       0.000205336-9.20794e-5im   0.00015378+8.4368e-5im    -0.00050337-0.000467099im
 -18.8173-4.67237im  -12.5232+5.21893im     -18.6118-8.85043im    93.5643+2.37441im  -29.5834+0.57144im        2.46076e-5+1.66109e-5im   2.63081e-5-5.43567e-5im  -0.000178784+0.000304898im
 -7.54052+3.15237im  0.167209-0.0400003im  -0.422192-0.138921im  -1.37419+1.25744im   15.4744-1.82328im       -7.47832e-6+9.26464e-6im   1.48767e-5-9.6255e-6im     2.90156e-5-5.24112e-5im

[:, :, 2] =
 -335.014-5.54812im    6.25742+10.4698im   -1.12686-2.98079im    7.85693+5.08401im   54.3551-1.97966im  …   -7.78279e-5-2.28754e-5im     1.70087e-5-3.55935e-5im   -0.000263484-7.29038e-5im
  81.7922+14.18im     -322.121-7.29761im    38.0407+23.4194im    21.1331-7.96023im   47.4066-2.59349im      0.000407158+0.000171592im   -7.27589e-5+0.000160077im    0.00204947-0.000606901im
  5.50565+7.43633im    85.6177-9.33061im   -376.441+2.63097im      105.4+6.75095im  -111.103-21.164im      -0.000484024+0.000195714im  -0.000336194-0.000181603im    0.00111331+0.00107763im
  33.4929+6.20911im    21.5975-9.6453im     34.3356+16.3786im   -133.735-4.83493im   57.0799-1.45503im      -5.82918e-5-4.0165e-5im     -5.28103e-5+0.000112762im   0.000371612-0.000638294im
  13.6708-5.69981im  0.0438575-0.165821im  0.239534+0.268036im   2.56169-2.29105im  -15.5818+3.05999im       1.83649e-5-1.97851e-5im    -3.06442e-5+1.74643e-5im    -5.98829e-5+0.000109941im

[:, :, 3] =
  109.178+2.16597im    -3.02884-3.68514im   1.38627+0.587721im   -2.99683-1.60381im   -22.5612+1.0175im    …    3.96631e-5+1.59496e-5im   -9.25991e-6+1.79891e-5im   0.000148157+5.81483e-5im
 -32.4011-4.80283im      94.842+1.50681im  -11.3308-6.56907im    -5.28281+3.61554im   -19.9337-0.700147im     -0.000267964-0.000109665im   1.85949e-5-6.87892e-5im   -0.00122774+0.000200825im
  2.33852-5.95755im    -29.7785+5.4455im    116.624-1.90199im    -39.6766-2.75925im    53.7021+8.60257im       0.000287575-9.57112e-5im   0.000168956+8.78618e-5im    -0.0005647-0.000592959im
 -12.1833-0.999723im   -6.49452+3.88633im  -12.8822-6.23349im     36.4346+2.15847im    -23.086+1.07671im        3.62764e-5+2.40524e-5im    2.31506e-5-5.38484e-5im  -0.000163581+0.000315503im
 -5.09607+2.12077im   -0.315539+0.2054im   0.298695-0.105582im  -0.953839+0.835944im   2.34597-1.03844im       -1.12161e-5+9.7458e-6im     1.38746e-5-6.80288e-6im    2.65954e-5-5.27792e-5im

...

[:, :, 18] =
  -0.676257-0.0177394im   -0.039096-0.0115808im  0.00113284-0.0148618im    0.0116437+0.024251im    -0.0996375+0.0240894im    …  -1.56918e-7-2.31888e-7im  2.76601e-8+7.84491e-8im   -9.7065e-7-1.09556e-6im
  -0.115473-0.0265823im   -0.777681-0.0926655im   -0.192563-0.0314404im    0.0882955-0.00559142im  -0.0615918-0.0824869im        1.91192e-6+1.58769e-6im  1.32607e-6-6.10974e-7im   9.40136e-6+3.17597e-6im
 0.00498817-0.173352im    -0.497152-0.0706785im   -0.432296-0.139628im    -0.0354782-0.0390202im     0.508442+0.1139im          -2.33994e-6+3.14107e-7im   1.8756e-6+1.87428e-6im   2.75766e-6+6.09196e-6im
  0.0790076+0.0665563im   0.0188119+0.0257796im   0.0160366-0.0267008im    -0.534568-0.00398197im   0.0139941+0.00776254im      -2.18662e-7-2.216e-7im    4.32164e-7+1.98308e-8im   6.35806e-7-1.11448e-6im
 0.00290916+0.00154698im  0.0020799+0.001012im    0.0166985-0.00112835im   0.0116094-0.00227187im  -0.0939913-0.000645764im      8.30121e-8+1.38791e-8im  1.24766e-7-5.79882e-8im  -6.51568e-8+2.01864e-7im

[:, :, 19] =
  -0.107334-0.00355509im  -0.0529435-0.027563im     -0.0147926-0.00291193im  -0.00172675+0.0218792im    …  -3.79481e-7-2.83214e-7im  9.47952e-8-4.54784e-8im   -1.5025e-6-1.10334e-6im
  -0.263788-0.0617417im   -0.0754408-0.0267849im      -0.31465-0.0713067im    -0.0181204+0.00546316im       2.03746e-6+1.8454e-6im   1.09922e-6-4.02778e-7im   1.23225e-5+8.02369e-7im
  -0.130141-0.142037im     -0.716424-0.106398im       0.334948-0.0960085im     -0.278967-0.0593893im       -2.98155e-6+7.29093e-7im  9.98121e-7+1.09904e-6im   5.24618e-6+9.14977e-6im
  0.0306453+0.038848im    -0.0643465+0.0274924im     -0.044519-0.0583951im     -0.458974-0.0084236im       -2.66422e-7-3.0629e-7im   3.21799e-7+3.78104e-7im   2.57649e-6-2.36572e-6im
 -0.0231636+0.0102376im   0.00475998-0.000219874im   0.0125941-0.00110165im   0.00742379+0.000188889im      1.26911e-7-4.26439e-8im  2.72924e-8-4.24447e-8im  -3.78427e-7+4.57873e-7im

[:, :, 20] =
   0.353679+0.00877568im  -0.0584782-0.0371111im   -0.0257736+0.0072946im    -0.0124649+0.0171081im   …   -4.4781e-7-2.5114e-7im    1.29961e-7-1.28745e-7im  -1.49055e-6-8.19709e-7im
  -0.351059-0.0805152im     0.494022+0.0322265im    -0.366643-0.093171im      -0.102371+0.0136774im       1.56947e-6+1.5431e-6im    6.76756e-7-1.57337e-7im   1.11852e-5-1.22135e-6im
  -0.222552-0.0956948im    -0.796046-0.116002im      0.893136-0.0416615im     -0.442057-0.0672614im      -2.66151e-6+8.40151e-7im   1.06746e-7+2.0814e-7im    5.65043e-6+9.04038e-6im
 -0.0113006+0.0123076im    -0.123663+0.024437im     -0.085855-0.0754089im     -0.355753-0.012094im       -2.32654e-7-2.91162e-7im   1.60514e-7+5.62765e-7im   3.39373e-6-2.60763e-6im
 -0.0417362+0.0159961im   0.00568305-0.00111876im  0.00736837-0.000758052im   0.0030219+0.00200095im      1.28269e-7-7.3675e-8im   -5.34588e-8-2.35208e-8im  -5.24836e-7+5.16665e-7im

```
Here are more output results:
```
Get_wf computation time: 281.961712 seconds (320.96 M allocations: 213.939 GiB, 1.44% gc time)

CKMS Computation time: 35.197431 seconds (58.60 M allocations: 125.967 GiB, 1.92% gc time)
Number of CKMS iterations: 2950
errK errR :  6.883824319059169e-11 2.153244140484029e-14
```
#### Experiment Dec 2, 2020 4 (job 172)
This experiment is exactly the same as *Experiment Dec 1, 2020 2* except that is runs on the data from "data/KSE_Data/KSE_sol_lin1.jld" (`seed = 2021`). the key parameters are
```julia
nfft = 2^10
par = 100
```

The experiment was run of thelio. Here is the resulting WF:
```
julia> h_wf = load("KSE_wf12_02_20_4-Mo20.jld","dat_h_wf")
5×25×20 Array{Complex{Float64},3}:
[:, :, 1] =
    23.546+0.498655im  -0.0956897-0.482246im   -0.265339+0.34457im    -0.101936-0.401228im    -0.633091-0.0375427im  …   2.11699e-6-1.57125e-6im  -8.27577e-7+9.38101e-7im   7.09028e-6-2.52347e-6im
  -1.14436-0.490648im     26.3936+1.02799im     -1.06679-1.45091im     -0.44176+0.00454645im  -0.448255+0.159925im      -3.58909e-6-1.15393e-6im    4.7067e-6-9.22678e-6im  -3.83957e-5+6.49905e-5im
  -0.43627+1.28838im     -2.40048+0.370856im     33.9142+0.32603im     -2.21926-0.510821im    -0.673633+0.527134im       1.25062e-5-5.97576e-6im   7.05092e-6+8.03293e-6im  -4.55112e-5-3.81674e-5im
 -0.569979-0.720948im   -0.800453+0.401786im   -0.701277-0.537697im     14.5416+0.216571im     -0.38309-0.0769022im      1.55932e-6+2.61717e-6im   2.13921e-6-6.64353e-6im  -2.80985e-5+4.09094e-5im
 -0.714571+0.169574im  -0.0203683-0.0240834im   -0.17481+0.030526im  -0.0275641+0.00597321im    6.08081-0.128193im       1.22972e-7+9.73454e-7im   1.02773e-6-1.41853e-6im   6.19438e-6-7.53156e-6im

[:, :, 2] =
 -28.2528-0.579048im  0.441721+0.650405im    0.468171-0.398349im  0.199241+0.479307im    1.38028+0.0223927im  …  -3.91893e-6+2.20355e-6im   8.16196e-7-1.50839e-6im  -1.58026e-5+3.05743e-6im
  2.67035+0.67333im   -31.3752-1.29624im      2.51876+1.85753im    1.07004-0.0259344im   1.37569-0.219438im       8.84472e-6+2.37315e-6im  -8.83629e-6+1.24034e-5im   7.60318e-5-9.14546e-5im
  1.59642-1.55654im    5.24447-0.426688im    -40.5984-0.459205im   4.47452+0.567925im   0.830111-0.708277im      -1.75844e-5+8.50227e-6im   -1.5203e-5-1.17454e-5im   9.75466e-5+6.0353e-5im
  1.15756+0.83419im    1.60132-0.513073im     1.77967+0.716967im  -15.2334-0.281743im    1.31443+0.140887im      -9.18402e-7-4.12236e-6im  -4.23264e-6+9.34779e-6im   5.17511e-5-5.66377e-5im
  1.22146-0.23941im   0.160803+0.00381285im  0.247715-0.027118im  0.143065-0.0472409im   -5.5307+0.153203im      -3.34851e-7-1.47684e-6im  -1.71295e-6+1.73991e-6im  -1.20308e-5+1.0815e-5im

[:, :, 3] =
   7.00068+0.13601im    -0.427974-0.205472im     -0.189307+0.0755556im  -0.122963-0.114726im   -0.956957+0.00885785im  …   2.05228e-6-5.261e-7im    1.23195e-7+6.10905e-7im   1.11214e-5+1.11467e-7im
   -1.8438-0.195886im     6.60935+0.278521im      -1.68526-0.490495im   -0.727017+0.0361389im   -1.08682+0.0696365im       -7.4497e-6-1.89473e-6im  4.40475e-6-3.45672e-6im  -4.53043e-5+2.38541e-5im
  -1.13275+0.342035im    -3.37339+0.114328im       7.92556+0.116431im    -2.69759-0.114599im     0.16717+0.20543im         4.30121e-6-2.91591e-6im  1.00513e-5+3.7487e-6im    -6.5233e-5-2.15355e-5im
 -0.749299-0.180559im   -0.910001+0.148108im       -1.3238-0.215457im     2.78656+0.0719681im   -1.21377-0.0121646im      -9.43934e-7+7.9557e-7im   2.45201e-6-2.84739e-6im  -2.50535e-5+1.47419e-5im
 -0.626192+0.0657081im  -0.172178+0.00023491im  -0.0513083+0.0100668im  -0.142084+0.0248714im   0.996028-0.032403im         3.5531e-7+4.6108e-7im   6.04527e-7-4.61056e-7im   7.22207e-6-2.72888e-6im

...

[:, :, 18] =
  -0.0186592+0.00155756im   -0.000309672-0.000769564im  0.000634836+0.000715835im    4.0436e-5-0.0012024im   …  -1.32685e-10+2.48362e-9im   4.78653e-9+7.46048e-10im   3.26298e-9+1.25019e-9im
 0.000762495-0.000932682im    -0.0342624+0.000651112im   0.00010868-0.00255335im   0.000360808+3.75209e-5im       9.02631e-9+1.71519e-9im  -8.87386e-9+7.8338e-9im      1.5974e-9+4.93954e-8im
  0.00290174+0.00259233im    0.000830969+0.000553168im   -0.0570539-0.0018783im     0.00137425-0.00208126im      -4.07658e-8-1.1632e-8im   -4.57514e-9-1.71221e-8im    1.70822e-8-4.5864e-8im
  0.00110436-0.00230532im    0.000715949+0.000864973im  0.000295619-0.00061463im    -0.0130203-9.48622e-5im         -1.59e-8+2.11141e-9im   -3.3263e-9+4.03287e-9im    2.01649e-8+1.26062e-8im
 0.000617197+0.000334102im    4.18126e-5-0.000343164im  0.000208557+0.000159429im   7.49906e-5-2.4495e-5im        1.41855e-9+1.08974e-9im  -4.30631e-9-5.384e-10im    -1.03334e-8-7.58617e-9im

[:, :, 19] =
  -0.0252061+0.00116553im   -0.000534437-0.000599186im    0.00063461+0.000534265im   3.57286e-5-0.00106156im   …   1.89327e-9+7.57791e-10im    4.0564e-9+2.94479e-9im    7.63527e-9-2.22514e-9im
 0.000398458-0.000707658im     -0.043868-0.000224428im  -0.000153162-0.00196723im   0.000153889+2.55045e-5im       6.18271e-9+8.538e-10im    -5.62427e-9+7.20665e-9im   -4.07979e-8+1.33577e-7im
  0.00218966+0.00206263im    0.000369033+0.000435617im    -0.0710326-0.00275449im    0.00139762-0.00196347im      -1.93323e-8-1.23112e-8im   -4.2104e-10-1.26522e-8im    1.40611e-8-9.55051e-8im
 0.000986571-0.00197427im    0.000705452+0.000691747im     3.1808e-6-0.00038234im    -0.0127061-0.000169099im     -1.24015e-8+9.41134e-9im   -2.70617e-9-1.73691e-9im   -1.09804e-8+6.03385e-8im
 0.000789891+0.000267544im   -7.57993e-5-0.000279372im   0.000257545+0.000129672im   8.90378e-5-2.06501e-5im       8.3626e-10+2.28326e-9im   -3.29751e-9-4.56214e-10im  -2.54271e-9-1.70961e-8im

[:, :, 20] =
  -0.0301066+0.000807787im  -0.000692312-0.000472727im   0.000595511+0.000398964im   3.71118e-5-0.000931234im  …   3.05647e-9-8.20312e-10im   3.21626e-9+4.55427e-9im    9.36364e-9-5.17686e-9im
 0.000106323-0.000509606im     -0.050235-0.000900326im  -0.000345197-0.00154786im   -2.74729e-5+1.71119e-5im        3.3816e-9-9.90492e-11im  -2.43102e-9+4.55448e-9im   -6.70837e-8+1.92628e-7im
  0.00138953+0.00164384im      8.9548e-6+0.00036386im     -0.0797292-0.00343084im    0.00137845-0.00183741im      -3.02809e-9-1.2458e-8im     2.47938e-9-7.04826e-9im    1.68706e-8-1.25903e-7im
 0.000882195-0.00167997im    0.000661066+0.000559051im  -0.000230106-0.000235372im   -0.0123345-0.000242847im     -9.22636e-9+1.43144e-8im   -1.92441e-9-6.72142e-9im   -3.11923e-8+9.49403e-8im
 0.000917905+0.000221184im  -0.000170187-0.000206939im   0.000289637+9.86601e-5im    9.89973e-5-1.66437e-5im      3.01549e-10+3.22569e-9im   -2.43395e-9-5.68422e-10im   2.98964e-9-2.3537e-8im
```
Here are more output results:
```
Get_wf computation time: 250.268788 seconds (275.21 M allocations: 98.594 GiB, 1.36% gc time)

CKMS Computation time: 5.856998 seconds (14.29 M allocations: 10.694 GiB, 2.16% gc time)
Number of CKMS iterations: 1197
errK errR :  8.813444296750153e-11 2.8987032740866613e-13
```

1:05 PM - Now, I want to
1. get some of the scripts in place to run the reduced models. And
2. look at enforcing the fact that distinct models should be uncorrelated.

5:46 PM - I Wrote the function `modredrun` in the module `Model_Reduction_Dev` in the file "Tools/Model_Reduction_Dev.jl" This should run a generic reduced model given the following
```julia
function modredrun(;
   sig,              # Here only the first M_out = size(h_wf,3) are needed
   sig_m,            # The mean of the signal process
   pred_m,           # the mean of the predictor process
   h_wf,             # The Wiener filter
   Psi,              # The basis functions of the reduced model
   steps,            # How many steps you want to run the RM
   discard,          # How many steps we discard
   noise = false,    # true includes the noise term
   Nosie_dist        # The distribution of the noise term
   )
```

Done for the day.


# Thursday, December 3, 2020

1:15 PM - Tested the new function `redmodrun` in "Tools/Model_Reduction_Dev.jl" on the DWOL and it looks OK. I started investigating this function more closely. I want it to reproduce a generated time series given the right Wiener filter and noise. The way I can verify this is to set the seed and collect the noise sequence. I did this and fond that when the noise values are the same the two series are almost identical observe:
```julia
using Random
using Distributions

# Get software to generate model
gen = include("../Nonlinear Langevin/DataGenDWOL.jl")

# Get model reduction software being tested
mr  = include("../../Tools/Model_Reduction_Dev.jl")
at  = include("../../Tools/AnalysisToolbox.jl")

# Model run Parameters
sigma    = [.4]
sig_init = [1.5]
# Numerical estimate parameters
scheme   = "FE"
steps    = 10^5 #(10^7) # Number of time steps (not including those discarded)
h        = .01
discard  = steps # Number of time steps discarded
gap      = 1

V_c_prime  = x -> -x.*(x.^2 .- 1)

h_wf = zeros(1,3,1)
h_wf[1,:,1] = [1.01 0 -0.01]

Psi(x)  = [x; x.^2; x.^3]

noise = true
noise_dist = MvNormal(zeros(1),sqrt(h)*sigma)
rand(noise_dist)

Random.seed!(2014)
X_rm_c = real(mr.redmodrun(reshape(sig_init,1,:), h_wf, Psi; steps,discard, noise, noise_dist))

Random.seed!(2014)
steps_tot = steps + discard
e = 1/sqrt(h)/sigma[1]*rand(noise_dist,steps_tot)

X_c = @time gen.DataGen_DWOL(; sigma, V_prime = V_c_prime, sig_init, scheme, steps, h, discard, gap, ObsNoise = true,e)

sum(abs.(X_rm_c[1,2:end] - X_c[1][1,1:end-1]))
```
the output is:
```
julia> sum(abs.(X_rm_c[1,2:end] - X_c[1][1,1:end-1]))
3.7973773171079417e-10
```
for `steps = 10^5` (`discard = steps`) and
```
julia> sum(abs.(X_rm_c[1,2:end] - X_c[1][1,1:end-1]))
1.6955452561287845e-8
```
for `steps = 10^7`. So, the reduced (reproduced, rather) model did exactly what the true model did.

11:15 PM - Today I also met with Dr. Lin and we looked at the data (KSE) and the resulting Wiener filter for parameters values (`par = 1500` and `nfft = 2^12`). A few things were noted:
1. I am scaling the Fourier transform differently from the way Dr. Lin did it and he suggested that his way is better for our purposes. I will need to investigate this issue. This was surmised from the fact that the variances of the model were on the order of 100 000.
2. In consideration of the above it seemed that the computed Wiener filter may be worth running in a noise-free reduced model run. Before that though I will investigate at which point the WF coefficients decay enough that truncation is not a big issue. May ways of doing this were discussed.
3. After a suitable cut off was purposed I would run the noise-free reduced model.


# Friday, December 4, 2020

11:50 AM - The first thing I want to work on today is to check the scaling (or lack thereof) associated with the DFT I am using.  Here I have a slow but intuitive implementation of the same transforms that `FFTW.jl` uses, and demonstrate a little test.
```julia
using FFTW
using Statistics: mean

function my_dft(u)
	N = length(u)
	v = [sum(u[j]*exp(-im*2pi/N*(j-1)*(k-1)) for j = 1:N) for k = 1:N]
end

function my_idft(v)
	N = length(v)
	u = [mean(v[k]*exp(im*2pi/N*(k-1)*(j-1)) for k = 1:N) for j = 1:N]
end

nfft = 2^10
Theta = 2π*(0:nfft-1)/nfft
Z = exp.(-im*Theta)

f(z) = 2z^(-1) + 6 + 2z
R = f.(Z)
S = fft(R)/nfft
S_my = my_idft(R)
```
The result was
```
julia> sum(abs.(S-S_my))
1.9216778795533907e-11
```
So, Julia (FFTW) puts the mean with the `ifft`. This means the `vv` in the function `my_KSE_solver` is unscaled and then scaled when it is converted back to `uu`. This make `vv` very large so I would just as soon not have that. So, today, I switched each instance of `fft` with `ifft` and vice versa. Now to run it on thelio. This is just a repeat of what was done at the beginning of Nov 24, 2020, for reference I repeat my self a little: I ran "Examples/KSE/KSE_data_gen.jl" (job 174 at Fri Dec  4 14:41:00 2020) to generate KSE data and save it as "data/KSE_Data/KSE_sol_linn.jld". Here, `gen = "linn"`. So, now we have a copy of the data of thelio, with the standard `gen = "linn"` parameters:

```julia
gen = "linn"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2 # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100
seed     = 2020
```
The goal will be to view this data on the "KSE_data_analyzer" notebook.

It finished at about 3:29 PM. And seems not to have worked. So, I found a typo in the code, I had missed an instance of fft to be changed to ifft. I fixed it and checked it in a new notebook "KSE_data_gen_test". The data looked good. So, I just sent it to thelio as job 175 at Fri Dec  4 16:16:00 2020.



# Saturday, December 5, 2020

2:14 PM - Thelio is down so I will run the data on my computer. I ran "KSE_data_gen.jl" in Atom, in it's own window. I started it at 2:16 PM.
8:54 PM - I don't know when it finished, but I now have a set of data on my computer to work with. I ran it through the "KSE_data_analyzer" notebook and saw the the data was much smaller in absolute value but I also found that all the predictors were of the around the same order, that is they were order unity or smaller. This is because before the inertial manifold predictors were much bigger in magnitude than the signal itself. These terms are quadratic in the terms of the signal process.

```julia
gen = "linn"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

T        = 10^5 # Length (in seconds) of time of run
T_disc   = 10^5 ÷ 2 # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = 1e-3 # Timestep
g        = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap  = 100
seed     = 2020

Random.seed!(seed)
uu, vv, tt =  kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)
```
The data was saved in "KSE/Data/KSE_sollinn.jld" (on my desktop).

#### Experiment Dec 5, 2020 1
Here I compute the Wiener filter for the data generated on my computer. This data was then loaded and a Wiener filter was computed. Here is the code I used along with the solution.
```julia
M_out = 100
Δt = h*obs_gap

par = 1500
nfft = 2^12

Psi(x) = kmr.PSI(x; h, obs_gap, P, N)

@time h_wf = mr.get_wf(X, Psi; M_out, par, nfft, rl = false, PI = false)
```
The result was:
```
5×25×100 Array{Complex{Float64},3}:
[:, :, 1] =
  9349.87-1255.54im  -233.103+169.473im  …  -11.5853-5.06667im
  133.295-930.953im   7147.62-594.171im     -7.12174-3.59719im
 -747.302+394.273im  -1655.48-932.133im     -18.7627-27.4841im
  2516.65+2272.37im   540.984-320.928im     -18.2888+12.4662im
 -86.6757+407.361im  -36.4745+52.2861im      7.88715+5.72918im

[:, :, 2] =
  1879.78+1938.71im   178.546-321.55im   …    34.578+13.6545im
  111.045-1233.06im   2284.96+996.729im      21.5259+8.60318im
 -668.679+758.057im  -1165.18-620.042im      51.1214+71.2045im
  295.955+1087.76im    440.94+218.167im      54.6638-38.8932im
  -161.55+62.4245im  -9.80719+35.5088im     -22.3996-15.9665im

[:, :, 3] =
 -9517.05+166.643im   334.363-116.327im  …  -29.6608-10.2677im
  1353.09+2593.49im  -6123.03-584.627im     -20.5444-2.09551im
 -1067.03-980.211im   2071.81+1069.72im     -40.1492-43.7478im
 -1588.29-2321.62im  -327.739-104.45im      -44.0172+36.1104im
 -47.1073-401.068im   97.1253-101.104im      17.7296+12.1839im

...

[:, :, 98] =
 -79.0629-588.198im  -66.3832+756.56im   …  -1.36171+0.335836im
 -977.496+1788.77im  -63.1384-539.863im      2.50092-4.33989im
  1676.07-3361.66im   394.014-14.5664im     -2.00886+3.18274im
 -323.927-822.25im   -268.681+60.2782im     0.414216-0.477683im
 -184.819-129.395im  -39.8059-17.3266im      0.22577-0.0295624im

[:, :, 99] =
 781.357-542.832im    -326.9+432.729im  …    1.51566+0.05337im
 1508.15+691.723im   297.618+455.677im      -4.48108+3.61844im
 2046.35-645.228im   305.331-472.463im      0.315305+2.43262im
 1279.73+43.8159im  -100.379+262.281im      0.152594+0.712349im
 176.589+42.4591im   61.8087+12.7585im     -0.121087-0.0740854im

[:, :, 100] =
  127.238-75.906im   -380.392-160.87im   …   -1.90492-0.588557im
  -1162.0-2338.46im   -1229.4+1755.66im       5.47405-4.27488im
 -2681.35-1309.98im   1320.45-847.341im       1.19182-5.28527im
  1.89863-127.924im  -525.121+239.076im     -0.828088-0.968056im
 -40.1118+126.569im  -19.3036-44.9249im      0.146902+0.143395im
```
Also,
```
 Get_wf computation time: 710.222636 seconds (411.03 M allocations: 861.093 GiB, 9.44% gc time)

 CKMS Computation time: 507.801449 seconds (169.25 M allocations: 773.837 GiB, 11.57% gc time)
 Number of CKMS iterations: 6103
 errK errR : 8.148016231247954e-11 4.512741213213582e-14
```
I think this will not be a very good Wiener filter. However observe the following:
```julia
sig  = X[:,2:end];
pred = mr.get_pred(X[:,1:end-1], Psi);

sig_m  = mean(sig, dims = 2);
pred_m = mean(pred, dims = 2);
C = sig_m
```
```
julia> sum(h_wf[:,:,k]*pred[:,M_out - k + 1] for k = 1:M_out) + C
5×1 Array{Complex{Float64},2}:
   0.6817769290154864 + 0.2837080449283443im
  -0.9962927993676172 - 0.09952869236130044im
  -1.8758810133270523 + 0.6414698155927139im
   0.6861379865270217 - 0.1852998998267594im
 -0.23151584478867768 - 0.09224841253651911im

julia> sig[:,M_out:M_out+1]
5×3 Array{Complex{Float64},2}:
   0.49894+0.166077im   0.495678+0.162481im   0.491848+0.15849im
  -1.76475-0.555181im   -1.73796-0.562808im   -1.71021-0.569899im
 -0.753714+0.810655im   -0.75247+0.834937im  -0.749829+0.857893im
  0.636802-0.387014im   0.651121-0.343075im   0.664037-0.298592im
 -0.203709-0.126904im  -0.201219-0.119974im  -0.197976-0.112141im

julia> norm(sig_rm - sig[:,M_out])/norm(sig[:,M_out])
0.6257158457486914

julia> norm(sig_rm - sig[:,M_out])
1.463322098376953
```
Despite the outrageous looking coefficients the relative one-step prediction error is much smaller than I though it would be.


# Monday, December 7, 2020

1:35 PM - Today, I am investigating the Wiener filter after I have modified the KSE data generating code. This is all described above. Today, I am back on thelio. I will repeat the experiment above (*Experiment Dec 5, 2020 1*)

#### Experiment Dec 7, 2020 1 (Rerun on thelio)

This experiment was run in the notebook "Examples/KSE/KSE Model reduction.ipynb"

Here is the code
```julia
using JLD
using Dates
using PyPlot
using Statistics: mean
using LinearAlgebra: norm

mr  = include("../../Tools/Model_Reduction_Dev.jl")
at  = include("../../Tools/AnalysisToolbox.jl")
kmr = include("KSE_modredTools.jl")

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")
# sol_file = server ? "../../../data/KSE_Data/KSE_sol$Exp.jld" :
#    "Examples/KSE/Data/KSE_sol$Exp.jld"
# println("Sol save location: " * sol_file)
# wf_file = server ? "../../../data/KSE_Data/KSE_wf$Exp-Mo$M_out.jld" :
#    "Examples/KSE/Data/KSE_wf$Exp-Mo$M_out.jld"

# When I want the standard lin et al. (2017) data.
sol_file = server ? "../../../data/KSE_Data/KSE_sol_linn.jld" :
   "Data/KSE_sol_linn.jld"

 M_out = 100
 Δt = h*obs_gap

 par = 1500
 nfft = 2^12

 Psi(x) = kmr.PSI(x; h, obs_gap, P, N)

 X = vv[2:d+1, 1:end]

 @time h_wf = mr.get_wf(X, Psi; M_out, par, nfft, rl = false, PI = false)
 ```
Here is the result:
```
julia> h_wf
5×25×100 Array{Complex{Float64},3}:
[:, :, 1] =
  10545.7-794.949im   232.383+1179.6im   …  -20.1252-10.494im
  593.575-1799.06im   12685.9-1111.27im     -11.8238-10.7905im
 -911.765+3735.54im   200.106-1160.23im      39.9278-27.5552im
  57.7441-42.4492im  -398.934-200.039im      3.85142-38.91im
 -206.186+274.822im   14.6174-47.9232im     -9.89727+27.9757im

[:, :, 2] =
 -923.329+1049.9im   -257.06+35.9999im  …   43.9612+30.9441im
 -1283.12-2704.16im   881.87+2394.12im      24.7885+17.2973im
 -2448.36+2692.45im  36.5887-23.7307im     -103.841+26.8297im
  394.477+678.363im  -466.74+11.9714im     -18.2804+86.1575im
  110.848-135.379im  34.0694+22.3472im      27.7382-74.3678im

[:, :, 3] =
 -8715.85+832.049im   137.084-1023.18im  …  -22.6788-21.6486im
  1215.29+743.874im  -11344.4-356.744im     -5.70314-1.44347im
  453.665-2426.23im  -353.878+562.911im      74.4533+14.1772im
  293.217+391.919im   -293.16+443.571im      22.2058-48.2931im
  192.215-271.336im  -12.0341+29.414im      -21.1788+53.165im

...

[:, :, 98] =
   1051.4+710.149im  -129.466-34.2065im  …  -0.225742+0.490731im
 -939.966-859.626im   292.121+729.808im       8.53509+8.78829im
 -1233.31+1505.86im  -183.629+106.972im       1.87664+10.7836im
 -671.687+759.237im   169.179-247.71im       -1.39865+4.03962im
  83.2924+56.5083im   16.5698-35.9778im      0.108188-0.0874232im

[:, :, 99] =
  519.149+1146.91im  -628.011+399.565im  …   0.297174-0.903767im
 -3657.41-516.108im   1612.15-488.185im       3.30528+4.25188im
 -1156.27+3548.35im   596.615-958.239im      -1.36434-3.68165im
  581.055+906.052im  -519.433-531.896im      -1.51742-1.17972im
 -105.143-123.866im   12.2814+1.25605im     0.0790085+0.0356196im

[:, :, 100] =
  1307.79-233.191im  -789.205-301.516im  …   -1.33643+2.63357im
  1046.75-2831.8im    17.9295+579.418im      -6.28936-2.99582im
 -4160.22+249.454im   3269.37-798.956im       0.17423-8.82939im
  1434.13-635.662im  -464.613-525.253im     -0.353796-2.80418im
  -43.064+28.3919im  -14.6587-12.9168im     -0.217532+0.0821389im
```
Here are more output results:
```
Get_wf computation time: 441.084422 seconds (411.78 M allocations: 924.783 GiB, 1.92% gc time)

CKMS Computation time: 218.085101 seconds (169.31 M allocations: 837.491 GiB, 2.18% gc time)
Number of CKMS iterations: 6608
errK errR : 9.577601950902882e-11 7.172692752063925e-14
```


# Tuesday, December 8, 2020

10:16 AM - Today I am at my office hoping to get a lot done!
The first thing I want t do is think about the data that I have just generated after the redefinition of the discrete Fourier transform. I will quickly compare some aspects of it with what Dr. Lin has published. This should not take more than 20 mins.

11:45 AM - This ended up taking about 90 mins. In the time I finally added `.ipynb` to my `.gitignore`. So, I will be using jupytext a lot more now. Anyway, I compared the KSE statistics and found some the same and some different. The energy looks off by a factor of ten. That is to say the 2017 series is 1/10 in energy, from mine. this suggests the omission of a multiple of a time-step perhaps.  While I think of something better to do I think I will have thelio make another times series.

3:40 PM - I did not make the time series, as I said I would. I have instead been looking at the code. I verified that the code (when not correcting for aliasing) behaves the same as Trefethen's MatLab code. I looked at the first and last line of both outputs and they were identical. When I corrected for aliasing there was only a small change in the output visually. Now, I switched the definition of the discrete Fourier transform (putting the mean in the forward transform rather than the backward as both FFTW and MatLab do). This gave me a very different solution. So, there is more I need to do when change I switch the definition of the transform. This is what I will focus on in the last hour before I have to go home. The only thing else I need to change is the sign of 𝑖 that comes down when I take the derivative with respect to 𝑥. This is because in julia (FFTW) the negative exponent is with `fft` (as it is with the our work here) so since I switched `ifft` to `fft` I needed to replace the coefficient to the transform of the derivative of 𝑢². Once I did that the code again agreed with Trefethen's. I have now run this code now in thelio to get the data that should more closely agree with Dr. Lin's.


# Wednesday, December 9, 2020

1:16 PM - Just got through talking a while with Eric. Before that I finished some class stuff. What I want to do now is get the data analyzer working and then run the reduced model.

# Thursday, December 10, 2020

8:26 PM - Yesterday, I ran some simulations and analyzed the data. I got the same autocorrelations as before the change. Which are different from the earlier data which more closely resembled Dr. Lin's results. I also cached the notebooks and added `*.ipynb` to the `.gitignore`. Hopefully this works well for me.

Today, I plan on investigating the solver to see why it does not give the other results. I think I will run a few experiments.
1. If I just supply the scalers myself what happens.
   ```julia
   dft(X)       = fft(X)/size(X,1)
   idft(X)      = ifft(X)*size(X,1)
   plan_dft(X)  = plan_fft(X)/size(X,1)
   plan_idft(X) = plan_ifft(X)*size(X,1)
   ```
2. Then try this again:
   ```julia
  dft(X)       = ifft(X)
  idft(X)      = fft(X)
  plan_dft(X)  = plan_ifft(X)
  plan_idft(X) = plan_fft(X)
  ```
  This one requires a change on the coefficient of 𝑖.

I don't think I like this, I need to think of something else, a fresh session sort of thing. So, When I use
```julia
dft  = ifft
idft = fft
```
and make the appropriate change in the coefficient 𝑖. Then the code behaves just the same as the original solver on the Trefethen data. (not correcting for aliasing). I will verify this again.

I already have the module `Model_KSE` in "Examples/KSE/Model_KSE.jl" which I had altered and then reverted back to as it was before. I have verified that this script (when not correcting for aliasing) preforms the same as Trefethen's MatLab code, I compare the final term in the series and they were the same. Then I duplicated this modele and made `Model_KSE_Dev` in "Examples/KSE/Model_KSE_Dev.jl". I then changed the definition of `dft` used using
```julia
dft(X)       = ifft(X)
idft(X)      = fft(X)
plan_dft(X)  = plan_ifft(X)
plan_idft(X) = plan_fft(X)
```
together with a change in the sign of the coefficient of 𝑖.

I ran both solvers on the Trefethen problem, in a new session, both not accounting for aliasing. And they are identical here is the script with result:
```julia
using PyPlot
using Statistics: mean, var

kse = include("../Model_KSE.jl")
ksed = include("../Model_KSE_Dev.jl")

T        = 150 # Length (in seconds) of time of run
T_disc   = 0    # Length (in seconds) of time discarded
P        = 32π  # Period
N        = 128  # Number of fourier modes used
h        = 1/4  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
obs_gap  = floor(Int, T/h/100)

Δt = h*obs_gap
uu, vv, tt    =  @time kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)
uud, vvd, ttd =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)
```
The result:
```
julia> sum(abs.(uu - uud).^2)
0.0
```

Now, When I do correct for aliasing in both algorithms and run the entirely same script, I get
```
julia> sum(abs.(uu - uud).^2)
259442.24694821244
```

I conclude that the dealiasing part needs to be adjusted a long with the change in the dft, though I don't understand why it should.

12:42 PM - I am taking a break from this (since it is becoming unproductive) to grade.

### Meeting with Dr. Lin
The issue I wanted to discuss surrounded the following facts.
1. When I do not my solver (from "Examples/KSE/Model_KSE") with out correcting aliasing my solution matches very well Trefethen's code, it reproduces his solution.
2. When I switch for `dft = ifft` (in every instance) and change the coefficient of the first derivative this code with aliasing also matches.
3. When I correct for aliasing in both implementations. The solutions are very different.

We discussed this at length and concluded that the de-aliasing lines of code where not correct.   


# Friday, December 11, 2020

2:23 PM - I have been working on this dealiasing code. But it has proven very enigmatic to my methods which are perhaps too sloppy. When I tried rerunning, what I thought was exactly the same code as yesterday, the solution was populated entirely of `NaN`s. I went back to the Trefethen problem and was able to reproduce those results. For the past few hours I have been investigating this issue and feel no closer to it's resolution. My methods of investigation are way too causal. This is why I am grateful for this journal. I need it to help me systematically investigate the problem. Without it, I seem to just run a bunch of things and try to remember it all. I am taking a break now to grade and will return with fresh eyes. I think it will be really good and worth while to write a tutorial about aliasing with regards to KSE (geared to an undergraduate audience I guess)

4:54 PM - Very interesting! I have been struggling with a coding problem I had coded, with in my solver:
```julia
if aliasing
    Nh    = ceil(Int,N/2)
    v_pad = [v[1:Nh]; zeros(2N); v[Nh+1:end]]
    F     = plan_ifft(v_pad)        # julia's ifft is my fft for this problem.
    iF    = plan_fft(v_pad)         # julia's fft is my ifft for this problem.

    function NonLin(v)
        v_pad = v_pad = [v[1:Nh]; zeros(2N); v[Nh+1:end]]
        nv    = F*(real(iF*v_pad)).^2
        nv[1:N]
    end
else
    ## Not correcting for aliasing
    F = plan_ifft(v)          # julia's ifft is my fft for this problem.
    iF = plan_fft(v)          # julia's fft is my ifft for this problem.
    NonLin(v) = F*(real(iF*v)).^2
end
```
the problem was that things were behaving very funny if I ran `aliasing = true` I get:
```
julia> uu_a, vv_a, tt    =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap,aliasing = true)
ArgumentError: FFTW plan applied to wrong-size array
```
Even in a fresh session. So for some reason the function definition in the else statement is beign referenced on `NonLin` calls. Similarly if I run: `aliasing = false` I get:
```
julia> uu_a, vv_a, tt    =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap,aliasing = true)
UndefVarError: NonLin not defined
```
I think that julia doesn't know where to get the function definition which is so odd. Anyway, I looked on the discourse and found a solution that seems to work (though I think multiple dispatch might be an issue with this fix, but that's fine with me). The key is using anonymous functions as in this:
```julia
if aliasing
    Nh    = ceil(Int,N/2)
    v_pad = [v[1:Nh]; zeros(2N); v[Nh+1:end]]
    F     = plan_ifft(v_pad)        # julia's ifft is my fft for this problem.
    iF    = plan_fft(v_pad)         # julia's fft is my ifft for this problem.

    NonLin = function (v)
        v_pad = v_pad = [v[1:Nh]; zeros(2N); v[Nh+1:end]]
        nv    = F*(real(iF*v_pad)).^2
        nv[1:N]
    end
else
    ## Not correcting for aliasing
    F = plan_ifft(v)          # julia's ifft is my fft for this problem.
    iF = plan_fft(v)          # julia's fft is my ifft for this problem.
    NonLin = v -> F*(real(iF*v)).^2
end
```
The source is https://discourse.julialang.org/t/defining-a-function-inside-if-else-end-not-as-expected/13815

I will have to verify this more completely later right now I have to go. But it worked in my simple tests.



# Wednesday, December 16, 2020

3:12 PM - I am very close to being done with teaching this semester. So, I can start ramping up my research.

Now, the first thing I want to do today is to test the new solvers that have optional aliasing. I need to make sure that code is working. The I will test the different de-aliasing algorithm. So, here is a goal:
##### Goal:
Get a run of data that agrees with the 2017 paper.



# Thursday, December 17, 2020

The thing to look for is that the energy decays exponentially.

3:58 PM - I have been trying to fix the KSE solver.



# Friday, December 18, 2020

Today I have been relearning all about aliasing and de-aliasing. So, now I will attempt to summarize what I have learned so far. I wish this had a latex editor.



# Monday, December 21, 2020

12:21 PM - Today I will rewrite the KSE solver a little bit. I start with an old version that predates recent iterations the idea will be to document all changes made and they worked.

One thing I noticed is that I was actually only using half as many modes as I thought I was. This is because I double count the modes by there conjugates with opposite wave number. These provide no new information since the original  signal is real.



# Monday, December 28, 2020

12:00 PM - The goal for today is to get my KSE solver working. I'll start by running the the oldest version that worked then move towards where I want it to go making sure there are no NaN's. I will mostly just focus on the 2017 paper parameters. I will being testing them using the KSE script "Aliasing_KSE.jl".

12:10 PM ran current `my_KSE_solver` from "Model_KSE.jl".
```julia
kse  = include("Model_KSE.jl")

T        = 1000 # Length (in seconds) of time of run
T_disc   = 0    # Length (in seconds) of time discarded
P        = 2π/sqrt(0.085)  # Period
N        = 96  # Number of fourier modes used
h        = .001  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
obs_gap  = 100 #floor(Int, T/h/100)

Δt = h*obs_gap

uu_a, vv_a, tt   =  @time kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap);
```
And I get
```
julia> findfirst(isnan,sum(uu_a[:,:],dims = 1))
CartesianIndex(1, 3635)
```
I would like to checkout and old working version. So I went back to just before I made the file a `module` this was commit: `1be8fb7c36cfb96562020ae8cf6a8b09f7db2b7c`. I run exactly the same above and get
```
julia> findfirst(isnan,sum(uu_a[:,:],dims = 1))

julia>
```
Which I believe means it returned `Nothing`. That is, there were no `NaN`'s. I also checked it by looking at `sum(uu_a[:,end])` which returned `≈103.0900`.

So, here there is my starting point.

1:35 PM - I have been comparing the performance of my code with that of Trefethens (in MatLab). Here is a summary:
When I run the original parameters:
```julia
T = 150
P = 32π # Period
N = 128 # Number of fourier modes used
h = 1/4 # Timestep
g = x -> cos(x/16)*(1 + sin.(x/16)) # Initial condition function
T_disc = 0
n_gap = 6 # 1 +  No. of EDTRK4 steps between reported data
```
The solutions are very close. But when I run the following:
```julia
T = 1500
P = 32π # Period
N = 128 # Number of fourier modes used
h = 1/4 # Timestep
g = x -> cos(x/16)*(1 + sin.(x/16)) # Initial condition function
T_disc = 0
n_gap = 60 # 1 +  No. of EDTRK4 steps between reported data
```
My solution is all `NaN`'s after 23 steps and the Trefethen MatLab code does not get any `NaN`'s. So, there is something different about the implementations. First thing I want to do is run this experiment of Thelio. Then I think I will see if Julia (FFTW) and MatLab do `fft` the same way. I ran the code in Thelio and got the same results.

```
julia> kse = include("Model_KSE.jl")
Main.Model_KSE

julia> uu, vv, tt = kse.my_KSE_solver(;T_disc = 0);

julia> uu
128×101 Array{Float64,2}:
 1.0478    0.96252   0.882964  …   1.2048      1.18154     1.304
 1.09273   1.00902   0.92831       1.96644     1.79324     1.55974
 1.13432   1.05366   0.972725      2.017       1.55383     0.809059
 1.17213   1.09608   1.01596       0.759738    0.155579   -0.761331
 1.20573   1.13591   1.05776      -1.16615    -1.52324    -2.03811
 1.23473   1.17277   1.09786   …  -2.28915    -2.26304    -2.25764
 1.25874   1.20627   1.13595      -2.10622    -1.89944    -1.67577
 1.27743   1.23602   1.17174      -1.23856    -1.04149    -0.854294
 1.29049   1.26161   1.20491      -0.321678   -0.155334   -0.0508895
 1.29766   1.28263   1.23508       0.41455     0.645473    0.756064
 ⋮                             ⋱                           ⋮
 0.570326  0.514955  0.469216     -2.04773    -1.78197    -1.4788
 0.624347  0.56359   0.513247  …  -1.34651    -1.10752    -0.908016
 0.679155  0.613166  0.558214     -0.514224   -0.305007   -0.221851
 0.734333  0.663414  0.603929      0.0714131   0.245902    0.252322
 0.789444  0.714051  0.650203      0.3116      0.413137    0.381143
 0.844034  0.764787  0.696842      0.264425    0.260785    0.233713
 0.89764   0.815318  0.743646  …   0.114244    0.0308104   0.0628739
 0.949787  0.865334  0.790409      0.116605    0.0314568   0.161842
 1.0       0.914511  0.836921      0.47624     0.444077    0.646349

julia> uu, vv, tt = kse.my_KSE_solver(1500;n_gap = 60,T_disc = 0);

julia> uu
128×101 Array{Float64,2}:
 1.0478    0.512003  0.323489  …  NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.09273   0.539104  0.33607      NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.13432   0.56586   0.351116     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.17213   0.592363  0.371264     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.20573   0.619406  0.398449     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.23473   0.647726  0.432208  …  NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.25874   0.677948  0.468828     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.27743   0.709314  0.501755     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.29049   0.740014  0.524222     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.29766   0.767291  0.532768     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 ⋮                             ⋱         ⋮                        ⋮
 0.570326  0.275668  0.190243     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.624347  0.300571  0.206365  …  NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.679155  0.325793  0.222428     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.734333  0.35138   0.238324     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.789444  0.377251  0.253939     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.844034  0.403531  0.269229     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.89764   0.430172  0.284036  …  NaN  NaN  NaN  NaN  NaN  NaN  NaN
 0.949787  0.457264  0.298114     NaN  NaN  NaN  NaN  NaN  NaN  NaN
 1.0       0.484588  0.311193     NaN  NaN  NaN  NaN  NaN  NaN  NaN

julia> findfirst(isnan,sum(uu,dims = 1))
CartesianIndex(1, 23)
```
And so, I get the same results on thelio. I investigated the MatLab `fft` function and found that there algorithm is based off of none other than the `FFTW` library, just like Julia and that the definition for the DFT is the same. So, it is neither of the two problems I suggested earlier. Let me look at the energy spectrum. I found a few interesting things which I hope to communicate here.

#### Experiment Dec 28, 2020 1

For this experiment I use the Trefethen parameters except that I run it for `T = 1500`. This is just as described above. The interesting piece is when I plot the energy at various times. I did this using the following code.
```julia
findfirst(isnan, sum(uu_a,dims = 1))

EE = abs2.(vv)

for i in 100:10:216
    semilogy(EE[:,i],label = i,lw  = 1)
    legend()
end
```
The blow-up occurred at index 214 with equates to time (214 × `n_gap` × `h` = 214×6/4 = 321). So, we plot the energy profile (|v_k|²; k in [0:N÷2-1; 0; N÷2-N+1:-1]), at each time and what I find is that the energy begins to accumulate at around `k = 13` This seems to be what contributes to the blow up.

Is there supposed to be a conservation of energy? Where is this energy coming from?


# Tuesday, December 29, 2020

Today I reached way back to get a solver that produces no NaN's. I ran it on the full 2017 parameters and it came out with no NaN's. This is with the naïve de-aliasing. I feel like I have a good understand of the state of the code. One question I still have is how does Trefethen's code do so well with out de-aliasing? His code run in MATLAB is stable and does not seem to have the aliasing artifacts I have noticed in my aliased data.   

2:41 PM - Now I attempt to change the scaling, this is basically a change in definition of the DFT I am using. I changed the scaling and it all works just find testing it on the Trefethen code I get the same solution.

3:57 PM - Now I will change the de-aliasing routine. Before I did that I wanted to address the issue of saying I am using so many Fourier modes but because the process is real there half on the modes are just the conjugates of the other so, there really is only half the information in it.

So, what I did was set a new parameter `n` and set `N=2n`.


# Wednesday, December 30, 2020

12:03 PM - Today I will start with verifying that the code is de-aliasing in an appropriate way. So, here is the setting:  


# Monday, January 4, 2021

10:00 AM -  I have been writing a note on aliasing and de-aliasing. I feel like I understand it pretty well. However, I can't explain why Trefethen's MatLab code is more stable than mine is. I think I will give myself 2 hours (till about noon) and then I will look at Dr. Lin's code to see how I might tweak mine.

3:00 PM - I worked on the tutorial and did some interesting calculations. But I did not achieve what I set out to do. Though I do feel more comfortable with the complexities in the indexing between the mathematical DFT and the numerical `fft`.



# Tuseday, January 5, 2021

12:36 PM - Though I have been going through the code very carefully, which processes has been facilitated by writing the de-aliasing tutorial and the FFTW tutorial I have not been able to identify why the solver is unstable, when it seems to mimic exactly Trefethen's solver in MATLAB.

4:25 PM -  I have been working on a large communication to Dr. Lin in a google doc. Here is how I got the first graphic
```julia
uu, vv, tt = ksed2.my_KSE_solver(150;
                               N = 128,
                               T_disc = 0,
                               n_gap = 6)
uu

plot(32π*(0:127)/128,uu[:,end], label = "me (Julia)")
title("KSE solution at T = 150 (both with aliasing)")
```


# Thursday, January 14, 2021

Today the goal is to write out algorithm by hand. Just, go through the whole routine, almost from scratch.



# Friday, January 15, 2021

7:22 AM - The plan today is to completely rewrite a new KSE solver. I hope to have a complete first draft by noon. Basically, this is nothing mare than a system of ODE's though nonlinear. So, first (1) I will write my ODE solver using ETDRK4 using the Kassam-Trefethen approach of writing the division as a contour integral. The solver will be of the usual general form of `d/dt(x) = f(x)`. The principal inputs will be the function `f` (in-place function, like Dr. Lin does, I think), `h` time step. Then I will write a function for the RHS, the "stepper" as it were, complete with de-aliasing. So, that is really all, I think. two main parts.

12:16 PM I have a first complete draft of the ODE solver capable of using ETDR4. I am testing this. Hopefully tomorrow I can have it working.



# Saturday, January 16, 2021

4:58 PM - I am going to test the code and hopefully get the ODE' solver working.

# Monday, January 18, 2021

1:19 PM - Tested the "myODE_solver.jl" using a diagonal linear example, a stiff example, and Lorenz63. It looks like it is working in all cases. The task now is to us it KSE.

3:39 PM -  I tested "myKSE_solver.jl" and got the same results as the oldder code, "Model_KSE.jl". I think that there is something then wrong with my dealiased nonlinear part. I will look at that tomorrow.

Tomorrow I will look at Dr. Lin's dealiased nonlinear part. This will involve looking up functions like `make_r2r` in FFTW. If I use that does it work? Also, I wonder if my ETDRK4 solve is no good because it dosen't seem to have much of an advantage over RK4. Maybe I will throw something on stack overflow about ETDRK4.



# Saturday, January 17, 2021

2:20 PM - I have been comparing Dr. Lin's implementation and my own of the KSE solver. Here is what I conclude so far:

Dr. Lin's nonlinear implementation is different from mine. Here is the evidence:
```julia
using FFTW
ks   = include("C:/Users/JaredMcBride/Desktop/DDMR/klin/ks.jl")


## Set up
T       = 150
P       = 32π
n       = 64
h       = 1/4
g       = x -> cos(x/16)*(1 + sin(x/16))
T_disc  = 0
n_gap   = 6
aliasing= true

N = 2n+1

# Spatial grid and initial conditions:
x = P*(0:N-1)/N
u = g.(x)
v = fft(u)/N        # The division by N is to effect the DFT I want.

# Now we set up the equations
q = 2π/P*[0:n; -n:-1]
ℓ = -0.5im*q


## Dr. Lins function
fillout(v::Array{Complex{Float64},1}) = [0; v; reverse(conj(v))]

NonLin! = ks.make_ks_field(n; alpha = 0, beta = 0, L = P)

NonLin = function (v)
    u = zeros(Complex{Float64},n)
    NonLin!(v[2:n+1],u)
    fillout(u)
end

## My function
pad = (3N+1)÷2
K = N + pad
NonLinNA = function (v)
    v_pad = [v[1:n+1]; zeros(pad);v[n+2:N]]
    nv = fft(bfft(v_pad).^2)/K
    Nv_dealiased = ℓ .* [nv[1:n+1]; nv[end-n+1:end]]
    # ifftshift(conv(fftshift(v),fftshift(v))[N-(n-1):N+n])/N
    # v_pad = [v[1:n]; zeros(pad);v[n+1:N]]
    # nv = F*(real(iF*v_pad)).^2*K/N
    # [nv[1:n]; nv[end-n+1:end]]
end

# simple little diff function
Diff(x,y) = maximum(abs.(x - y))

Diff(NonLinNA(v),NonLin(v))

#However

v = fft(randn(N))/N
v = fft(u)/N

Diff(NonLinNA(v),NonLin(v))
findall(x -> x>1e-15, abs.(NonLinNA(v) - NonLin(v)))
```

Then I use Dr. Lins function in my ODE solver and get the same thing as with my own function.



# Tuesday, January 26, 2021

8:12 AM - Today I want to put my nonlinear implementation into Dr. Lin's ODE solver.

12:32 AM - I got my function to run in Dr. Lin's ODE solver. I also discovered that the above is can be corrected by the following: (observe the difference in the `v_pad = ...` line)

```julia
## My function
pad = (3N+1)÷2
K = N + pad
NonLinNA = function (v)
    v_pad = [0; v[2:n+1]; zeros(pad);v[n+2:N]]
    nv = fft(bfft(v_pad).^2)/K
    Nv_dealiased = ℓ .* [nv[1:n+1]; nv[end-n+1:end]]
    # ifftshift(conv(fftshift(v),fftshift(v))[N-(n-1):N+n])/N
    # v_pad = [v[1:n]; zeros(pad);v[n+1:N]]
    # nv = F*(real(iF*v_pad)).^2*K/N
    # [nv[1:n]; nv[end-n+1:end]]
end
```



# Monday, February 2, 2021

9:16 AM - On Friday of last week Dr. Lin and I mad a big breakthrough in the stability of my solver. The issue with stability seemed to be that the symmetry (conjugate symmetry) in the Fourier modes was not being persevered. This had the effect of creating a pile-up (of energy) in the linear modes. I fixed this by enforcing the symmetry in the solver the ETDRK4 stepper. as shown here:

```julia
function step!(u,v)
       Nu = NonLin(u)
       @.  a  =  E2*u + Q*Nu
       Na = NonLin(a)
       @. b  =  E2*u + Q*Na
       Nb = NonLin(b)
       @. c  =  E2*a + Q*(2Nb-Nu)
       Nc = NonLin(c)
       @. V =  E*u + alpha*Nu + 2beta*(Na+Nb) + gamma*Nc
       v[:] = [0; V[2:(N-1)÷2+1]; reverse(conj.(V[2:(N-1)÷2+1]))]
end
```

However, after I did this some other things appear to have gone wrong, which I am currently investigating. For one things the autocorrelations seem to have steeper decay than that of Dr. Lin's. Another is that the energy spectrum flattens out at 1e-7 rather then 1e-37 for Dr. Lin's. The last observation is with regard to a time series of 1000 observations. The steps size is 1e-3 and we observe every 100 steps.

11:15 AM - It looks like I found the bug. it read
```julia
function NonLin(v)
    v_pad[2:n+1]        = v[2:n+1]
    v_pad[end-n+1:end]  = v[n+2:end]
    nv = Fp*(real(iFp*(v_pad)).^2)/K
    return ℓ .* [nv[1:n+1];v_pad[end-n+1:end]]
end
```
but it should have read (and now does)
```julia
    return ℓ .* [nv[1:n+1];nv[end-n+1:end]]
```
I made this change and now things seem to be working much better. It resembles the Kassam-Trefethen code in the transient part and is stable. An interesting note: when I only enforced symmetry in the NonLinear function the code was still unstable. So, I returned the symmetry enforcement to the stepper and that fixed the problem.
