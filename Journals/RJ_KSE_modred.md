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


# Tuesday, November 24, 2020

9:55 AM - Yesterday, I got a result for the Wiener filter using KSE data and short psi predictors. The Wiener filter was suspect, though, as it had coefficients in the thousands. Today, I want to vary the parameters a bit and see if I can get a more reasonable output for the Wiener filter. I will run these now on thelio because it takes my computer a very long time.

So, first I will rerun *Experiment Nov 23, 2020 1* on thelio.

10:31 AM - Ran "Examples/KSE/KSE_data_gen.jl" (job 157) to generate KSE data and save it as
"data/KSE_Data/KSE_sol_lin.jld". Here, `gen = "lin"`. So, now we have a copy of the data of thelio, with the standard `gen = "lin"` parameters:

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

#### Experiment Nov 30, 2020 1 (Langevin data with smallest possible `nfft`)
