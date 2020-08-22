# July 24, 2020

10:58 AM - I am working on investigating the effects of varying the parameter `Nex` in the Wiener filtering functions. The first thing to do will be to add controls for `Nex` in the functions.

* added `Nex` as function argument to `KSE_modredrun` in `KSE_moderedrun.jl` and adjusted it's call in `runme_g3.jl`

* added it to `get_wf` in `Model_reduction.jl` and it's call in `KSE_modredrun`

* `Nex` is already a parameter of `vector_wiener_filter_fft`. 

Testing work with run on the server. Earlier I tested the "adjusted" scripts which set correlations between district F. modes to zero. 

11:20 AM - Tested `Nex` dependent functions.  Test failed because of incomplete coding

11:49 AM - Corrected problems and reran test on thelio. This time it worked and the results notebook also started to work though it is going very slowly. 

12:19 PM - Ran runme_g3_nex.jl script this will run three runs of `Nex` = 512, 128, 32. It loads the solution each time so It  may be able to be run on my computer. Didn't work because I asked for `M_out = 1024` when `Nex = 512`, so I included the code `M = min(M_out, Nex)` just before the output `h_wf` is sent. 

12:27 PM - Reran `runme_g3_nex.jl` (job 104). Another problem when running the reduced model it still had `M_out > Nex` and tried to access `h_wf` at a index higher than `Nex`. So, I add the code: `M_out > Nex && ( M_out = Nex )` at the beginning.

12:35 PM - Reran `runme_g3_nex.jl` (job 105). Another problem in padding the coefficients for the spectral factorization sicne `L > Nex`. so I add the code:   
```julia
l_pad_minus = Nex >= L+1 ? cat(dims = 3,l,zeros(nu,nu,Nex - L - 1)) :  l[:,:,1:Nex] 
```
1:11 PM - Reran `runme_g3_nex.jl` (job 106). Another problem. Time to meet with Dr. Lin so I ran `runme_g3.jl` at 1:30 PM (108) hope fully I will have results to show. Didn't get results fast enough. But it looks like it went through. 

1:30 PM - Met with Dr. Lin. He started with a request that I look at how the Wiener filter changes as the `obs_gap` decreases. So, I said I would have that done by next Wednesday. The thing is that when the `obs-gap` changes the resultant Wiener filter does a slightly different thing. For instance, when `obs_gap` is 100 the wiener filter gives you a prediction in 100 steps (whatever the time equivalent of that is). So, that would have to be considered. 

Next we talked about what I was working on just now, which was adjusting the code, to allow for the z-sepectrum to be computed using fewer evenly spaced grid points on the unit circle. This lead to a discussion on how the z_spectrum was approximated. I need 

3:23 PM - Ran script `runme.jl`:
```julia
using Dates

include("Model_KSE.jl")

# Parameters for KSE model

T = 10^5 # Length (in seconds) of time of run
T_disc = Int(T/2) # Length (in seconds) of time discarded
P = 21.55  # Period
N = 96  # Number of fourier modes used
h = 1e-3 # Timestep
g = x -> randn()*cos(2π*x/P)*(randn() + sin.(2π*x/P))
obs_gap = 100

uu, vv, tt =  my_KSE_solver(T,
    T_disc  = T_disc,
    P = P,
    N = N,
    h = h,
    g = g,
    n_gap = obs_gap)

# set save destinations
sol_file = "../data/KSE_Data/KSE_sol_Lin.jld"
paramaters = Dict(
      "T" => T,
      "T_disc" => T_disc,
      "P" => P,
      "N" => N,
      "h" => h,
      "g" => "x -> cos(π*x/16)*(1 + sin.(π*x/16))",
      "obs_gap" => obs_gap,
      "tm" => now())
dat = Dict(
     "dat_uu" => uu,
     "dat_vv" => vv,
     "dat_tt" => tt)

Data = merge(paramaters,dat)
save(sol_file,Data)
```
to produce a run of Dr. Lin's data from the paper (job 109). This didn't work because I forgot to include using JLD

# Monday July 27, 2020

1:01 PM - Today I got to a late start because Erin and Daniel had doctors appointments, so me and Emeline watch Tinker Belle. Anyway, what little I have been doing this morning has been reading Chapter 23: Spectral Estimation, from 

* Pollock, David Stephen Geoffrey, Richard C. Green, and Truong Nguyen, eds. Handbook of time series analysis, signal processing, and dynamics. Elsevier, 1999.

It has been going pretty well. The problem here is that it is all done with real times-series. I think I read the Chapter through and try to understand it best I can in this setting and then I'll extend it to complex times-series myself. 

Actually, I am going to try and take a break from that to get a run of the KSE with Dr. Lin's parameters from

* Lin, Kevin K., and Fei Lu. "Data-driven model reduction, Wiener projections, and the Mori-Zwanzig formalism." arXiv preprint arXiv:1908.07725 (2019).

* Lu, Fei, Kevin K. Lin, and Alexandre J. Chorin. "Data-based stochastic model reduction for the Kuramoto--Sivashinsky equation." Physica D: Nonlinear Phenomena 340 (2017): 46-57.

1:19 - Ran `runme.jl` (I will rename this one `runme_Lin.jl`) with the following parameters
```julia
## Parameters for KSE model
T = 10^5 # Length (in seconds) of time of run
T_disc = Int(T/2) # Length (in seconds) of time discarded
P = 21.55  # Period
N = 96  # Number of fourier modes used
h = 1e-3 # Timestep
g = x -> randn()*cos(2π*x/P)*(randn() + sin.(2π*x/P))
obs_gap = 100
```
This is job 124. The email never came for some reason but the job completed anyway. 

3:15 PM - I have already downloaded the data made by the above job (`KSE_sol_Lin.jld`) and am now writing the spreadsheet to discuss it. I also read more of chapter 23, it is interesting.

# Tuesday, July 28, 2020

11:30 AM - Started late. Right now I am coding the smoothers suggested in Pollack (p. 711-12). 

12:32 - Lunch time.

1:18 - Back to work. I am working out the computation of the smoothers using the polynomials package. Worked on `z-spect` function pretty much the whole time. 

5:17 - It looks like using a periodogram with many more points and then smoothing is going to give me a lot better estimates of the z-spectrum. I will continue to pursue this tomorrow. Done for the day.

# Wednesday, July 29, 2020

10:37 AM - Starting research. I will start by reading Pollack on periodograms. 

12:17 PM - Seeing if I can reproduce the periodogram from `DSP.jl`

1:16 PM - I was able to reproduce the periodogram function from `DSP.jl`. but I don't really know what that does form me at this point. 

2:00 PM - Meeting with Dr. Lin and group. Use better approximation techniques in code.

# Thursday, July 30, 2020

2:00 PM - A good day so far I did a lot with the new spectral density estimators and am contemplating how to incorporate them into the code. I wrote a small tutorial/ reference sheet about the DFT and `fft`. There are two things I want to work on now. (1) use a periodogram analogue, some "cross-periodogram" to approximate the cross spectral density. (2) consider rational approximation of the smoothed spectral density of the predictors and preform the spectral factorization on the quotient. 

4:00 PM - Finished writing the function `z_crsspect_scalar` in `Model_reduction.jl` it gives much smoother approximations to the cross spectrum.

# Friday, July 31, 2020

12:00 PM - The goal today will be to use the smoothed periodogram to improve performance for spectral factorization. 

12:54 PM -  I realized if I just increase `Nex` and `L` considerably the approximation becomes a whole lot better. The problem I run into now however is memory. So, I will look into sparse array computation. 

I put in sparse array computation and it worked out, that is I got a computation without an error. I talked to Dr. Lin for about an hour and a half. I feel a lot like I am finally fitting into the research area. 

Anyway, I think I will work a little tomorrow. I really need to try and work harder each day. 

* Test the new Wiener filter is it stable, how does it do. 

* Think about rational approximations!!

Worked: 4 hrs.

# Monday, August 3, 2020

9:41 AM - Today I plan on testing the new and improved filter in many ways. Then I want to read up of rational approximations over the unit circle. 

To start with I will prepare a script to run this on thelio. 

11:14 AM - I just set the new code to run. It occurred to me to compute the storage it will need. I will do that now. So, I ran a few things and debugged working right on the terminal. The job I just sent was in batch. 

Remember to symmetrize autocovariance sequence of nonnegative lags: Done

I realized that the reason the process was being killed was it that I forgot to change PI ( the control parameter for which factorization method used) was not changed use the `SparseArray.jl` one. I ahev added `SparseArray.jl` to both routines now and have kept the "/" routine rather than the `pinv`. 

11:30 AM - Met with Stacy she said that I can probably expect funding in the spring, it's just that when I get funding from the graduate college, they only ever approve on a semester by semester basis. 

2:05 PM - I got the code to run. So Now I will put it on the server and see what comes out. 

3:36 PM - The server keeps killing the job. It was last run on the server at T = 10^5. Since then I was able to run it on my computer with T = 10^4. I guess I will try to run it on my machine at 10^5 while I am reading and see what happens. If it works then I will have to figure out why it doesn't run so well on the server. Then once I have done that I will need to run an experiment with par = 500, 1000, 1500, 2000 or so. 

3:41 PM - I now begin a literature review of rational approximations on the unit circle. 

5:30 PM - I am done for the day. I thought about the approximation problem. One thing I've got is that it should be a sort of least-square things since we have many more points then degrees of freedom in the approximation. 

10:15 PM - This evening the calculation ran on my computer and here is the output:
```Sol save location: Data\KSE_sol_Lin.jld
WF save location: Data\KSE_wf_Lin-Mo10000.jld
redmodrun save location: Data\KSE_rmrun_Lin-Mo10000.jld
the Parameters ===================
T : 100000
P : 21.55
short : false
h : 0.001
gen : _Lin
N : 96
loadsol : false
loadwf : false
q : 0.0:0.29156312330299705:27.69849671378472
d : 5
par : 1500
g : x -> cos(π*x/16)*(1 + sin.(π*x/16))
T_disc : 50000
tm : 2020-08-03T15:40:37.076
M_out : 10000
obs_gap : 100
p : 1500
n : 3
==================================
data saved
Get_wf computation time: 10387.237116 seconds (27.84 G allocations: 1.112 TiB, 29.33% gc time)
Wiener filter saved
Reduced Model Run Time:   0.535660 seconds (678.95 k allocations: 170.856 MiB, 8.72% gc time)
Reduced Model Run saved
```

# Tuesday, August 4, 2020

9:30 AM - The first thing I did today was run the job that ran on my computer yesterday evening (job 128). It completed just fine on my computer and should have ne reason to die. The results from last night were not good at all so today I will be investigating that. I am glade I have all the files for this on my computer. I will need to really open it up. 

Everything revolves around three scripts: `runme.jl`, `KSE_modredrun.jl`, and `Model_reduction.jl`.

1:31 PM - I really opened it up and took a look. I was analysing the autocovariance sequence of the predictors, and then wanted to inspect the difference between the "true" spectral density and the one approxiated by the Laurent polynomial (the one that is theoretically feed into the factorization function). Then I computed the factorization and multiplied them pointwise to recover the approximated  spectral density. I found that the closeness of the recovered approximation improved as I increased 200 iterations of the CKMS filter to 500, The approximation is still rather poor. When I test N_ckms = 1000, the result was very poor, wosre then all previous this was repeated. I am now trying N_ckms = 500 again. Same as before.

While these calculations are happening I am reading Zwanzig I read about the first four pages of CH. 8 and now dicided to read about statistical mechanics from Chorin

* Chorin, Alexandre Joel, and Ole H. Hald. Stochastic tools in mathematics and science. Vol. 1. New York: Springer, 2009.

Wednesday, August 5, 2020

12:38 PM - Getting ready for the research meeting. 

2:14 PM -  Today I am investigating what is going on with the CKMS filter. Currently, I a running it with par = 2000. 

I am trying to run things on the server using a jupyter notebook but it is not going so well.

# Monday, August 10, 2020

The past weekend I have been involved in the Integration Workshop. Which I think went pretty well. It was certainly very interesting. 

10:47 AM - Today I will continue investigating the wiener filtering code. Ideally I do this with a jupyter notebook on the server. The idea will be to take it nice and slow and be very, very careful. 

12:40 AM - Just finished meeting with Dr. Lin and get the impression he says too much and listens too little. I feel that the advise given could be improved if the time was taken to understand what I was doing better. At any rate much of what was said are appearently good ideas. 

Suggestions from meeting:

* Investigate stopping criterion of CKSM factorization algorithm.

* Investigate regularization of the Wiener Filter (Maybe do this for tomorrow mornings reading) 

* Try the new code on a smaller problem.

# Wednesday, August 12, 2020

Today I went back to a linear problem. I have an ODE I am working with. 

10:37 PM - Created the file 'Linear ODE Exploded View.jl'

# Tuesday, August 18, 2020

Today I am testing the code on a linear SDE. I have found a few bugs already. 

1.  in `get_wf` when the function `vector_wiener_filter_fft` is called the `sig` is not adjusted to be temporally offset with the `pred`. I put in place of `sig` as the argument `signal[1:end-1]`. This way `sig[:,i] ~ sum(psi(sig(i-1),..., psi(sig(1))`.

2.  in `z_crossspect_fft` when the function 'z_crosspect_scalar' is called the arguments given were: `sig[i,1:steps]` and `pred[j,1:steps]` This was done to ensure that the arguments would have the same length. But the function `z_crosspect_scalar` already handels input of varying length by truncation (exactly what I tried to do in the first place) and since in `z_crossspect_fft` `steps` is possibly bigger then the length `steps=nfft=nextfastfft(steps)` to let do the truncation and padding is better. 

3.  in `z_crsspect_scalar` I changed the output length from `l` to `nfft` since `nfft` is the standard length when dealing with `fft` and that is what the function `z_crssspect_fft` wants when it is called.

With these changes, it compiles.

# Wednesday, August 19, 2020

8:40 AM - Yesterday, I saved data and to day I loaded it and analyzed it. The data was a reduced model run of length 5e5 - 1e4 run on a Wiener filter computed from a run of length 1e4. The reduced run (no reduction at all actually) blow up at t = 208, it grows exponentially with a factor of 45.44.

Then I began to work on state space model generator.

2:35 AM - Meeting with Dr. Lin it was decided that I would focus on the stability of the factorization algorithm. The idea is that we need the to check and make sure the factorization algorithm is converging. So, here is the plan:

1.  read up on algorithm

2.  go back to all the old studies I did and verify the new code. This will require me to clean up my github.

Thursday, August 20, 2020

9:33 AM - I arrived at the office ready to work at 9:00. Now, I have been relearning git and am preparing to move everything over and organize it.

# Friday, August 21, 2020

10:57 AM - This morning I got it all working the plan will be to leave my laptop at the school now and do all the work I need to do at home. This will confine me to the living room, but saves me transporting my laptop. I plan to basically work at school these days, anyway. 

Now, The goal is to start from the bottom and work my way up to the top. Basically, I do what before took months in a manner of days. First, I'll start with a scalar linear process and use the software in `

I moved my research journal to the repository and converted it to markdown. Now it will be travel with the code and is easily readable on github.

12:00 PM - Now, I get to work on the code. Here is the

## Todo
### Examples
1. scalar linear
2. vector linear
3. scalar nonlinear Langevin (double well)
4. general more general state space models (with easily computed spectral densities)
5. KSE!

### For each example
1. reproduce the model
2. apply a form of model reduction
3. verify results by computing and comparing autocovariances

3:18 PM - Working on the scalar linear SDE (LSDE).

4:18 PM - Verified numerical implementation of LSDE using `Scalar LinearSDE Model Tester.ipynb`. Now I will preform model reduction on it, I will use `Scalar LinearSDE Model Reduction.ipynb`.

(I inserted some latex in `Scalar LinearSDE Model Reduction.ipynb`)

5:00 PM - Done for the day. I computed the wiener filter by hand and got exactly (1 + hA) for h_0 and 0 everywhere else. Monday, I will continue to work in this direction. 