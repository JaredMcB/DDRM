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

5:00 PM - Done for the day. I computed the wiener filter by hand and got exactly (1 + hA) for h_0 and 0 everywhere else. Monday, I will continue to work in this direction. This is not however what the notebook says. I think that the autocorrelation time for my example is too, big so I don't have enought independent points.


# Monday, August 24, 2020

11:04 AM - Last Friday I computed the WF by hand and got exactly what makes sense. I will start there today and open up the program to make sure I get this result numerically. Goal: have this working by noon.  

The particular case I start with is
```julia
A = reshape([-0.05],1,1)
σ = reshape([1],1,1)
Xo = [1]
t_disc = 1000
gap = 1
scheme = "EM"

t_start = 0
t_stop  = 1e4
h       = 1e-2
```

11:38 AM - I rearranged this a trifle, I moved the functions
* `_crosscov_con`
* `_crosscov_dot`
* `my_crosscov`
* `my_crosscor`
* `z_crossspect_fft`
* `z_spect_scalar`
* `z_crossspect_scalar`

from `Model_Reduction_Dev.jl` to `AnalysisToolbox.jl`.
I also updated `auto_times` in the `AnalysisToolbox.jl`. There was a problem with me lagging the length of the series. Inserting the  "- 1" here was required, in `L = minimum([lx - 1, 10^6])`.

3:00 PM - I continue to investigate the LSDE model. It seems that the number of effective samples was a role to play in the accuracy of the wiener filter. So, I will try to quantify this.

3:58 PM - It is apparent that this issue I began the day investigating is a problem with the effective number of samples. When A = -0.05 the autocorrelation time I computed with `auto_times` was about 20 seconds. run for 9,000 seconds this made the effective number of samples about 500, which would seem unacceptable. so I ran it with `t_stop = 1e5` and went from `h_wf[1,1,1] = -.56` to `h_wf[1,1,1] = -1.2` which is better. (the WF took 81 secs to compute with a series of size 9,920,232).

I started a new session and now what to try `t_stop = 1e6` but with `gap = 10` this should take up the about the same memory but stretch further, as in have about 10 times more effective samples. This got the wiener filter even closer to the solution.  Now I try `gap = 100` and `t_stop = 1e7`

I think I want to give thelio a job. I would like to investigate the convergence of the Wiener filter.

4:49 PM - I am done for the day. Tomorrow I will see how these wiener filters run.


# Tuesday, August 25, 2020

8:45 AM - Started research for the day. The goal for today in the roughly three hours I have to research will be to learn all I need to from the scalar linear SDE experiments.

Yesterday, I found that the effective number of samples had a effect on the accuracy of the WF. I would like to study that effect and (1) find how the effective number of samples scales with the errors. (2) I can compute the autocorrelation time by hand, I would like to see how it compares to what I have estimated. (3) Then I would like to run the reduced models, perhaps with an optimal noise (?).

8:54 AM - How does the noise scale with the effective number of samples?

### Parameters:
```julia
A       = reshape([-0.5],1,1)
σ       = reshape([1],1,1)
Xo      = [1]
t_disc  = 1000
gap     = 10
scheme  = "EM"
d       = size(A,1)
t_start = 0
t_stop  = 1e6
h       = 1e-2
Δt      = h*gap
M_out   = 100
```

### Experiment Table
| `A` | `t_stop` |`gap` | `N_eff`     | First comp err | Tail err |
|---  |---       |---   |---          |---             |---       |
|-0.5 | 1e7      | 10   | ≃ 4.98e6   | ≃ 4.97e-3     | 4.57e-5   |
|-0.5 | 1e6      | 10   | ≃ 4.9?e5   | ≃ 5.95e-3     | 4.04e-5   |
|-0.5 | 1e6      | 10   | ≃ 4.94e5   | ≃ 4.94e-3     | 4.04e-5   |
|-0.5 | 1e5      | 10   | ≃ 5.14e4   | ≃ 5.47e-3     | 7.41e-4   |
|-0.5 | 1e4      | 10   | ≃ 4.34e3   | ≃ 5.46e-2     | 9.24e-2   |


3:19 PM - Back to work.

4:29 AM - Worked on the thelio script. Had to debug a little.


# Wednesday, August 26, 2020

11:49 AM - Got started.

12:23 PM - Ready to run thelio_runme.jl on thelio. I cloned my github repo on to thelio.
Here is thelio_runme.jl:
```julia
include("modgen_LSDE.jl")
include("../../Tools/Model_Reduction_Dev.jl")

# include("modgen_LSDE.jl")
# include("..\\..\\Tools\\Model_Reduction_Dev.jl")

using JLD
using PyPlot
using DSP: nextfastfft

function runner(;
    A       = reshape([-0.5],1,1),
    σ       = reshape([1],1,1),
    Xo      = [1],
    t_disc  = 1000,
    gap     = 10,
    scheme  = "EM",
    d       = size(A,1),
    t_start = 0,
    t_stop  = 1e6,
    h       = 1e-2,
    Δt      = h*gap,
    M_out   = 100)
    println("===========t_stop = $t_stop===========")
    println("Time to get data: ")
    @time X = modgen_LSDE(t_start,t_stop,h,
        A = A,
        σ = σ,
        Xo = Xo,
        t_disc = t_disc,
        gap = gap,
        scheme = scheme)

    N       = size(X,2)
    nfft    = nextfastfft(N)
    X = [X zeros(d,nfft - N)]

    τ_exp, τ_int    = auto_times(X[:])*Δt
    N_eff           = N*Δt/τ_int

    println("Time to get h_wf: ")
    Psi(x) = x
    @time h_wf_num = get_wf(X,Psi, M_out = M_out)

    h_wf_ana = zeros(1,1,M_out)
    h_wf_ana[1,1,1] = (1 .+ h*A)[1]

    err     = abs.(h_wf_ana - h_wf_num)
    comp1_err = err[1,1,1]
    tail_err = maximum(err[:,:,2:end])
    err_sum = sum(err)

    println("======================")
    println("N_eff : $N_eff")
    println("comp1_err : $comp1_err")
    println("tail_err : $tail_err")
    println("======================")
    [N_eff; comp1_err; tail_err]
end


T_stop = map(x -> 10^x, 4:.5:7)
data = zeros(3,length(T_stop))
M = 20

for i in 1:length(T_stop)
    for j in 1:M
        dat = runner(t_stop = T_stop[i])
        data[:,i] .+=  dat/M
    end
end

save("/u5/jaredm/data/LSDE_Data/NoiseVNeff.jld", "data", data)
# save("c:\\Users\\JaredMcBride\\Desktop\\DDMR\\Examples\\LinearSDE\\LSDE_Data\\NoiseVNeff.jld", "data", data)
```

Before I run it I need to analysis memory cost. So, the longest time series will be for `t_stop = 1e7` which since `gap = 10` and `t-disc = 1000` will be of size ≃1x1x999900, so about 8 MB. Not too much. Job 131.

Immediate error:
```ERROR: LoadError: InitError: PyError (PyImport_ImportModule

The Python package matplotlib could not be found by pyimport. Usually this means
that you did not install matplotlib in the Python version being used by PyCall.

PyCall is currently configured to use the Python version at:

/usr/bin/python3
```
I did not want to address it just now so I simply commented out `using PyPlot` and sent it back. Job 132.

2:22 PM - Job 132 still running, I think `M = 20` is what's taking so long.

Now to finally investigate the convergence of the factorization algorithm.

3:52 PM - I did a convergence study of the CKMS factorization algorithm. The main idea is the following

```julia
NN_ckms = map(x-> floor(Int, 10^x),2:.25:4)
LL = complex(zeros(1,1,L+1,length(NN_ckms)))
for i in 1:length(NN_ckms)
    LL[:,:,:,i] = spectfact_matrix_CKMS(A_v2_smoothed,
        N_ckms = NN_ckms[i])
end

Norm = map(i -> norm(LL[1,1,:,9] .- LL[1,1,:,i],Inf),1:9)
loglog(NN_ckms, Norm)
```

`NN_ckms` is the array of the various `N_ckms` used in the study. `LL` holds the factored coefficients. So, I did my convergence study using the data for the KSE model, the parameters were the same used in Dr. Lin's work. I used the modes 2,3,4 and found that the algorithm converged at about machine epsilon at 5.5e3.
in one dimension it was around 1.5e3. Here is another look:

```julia
pred = vv[3:6,:]
nu, stepsx = size(pred)

steps = stepsx
nfft = nextfastfft(stepsx)
nffth = Int(floor(nfft/2))

L = 1500
lags = -L:L

# Smoothed viewing window
lam = _window(L, win = "Par", two_sided = false)

R_pred_smoothed = zeros(Complex,nu,nu,length(0:L))
for i = 1 : nu
    for j = 1 : nu
        temp = my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
        temp = .5*(temp[L+1:end] + conj(reverse(temp[1:L+1])))
        R_pred_smoothed[i,j,:] = lam .* temp
    end
end
R_pred_smoothed
plot((0:L)*Δt,real(R_pred_smoothed[2,2,:]))

NN_ckms = map(x-> floor(Int, 10^x),2:.25:5)
LL = complex(zeros(nu,nu,L+1,length(NN_ckms)))
for i in 1:length(NN_ckms)
    LL[:,:,:,i] = spectfact_matrix_CKMS(R_pred_smoothed,
        N_ckms = NN_ckms[i])
end

Norm = zeros(9,0)
for i = 1:nu
    for j = 1:nu
        global Norm
        Norm = [Norm map(k -> norm(LL[i,j,:,9] .- LL[i,j,:,k],Inf),1:9)]
    end
end


loglog(NN_ckms, Norm)
```
The result is saved in `Figure8-26-2020.png` in the data folder in tools.

4:40 PM - Make the eps, which is in the stopping criterion it is a tolerance of the ∞-norm of the consecutive outputs, make this a tunable parameter. The paper suggested that `N_ckms = 200` usually work should work. I think that that is true when the number of terms is small. These
jobs, I am using it for have very many parameters. So, here is the table

| d = 1 | d = 3| d = 4 |
|---    |---   |---    |
|≃1.5e3 |≃5.65e3|≃5.65e3|

A finer grid of number of iterations may be helpful. Either way they all converged by 1e4.


# Thursday, August 27, 2020

8:36 AM - Getting started. Yesterday, I investigated the convergence of the CKMS factorization algorithm, and noted that the more variables involved, the more degrees of freedom, the longer the algorithm takes to converge. So, I think a stopping criterion should be imposed. So, I add that to today's agenda:

### To do:
* Add stopping criterion to CKMS implementation.
* Analyze the error message for job 132, fix problem and run it again. Continue with the convergence analysis of the WF. Then run the reduced models.

8:54 AM - Add stopping criterion with tunable parameter eps to CKMS implementation.

```julia
i = 0
errK = errR = 1
Err = zeros(0,2)
while (errK > ϵ || errR > ϵ) && i <= 10^5
    hL = h*L; FL = F*L

    # Stopping criteria stuff
    i += 1
    FL_RrhLt = FL/Rr*hL'
    hL_RrhLt = hL/Rr*hL'
    errK = norm(FL_RrhLt)
    errR = norm(hL_RrhLt)
    Err = [Err; errK errR]
    println("err : $errK and $errR")
    #

    K_new = K - FL_RrhLt
    L_new = FL - K/Re*hL
    Re_new = Re - hL_RrhLt
    Rr_new = Rr - hL'/Re*hL

    K = K_new
    L = L_new
    Re = Re_new
    Rr = Rr_new
end

println("i : $i")
k = K/Re
re = Re
```
This worked well and I was able to see that for the example the `Re` converged faster than the `K` did, that is, the estimated covariance of the noise converged faster than the Kalman gain vector. So, I think convergence of the Kalman matrix is what I want to use for the stopping criterion, of course I could check the other things but they do not directly contribute to the output as can be seem in the last two lines of the code snippet.

11:50 AM - I would like to test the factorization on a random function.

3:54 PM - I have set it up as follows:

The normal example converged very differently.

The `fft` gives coefficients but times `nfft` and `v_1` is the coefficient to `z/nfft`

4:56 PM - I really need to investigate this because CKMS is not working the way I think it should.


# Friday, August 28, 2020

1:57 PM - I am working at home this afternoon. The goal for today will be to update and correct the code for the CKMS implementation. I first noticed that something was internally wrong yesterday, when I tried to investigate the convergence of the algorithm. I included a dynamic stopping condition (one that depended on the difference of consecutive outputs). I began to realize that it was not getting out the right answer. Indeed, I put in the classic example

```julia
P = zeros(2,2,2)
P[:,:,1] = [6 22; 22 84]
P[:,:,2] = [2 7; 11 38]
```

The output should be
```julia
l[:,:,1] = [2 1; 7 3]
l[:,:,1] = [1 0; 5 1]
```
 But I did not get this out. Today I went back the earlier code:
 `spectfact_matrix_CKMS` from the file `DDMR\\Tools\\Model_Reduction_Dev.jl`.

 I put in the following
 ```julia
P = zeros(1,1,2)
P[1,1,1] = 10
P[1,1,2] = 3

L = spectfact_matrix_CKMS(P)
 ```
and got
```1×1×2 Array{Complex{Float64},3}:
[:, :, 1] =
 3.000000000000001 + 0.0im

[:, :, 2] =
 0.9999999999999998 + 0.0im
```
which is correct.


# Monday, August 30, 2020

8:00 AM - I started by read all of last week.

8:23 AM - It turned out, from what I saw Friday, that the function `spectfact_matrix_CKMS_SC` ('SC' stands for 'stopping condition') did give a correct answer but it was not the one found in the paper below:

G. Janashia, E. Lagvilava and L. Ephremidze, "A New Method of Matrix Spectral Factorization," in *IEEE Transactions on Information Theory*, vol. 57, no. 4, pp. 2318-2326, April 2011, doi: 10.1109/TIT.2011.2112233.

But these factors can be off by a unitary matrix. However, strictly speaking does not seem to be the case here, though it is close. I plotted the difference of the recovered spectrum, that is, I did this

```julia
P = zeros(2,2,2)
P[:,:,1] = [6 22; 22 84]
P[:,:,2] = [2 7; 11 38]

l1 = spectfact_matrix_CKMS_SC(P)[1]

l = zeros(2,2,2)
l[:,:,1] = [2 1; 7 3]
l[:,:,2] = [1 0; 5 1]
l2 = l

ll = size(l1,3)

S1_fun_minus(z) = sum(l1[:,:,i]*z^(-i+1) for i = 1:ll)
S2_fun_minus(z) = sum(l2[:,:,i]*z^(-i+1) for i = 1:ll)

res(z) = S1_fun_minus(z)*S1_fun_minus(z^(-1))' -
            S2_fun_minus(z)*S2_fun_minus(z^(-1))'

d= 2; nfft = 10^3
Res = complex(zeros(d,d,nfft))
for i = 1:nfft
    Res[:,:,i] = res(exp(im*2π*i/nfft))
end

for i = 1:d
  for j = i:d
    plot(abs.(Res[i,j,:]),label = "($i,$j)")
  end
end
legend()
```

12:01 PM - Everything's to be working as it should. I don't understand why the JLE example has so strong a difference between the CKMS solution and the analytic solution. Though this difference decreases and the number of iterates increases.

Now, I want to test higher order function (much higher order functions).


# Tuesday, September 1, 2020

8:17 AM - Today I want to clean up and focus the research I have done so far and draw some conclusions. That will be the goal of this morning summarize what I have been up to the past few working days, see what I have learned and see what more I need to learn before I can start working on the Wiener Filter again.

#### What I need:
I need to see how well CKMS factorization preforms with inputs of 1500 to 2000 coefficients  and of various sizes. I need to know how many iterations are necessary for effective wiener filtering. Right now my metric is qualitative, how well does the graph match.  

#### What have I done?
The project is all contained in `DDMR\Tools\CKMS_filter`. Here we find the files
* `AnalysisToolbox_scratch_ckms.jl` - this contains the CKMS implementation with the stopping criterion called `spectfact_matrix_CKMS_SC`         
* `CKMS with Approximated Spectra.ipyn` - this notebook was created yesterday to help me investigate the performance of CKMS when feed `fft` generated coefficients.
* `CKMS with Approximated Spectra.jl` - just a Julia file of the notebook
* `Basic functions.ipynb` - this is a notebook containing hither to successful tests of CKMS using natural coefficients (as opposed to `fft` generated coefficients)                   
* `Basic functions.jl` - just a Julia file of the notebook

In `Basic functions` I test the following cases:
|Type                        | No. of Terms|
|---                         |---          |
|scalar                                 | 2|
|2×2-diagonal                           | 2|
|2×2-nondiagonal                        | 2|
|2×2-nondiagonal (w/ unit circle zeros) | 2|
|d×d-matrix valued random               | m|

The above all worked very well visually.

The next thing to do is to look at `fft` generated autocovariances.


# Wednesday, September 2, 2020

9:09 AM - Yesterday, the CKMS filter was shown to behave well for certain examples, at least visually. Today I will investigate who it preforms on the KSE data. It preforms much better than it did when `L=55`, that is to say the plot of what I am taking to be the true spectral density is much closer to the recovered spectral factorization. However, there still seems to be a great disparity in the spectral density near 1 where it off by as much as a facto of ten. Anyway, lets see how it works now in the wiener filter. I think one problem was the number of iterates was too small.

So, now I keep all I did with CKMS very close and go back to the Wiener filter algorithm as a whole.

2:57 PM - I just concluded a meeting with the UQ Group. We mainly talked about estimation.

4:30 PM - I went to the ABD2PHD seminar. Here are some take aways:
* Find out how the graduate coordinator is.


# Thursday, September 3, 2020

9:08 AM - This morning I am working with the file `Scalar LinearSDE_ExplodedView.jl` which is in the folder `Examples\LinearSDE` the point is to return to the Wiener filter having a better understanding of the performance of CKMS factorization algorithm. the goal for the day is to get the code to return an approximation of the wiener filter (we know precisely what the exact WF is in this case) which is good enough to run a stable model.

# Friday, September 4, 2020

9:31 AM - Yesterday, I had some computer troubles. Today, I have the same goal. I control the factorization and cross-spectral estimation parameters to produce a good wiener filter. Something I did get from yesterday is that sans noise the reduced model did blow up. I am investigating that more now.

2:40 PM - I ran the LSDE and computed the Wiener filter, the filter was unstable.

3:35 PM - Also for some reason my machine is running julia very, very slowly. I think I will restart the kernel. Restart the kernel helped. Also, I took do the size of the series I was computing. They don't need to be the big for current purposes.
3:49 PM - I noticed there was a big difference between the analytical and numerical cross spectra I am investigating this now. I found an Error in the code. Not in the functions themselves but in the ancillary code in the file `Scalar LinearSDE_ExplodedView.jl`
The result of the error was that `sig` and `pred` rather than being off set by an index, as should be the case (`pred[:,n]` = `sig[:,n-1]`, since `Psi` is identity). The next problem was I reduced the stop time to 1e4 but and the real and imaginary parts of the crspect were matching up, so I increased the stop time to 1e5 and now these plots match beautifully. So, that is I went from roughly 4.99e3 effective samples to 5.01e4.

The early problem with not mismatching the `sig` and `pred` makes sense because the output of the wiener filter was very close to 1 in the first component and zero after. Which is just what we would expect with no shift at all. There is another problem because for the following data

```julia
A = reshape([-0.05],1,1)
σ = reshape([1],1,1)
Xo = [1]
t_disc = 1000
gap = 10
scheme = "EM"

t_start = 0
t_stop  = 1e6
h       = 1e-2

Δt      = h*gap
M_out   = 100

X = modgen_LSDE(t_start,t_stop,h,
    A = A,
    σ = σ,
    Xo = Xo,
    t_disc = t_disc,
    gap = gap,
    scheme = scheme)

d, N = size(X)

nfft = nextfastfft(N)
X = [X zeros(d,nfft-N)]


τ_exp, τ_int    = auto_times(X[:])*Δt
N_eff           = N*Δt/τ_int

Psi(x) = x

@time h_wf = get_wf(X,Psi);
```

I got the following Wiener filter.

```
1×1×20 Array{Float64,3}:
[:, :, 1] =
 0.0016677863376667796

[:, :, 2] =
 0.9995543511979313

[:, :, 3] =
 -0.0004237206933257932

...

[:, :, 18] =
 -0.00031437525459553626

[:, :, 19] =
 -0.00031135831703248747

[:, :, 20] =
 -0.0003042347678346185
```



# Tuesday, September 8, 2020

3:15 PM - Today I will look closely at the problem that I saw on Friday. First thing reproduce results. Then dive in. The strangeness of the solution, why it does not match what I expect, is because I am expecting the wrong thing. There is an observation gap that may be throwing things off.

3:24 PM - I have reproduced the "strange" solution.


# Wednesday, September 9, 2020

8:16 AM - Yesterday I was unable to determine what is wrong with the solution. I will run a few more tests real quick.

8:16 AM - I found a BIG bug in `get_wf` I don't know how long it was there. When it calls `vector_wiener_filter_fft` it now does it as follows:

```julia
h_wf = vector_wiener_filter_fft(sig, pred, M_out,
            n = n, p = p, par = par, PI = PI, rtol = rtol)

h_wf = rl ? real(h_wf) : h_wf
```

Before the I changed it just now, the `sig` and the `pred` were switched, these are non-kw-arguments, so that made a difference! This mistake must have occurred as I was removing type assertion for the function. I must have copied and pasted it, from a function that had these in a different order. With this fix accomplished the LSDE reproduced model now seems to be stable.

9:58 AM - I now put in the correct error term `sqrt(h)*σ*randn(d)` and find the statistics of the reproduced model.

4:17 PM - This worked pretty well, though not great. What worked better was when there were only 2 WF coefficint. Now that the code seems to be working the next step I think is tuning. So, tomorrow morning I will write out a script for thelio.


# Thursday, September 10, 2020

1:08 PM - So, it looks like the code is ready for parameter tuning. That is I feel the code is preforming as it should.  So, I would like to run it on a few examples.


# Friday, September 11, 2020

4:47 PM - Today the goal has been to get a big batch of scalar linear SDE's running do I can be sure the code works. And I have just barely done that. I am starting with a little thing that ran comfortably on Dell, this job is running now (job 134). This ran successfully so I went to the big job. And so far so good. I may need to work tomorrow to do analysis this and maybe start the higher dimensional linear model.

The thing to do will be to get the jupyter notebook up and running.  

###Experiment Table: Linear SDE
|gen| A | h   | `t_stop` |`gap`| `MM_out`           | `N_eff`|
|---|---|---  |---       |---  |---                 |---   |
|1  |-.5|0.01 | 1e5      | 10  |[2,4,6,10,15,20,100]|≃ 50k |


# Monday, September 14, 2020

8:14 AM - Today, I need to look over the data thelio made. So, I will open a notebook in thelio to do this I, in my thelio terminal type:

```
jaredm@thelio:~$ /usr/bin/jupyter notebook --port=8080 --no-browser
```

then in my windows command line I type

```
C:\Users\JaredMcBride>ssh -L 8080:Localhost:8080 jaredm@thelio.math.arizona.edu
```

Then I put the link from putty into a browser, and voil\`{a}!

8:54 AM - I tried to get a plot by PyPlot.jl but it gave me an error. Here is a smaller preproducable example

```julia
A = rand(10)
plot(A)
```
gave this
```
ArgumentError: hasproperty of NULL PyObject

Stacktrace:
 [1] pyhasproperty(::PyCall.PyObject, ::String) at /u5/jaredm/.julia/packages/PyCall/zqDXB/src/PyCall.jl:348
 [2] hasproperty at /u5/jaredm/.julia/packages/PyCall/zqDXB/src/PyCall.jl:354 [inlined]
 [3] plot(::Array{Float64,1}; kws::Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}}) at /u5/jaredm/.julia/packages/PyPlot/XHEG0/src/PyPlot.jl:174
 [4] plot(::Array{Float64,1}) at /u5/jaredm/.julia/packages/PyPlot/XHEG0/src/PyPlot.jl:174
 [5] top-level scope at In[14]:2
 [6] include_string(::Function, ::Module, ::String, ::String) at ./loading.jl:1091
```

So, I just switched to
```Julia
using Plots; pyplot()
```
and restarted the kernel and that worked fine. Will need to get `PyPlot` working later though!

###Experiment Table: Linear SDE
|gen| A | h   | `t_stop` |`gap`| `MM_out`           | `N_eff`|
|---|---|---  |---       |---  |---                 |---     |
|1  |-.5|0.01 | 1e5      | 10  |[2,4,6,10,15,20,100]|≃ 50k   |

10:32 AM - Ran `gen=2` scritp on thelio here the param`t_stop` is small but I run `reps = 10` runs of each of the `regs = 7`.


# Thursday, September 17

1:15 PM - Today I worked on class Math 122B from 7:30 - 12:00 this morning. The top priority for the day is to work on getting data for and write up my brown bag presentation. I start by running on my computer here (home comp) the 2-d linear system (reproduction).

Looking at the data


# Wednesday, September 23, 2020

Today I spent the morning on Math 122B and am now working on research.

12:15 PM - I read a lot of:

Zouzias, Anastasios, and Nikolaos M. Freris. *"Randomized extended Kaczmarz for solving least squares."* SIAM Journal on Matrix Analysis and Applications 34.2 (2013): 773-793.

## UQ Meeting Notes
First things to do for next time:
* In order to more deeply investigate the reliability of my implementation of the CKMS factorization algorithm, I will investigate the one step prediction error and test if that is orthogonal to the predictors. I will do this for all the models I've looked at so far especially the KSE.
* Test the WF for higher order ARMA models. *done: Sept 30, 2020*


2:54 PM - Goal for the day: **Clean up examples and run them with new code**
First I will try to get Julia 1.5 on Jupyter notebook

4:01 PM - Applying the new code to the old wiener filter tests. These are all in
`Tools\Wiener Filtering\Scalar Wiener Filter`.

The first one I looked at was `Scalar Wiener filter Test AR2N Sig add WN.ipynb`.
For this run I had: `r1, r2  = .99,-.5` Observe the zero very near the unit circle.
the parameter I played with was `par`

Here's a little table
| `par` | n | p  | Quantitative result |
|---    |---|--- |---                  |
|500    | 3 |1500|poorer that 'wiener_filter_fft.jl'|
|1500   | 3 |1500| better |
|2000   | 3 |1500| better still|
|1500   | 2 |1500| same|
|1500   | 3 |500| same|
|1500   | 2 |500| same|


# Sept 29, 2020

2:35 PM - listened in on RTG meeting. Very interesting.

Today, I continue to clean up and run basic experiments.

3:03 PM - The goal for the rest of the day is to write up a generator for general ARMA(p,q) processes. Compute it's wiener filter (tough).

3:27 PM - I already wrote cod for the first bit I think I will dispense with writing code for the analytic solution and just test that the predictors and orthogonal to the errors.


# September 30, 2020

12:29 PM - This morning in have been testing the Wiener filtering code on general ARMA processes. And it has been working well so far. I am convince the code is working right. Issues in the performance of the WF are likely due to insufficient effective samples. Since I am using randomly chosen values the autocorrelation time of the process can get long.

Now, It is time to go to us this tool (covariance of observations and errors) to test performance in model reduction settings

12:38 PM - *Nonlinear Langevin* I will do one example here and then go to KSE.

**Note:** There is a problem with the rejection sampling algorithm. I will need to get back to that later.

# October 1, 2020

11:32 PM - Today the task will be to get the Wiener filtering code to work for the Double-welled, overdamped Langevin equation (DWOL) to work.  

I am a bit suspicious of the time series generator so I will first verify that against the distribution derived by Fokker-Planck

11:52 AM - I checked this and the data generator does produce a series that has the desired distribution.

2:30 PM - Code is not working I found a little bug in the `dot` function from `LinearAlgebra.jl` small reproducable example is

```julia
x = randn(10^6) + im*rand(10^6)
dot(x,x)
```

That seems to immediately kill julia.

### Meeting with Dr. Lin
* We discussed the issue with `dot` form `LinearAlgebra.jl` and how it seems to be a installation error as it only seems to occur on my machine.
* The WF that the code is giving me is way off there must be a bug in the code. So, I looked at the estimated spectral density of the predictors first as recovered by multiplying the factorization and then by direct approximation using a smooth periodogram. There was a great discrepancy between the two.
* Dr. Lin provided a big picture plan for the very near future. That is here an technique to implement:
  1. we run Wiener Filtering code to get the approximate reduced model and use rational approximation to improve that as much as possible.
  2. Then use that as the initial estimate for a time domain optimization algorithm where the loss function has measures the kth step prediction error with noise. The parameters and form of the noise is approximated using the one step predication error given by the approximate initial WF.
  3. We use the time domain optimization step to improve the spectrally derived estimate.

Yesterday, there were some very strange plots as I was comparing the recovers spectral density of the predictors and the spectral density of the predictors approximated by the smoothed periodogram. I have discovered why the plots were in fact so strange. First, I recovered the spectral density incorrectly. Rather than S_x(θ) = S_x^-(θ)S_x^+(θ) = S_x^-(θ)S_x^{-\*}(θ) as it should have been I had S_x(θ) = S_x^{-\*}(θ)S_x^+(θ) = S_x^{-\*}(θ)S_x^{-\*}(θ). The periodogram was also in error, now this doesn't bare directly on the performance of the WF since the spectra density of the predictors is never actually estimate. It does make me suspect the function `z_crossspect_fft` from `Model_reduction_dev.jl` which is used in the production of the WF to compute the cross spectral density of the signals and the predictors.

The problem is in the smoothing.


# October 7, 2020

3:00 PM - Getting back into it. I am investigating how the current code preforms on the Linear SDE. The problem does seem to be in the smoothing parameters of the periodogram. I sometimes the `_smoother` function from `AnalysisToolbox.jl` would give bad smoother function (have negative values, an artifact of the `Polynomial.jl` package I suspect) this occurs only when I am asking for too much smoothing anyway.

I wrote a little smoother plotter function.
```julia
function smoother_plot(n,p,ty)
    μ = _smoother(n,p; ty)
    round(sum(μ);digits = 5) == 1.0 || println("bad smoother")
    plot(-n*p:n*p,μ)
end
```


# October 8, 2020

11:45 AM - I am still working on m=understand the affects of smoothing (both in the autocov of preds and the crossspectrum).

### Meeting with Dr. Lin
Some great suggestions:
* Automate the testing process. Ideally after I get the code working for one example some modifications may be necessary to run it on a different example. After those modifications are mad it is crucial that the code still works for the first example. So testing the old examples frequently is an important thing to do. So, I think what I will do is (1) make a list of examples that I feel are important for the code to work on:
  1. Linear SDE's ( of various dimensions and various numerical schemes)
  2. General ARMA and VRMA processes
  3. General state spaces models
  4. KSE
  5. Lorenz ('63 and '96)

  Then I'll make a short notebook for each in which I call a model simulator, run my model reduction algorithm, then generate the reduced model and compare various aspects of them with the original timeseries. Aspects to compare include
  - Autocorrelation
  - Spectral Density
  - Are the errors orthogonal to the predictors?
  - distribution of the marginals
* We talked about how the old code is working but the new code is not. By working I mean with regard to the DWOL. So, he suggested I swap out some of the new parts for the old and see if I can get it working again. This may be facilitated by creating a branch in git.

4:11 PM - Now, I'm going to work on teaching duties and write the code testers tomorrow when I start researching at 12:30.


# October 9, 2020

12:48 PM - Today I will get the automated testers ready. This should not take long and I expect to be done with the Linear SDE's by 2:00 PM (at least for the forward Euler scheme, the other schemes will not be hard to introduce later.) After that I will create a tester notebook for the DWOL model. Then I will go from there.

4:25 PM - worked on LSDE tester. It was a nice exercise. One thing I learned was that when I changed `get_wf` to include the functionality of outputting the predictor process I didn't realize that it was **not offset** with the input signal as it ought to be in our implementation of Data-drive model reduction. The Wiener filter I have programed does not build in the offset. It creates the estimator straight across. Including `X_n` in the causal prediction of `Y_n`. Anyway, I have not yet finished the tester notebook yet. Now I am going to work on my CV which is due Monday.


# October 13, 2020

1:33 PM - Today I want to finally discover what can account for the poor performance of the of Wiener filter with the new code. I will do as Dr. Lin suggested and interchange the algorithms for the various estimators.

Then after that has been accomplished and I can account for the error and fix it, I will continue to write the tester notebooks. And start to cover some good ground. The goal is still to submit a paper by Thanksgivingg.

4:44 PM - I set up two notebooks that ran the new and old code respectively the old code out preformed the new code by a lot (said differently the new code did not produce a Wiener filter while the old code did). So, I decided to open these up and look under the hood.

5:27 PM - Done for the day. Here is my report. I changed the new code (stuff in `Model_Reduction_Dev.jl`) very slightly to make `nfft` a tunable parameter (as it was in the old code under the guise of `Nex`). Then I moved the old cross spectral density estimator into `AnalysisToolbox.jl` and modified `vector_wiener_filter_fft` (in new code) to include this old estimator. I ran it and it worked. More on this tomorrow.

One more thing I want to add. I had trouble pulling the remote repository to thelio because I had change stuff on thelio and on my machine with out committing first. Anyway, there were a number of merge conflicts and it would let me pull to thelio. I didn't care about what was on thelio I don't use that for development just testing. Development always occurs on my machine. So, I just wanted the things on the remote repository to overwrite any conflict on thelio. As it turns out the thing to do was the following: (on thelio)

```
jaredm@thelio:~/DDMR$ git fetch
jaredm@thelio:~/DDMR$ git reset --hard origin/master
HEAD is now at 6d20134 Done for the day.
```


# October 14, 2020

1:01 PM - Today the goals today are as follows:
1. Does the current code ( the "new" code with the "old" xspectral estimate) work for Langevin (DWOL)?
2. Set up the poles analysis. What is happening with the poles.
3. How can I improve the xspectral estimator?

1:29 PM - code seems to be working for DWOL again. It is using the "old" xspectral estimator which I will refer to as the direct estimator as opposed to the smoothed periodogram.

## UQ Meeting Notes

Today I showed how the direct estimator when used in the wiener filtering function improved the overall preformance of the function. What I had done is described above. Now, this caused a little concerna dn so Dr. Lin has suggested I compare side by side these two estimators. They are stand alone functions so I can put in some time series of know spectral density, such as white noise and ARMA processes, as compare the out put. So, that is what I would like to do first. Here is the list.
1. Compare performance of the xspectral estimators. (This also goal 3 from today above.)
2. Run a study varying the parameters `nfft` the `par` to get a feel for what they can do.


3:14 PM - Starting on comparing the old and new spectral estimators.
#### Experiment:
The experiment is housed in Jupyter on my laptop. here is the code used. First we have
```julia
using PyPlot
include("../AnalysisToolbox.jl")

steps = 10^6
W = randn(1,steps)
```
then after running this one I run:
```julia
L    = 50
Nex  = 2^10
win  = "Par"

nfft = 0
n    = 3
p    = 1500
ty   = "ave"

spect_D  = z_crossspect_fft_old(W,W;
    L, Nex, win);
F_D = 2π*(0:Nex-1)/Nex

spect_SP = z_crossspect_fft(W, W;
    nfft, n, p, ty) ;
F_SP = 2π*(0:steps-1)/steps;

μ = _smoother(n,p;ty)
println("sum of μ: ",sum(μ))

plot(F_SP,spect_SP[1,1,:])
plot(F_D,spect_D[1,1,:])

plot([0, 2π],[1, 1])
```

5:07 PM - The results proved difficult to investigate I think I need to plot a bunch at a time. Anyway, that is what I have been doing.


# October 15, 2020

1:51 PM - I had a few thoughts last night.
1. The first was to subtract the analytic cross spectral density from the estimated one so I can see the noise better, this gives a visualization of the absolute error rather than relative error. This morning I thought about analytically computing the variance of the error of each estimator.
2. See how this error scales with sample size.
3. Vary the time step of the DWOL and see if the Wiener filter matches it. This is just what Dr. Lin suggested yesterday.
4. See how the err in the Wiener filter scales with effective sample size.

As for now, I will work on the first, so that I can compare the estimators.

### Meeting with Dr. Lin

In our meeting today we looked at the performance of the spectral estimators. It is inexplicable why the smoothed periodogram did not do well. I am thinking about setting a functionality that allows me to chose the which spectral estimator I want. to use. Either way It looks like I can get this thing to work now. Which much greater confidence.

#### Experiment
```julia
t_stop = 10^4

Random.seed!(2015)
X = DataGen_DWOL(
    steps;
    scheme, t_start, t_stop, discard,
    sig_init , sigma, V_prime,
    SM1, Obs_noise, d
    )

# Model reduction Parameters
M_out = 100
n = 2
p = 500
par = 5000 # 55
ty = "bin"
rl = true
Preds = true
PI = false
rtol = 1e-6

nfft = par*10 # 1024

@time h_wf, pred = get_wf(X, Psi;
    M_out, n, p, par, ty, nfft, rl, Preds, PI, rtol);

X_sig = X[:,2:end];

h_wf
```
Conpare these at large. as it is the commented values seemed to be beter because the tail decayed quick.
```1×2×100 Array{Float64,3}:
[:, :, 1] =
 1.01214  -0.0104821

[:, :, 2] =
 -0.00577268  0.00140345

[:, :, 3] =
 0.00251975  -0.0003896

...

[:, :, 98] =
 1.9403e-6  -4.75491e-7

[:, :, 99] =
 -4.8148e-6  1.1799e-6

[:, :, 100] =
 -7.86683e-6  1.92763e-6
```
 as opposed to
 ```
 1×2×100 Array{Float64,3}:
[:, :, 1] =
 1.01316  -0.0107767

[:, :, 2] =
 -0.00614669  0.0014868

[:, :, 3] =
 0.00276786  -0.000444076

...

[:, :, 98] =
 0.00172214  -0.000271382

[:, :, 99] =
 0.00129058  -0.000175618

[:, :, 100] =
 -0.00617988  0.00166548
 ```

The commented one took 2.70 sec and the other took 15.38 sec.

# October 20, 2020

4:56 PM - Just now getting to research.I ran the code on a 2-D LSDE and it looks good. The autocovariances look pretty close. Tomorrow I would like to do model reduction.


# October 21, 2020

1:24 PM - I am happy to report that the current WF code (with the old xspectral estimator, the "direct estimator") *seems* to work well in model reduction for the liner SDE. Here is what I did (in a nut shell):
1. Set up the full model
   - B = -[-0.5 1; 0 -0.2]*[-0.5 1; 0 -0.2]'/1.5`
   - then called `modgen_LSDE` from `Examples/LinearSDE/modgen_LSDE.jl` as follows:
   ```julia
   # Model run Parameters
    t_start = 0
    t_stop  = 1e4
    h       = 1e-2

    A       = B
    σ       = [1 0; 0 1]
    Xo      = [1; 1]
    t_disc  = 100
    gap     = 1

    # Get full model run
    Random.seed!(2016)
    Y = modgen_LSDE(t_start,t_stop,h;
        A, σ, Xo, t_disc, gap)
   ```
    - This gave me the full model.

2. Then I reduced the model to only the first and then second variable variables.

```julia
   @time h_wf_1, pred = get_wf(Y[1:1,:], Psi;
    	M_out, n, p, par, ty, nfft, rl, Preds, PI, rtol);
   noise_dist = MvNormal(h*I + zeros(1,1))
   Y_rm_1 = redmodrun(real(Y[1:1,:]), h_wf_1, Psi; noise_dist)


   @time h_wf_2, pred = get_wf(Y[2:2,:], Psi;
     	M_out, n, p, par, ty, nfft, rl, Preds, PI, rtol);
   noise_dist = MvNormal(h*I + zeros(1,1))
   Y_rm_1 = redmodrun(real(Y[2:2,:]), h_wf_2, Psi; noise_dist)
```

I then ran all the usual tests
1. plot time series
2. plot distribution
3. plot autocovariance
4. plot spectral densities

The spectral densities were pretty much right on top of each other, in both cases. The autocovariances in the first variables diverged around 80 lags. In the second variable they were parallel though numerically off by about 20% the actual starting near 100 the reduced model around 80. These stayed roughly parallel through out. The distributions crudely approximated each other and the timeseries bore no red flags.

## UQ Meeting Notes

Dr. Lin and I ran through what we talked about last week about the different spectral density estimators with Will. I demonstrated a little and we looked at the performance of the WF in the double well Langevin case. We discussed the oddness of the estimators them selves being so close and yet the resulting Wiener filters so different. It was noted that the wiener filter with a (potentially) "better" spectral density approximation (i.e. `par` = 1000, `Nex` = 10^6) decayed slower (or experienced greater noise) than that of one with a "sloppier" spectral estimator (`par` = 55, `Nex` = 1024). The fact is that the sloppier WF decayed to ≈1e-8 at the 100th term, where the tighter WF decayed to only ≈1e-4, around the same order of magnitude of the third coefficient. (while writing this I wonder if that is attributed to the noise, I wonder if that will scale with the number of samples and how it may scale). It was suggested that I do these two things.
1. Update the function `get_wf` from `tools\Model_reduction_Dev.jl` to be able to switch which spectral estimator it is using. This will allow me a little more ease in investigating what is making the difference.
2. Investigate the singular values of the matrices used in the division of the S_yx by S_x^+ and the later division by S_y^-. It was speculated that small singular values in the denominator may exacerbate the difference in the cross-spectral estimations.

5:02 PM - I will try and finish the first before dinner and test it a little.

6:08 PM - Done for the day. I added the DWOL processes to *Cross Spectral Density Estimator notebook* that line up with the *Tester DWOL Notebook* in the `Examples/Tester`.  


# October 22, 2020

2:27 PM - Today I finished the investigative code for the xspectal estimators with the DWOL data. this notebook is labeled *Cross Spectral Density Estimator notebook* and is in the folder `Tools/ToolBoxTest`. I found that  

|`L`|`Nex`|`n`|`p`|'ty'|xspect  |WFDM|WFSP|
|---|---  |---|---|--- |---     |--- |--- |
|100|2^16 |3  |200|bin |close   |yes |no  |
|500|2^16 |3  |200|bin |v. close|yes |no  |

### Meeting with Dr. Lin

The last line in the table is very unusual. The estimates should not be so close with all other things being equal and their resulting WF so different. So, I will look at each step after the cross spectral estimate is used.

4:00 PM - I added some outputs to the `vecter_wiener_filter_fft` code:
```julia
h_num_fft = [h_num_raw[:,:,1:M],
             matlog1,
             matlog2,
             z_spect_pred_minus_num_fft,
             z_spect_pred_plus_num_fft,
             S_sigpred_overS_plus_fft_num,
             S_sigpred_overS_plus_plus_num_fft,
             H_num] ###
```


# October 23, 2020

8:18 PM - Today I did a lot of exploration. A bit of time was spent in understand the current landscape of the problem. For instance when I ran the WF with the direct estimator with `L = 55` the xspectral estimator of the xspect of the signal and predictors.
The plot showed a fairly course estimate near θ = 0. This was opposed to the estimte (on the same data) with `L = 10000` (`Nex = 2^16`).
This estimate was more smooth and more closely agreed the SP estimate (`p=300`, `n=3`, `ty = "bin"`).
However there was very little difference in the WF that resulted using the same data and these parameters. I conclude that methods enjoys an insensitivity to the xspectral estimation.
Which is a little enigmatic since it preforms very poorly under the smoothed periodogram method. The issue though I thing is that the larger the number of points used in the xspectral estimation the greater the error.
This is the error scales at a slightly higher power with the number `nfft`, at some point getting worse as `nfft` increases beyond some point. This is something I would like to investigate more later (future work).  

The Direct method WF was observed to be highly reproducible while the smoothed periodogram WF was more fickle. There were times when I put in the same parameters as before and use the same data but got different wiener filters. I suspect my julia installation has a problem.

So, I when I ran the SPWF on thelio it got closer even within the margins of error. So, now I did a few tests to verify that this is working. So, I made an uneven double well overdamped Langevin simulation and the DMWF passed that test. This test cane be found in the following notebook:

`Examples\Testers\Tester DWOL.ipynb`

The spectral estimator comparisons can be found in the notebook:

`Tools\ToolBoxTest\Cross Spectral Density Estimator Comparison.ipynb`

The model worked on the overdamped uneven Langevin. It got the analytically computed values.

# Wednesday, October 28, 2020

4:13 PM - Now I am going to investigate the variations I get in the xspectral estimators.


# Friday, October 30, 2020

10:41 AM - Today I will rewrite the smoothed periodogram estimator. I will make it so that the grid on which the xspectral density is evaluated is user defined. This will be done by segmenting the data into subseries of length defined by the user, computing the periodogram of each subseries and then averaging all the segment periodograms.

4:48 PM - Here is the code (so far)

```julia
function z_crossspect_scalar_ASP(
    sig,
    pred;
    nfft = 2^10, # The length of each subseries
    n = 3,
    p = 10,
    ty = "bin",
    L = nfft,
    win = "Par"
    )

    # Check length of series
    l_sig = length(sig)
    l_pred = length(pred)
    l_sig == l_pred || println("sizes must be the same, taking min and truncating")
    l = min(l_sig,l_pred)

    # The total nuber of subseries
    R = floor(Int,l/nfft)
    # The windowing function
    lam = win == "none" ? ones(nfft) : _window(nfft-1; win, two_sided = false)
    # Computation of the average periodogram
    aperi = complex(zeros(nfft))
    for r = 1:R
        fftsig = fft(lam .* sig[(r-1)*nfft+1:r*nfft])
        fftpred = conj(fft(lam .* pred[(r-1)*nfft+1:r*nfft]))
        aperi .+= fftsig .* fftpred
    end
    aperi ./= nfft*R

    # Smoothing it too.
    if ty != "none"
        aperi_pad = [aperi[end - p*n + 1 : end]; aperi; aperi[1:p*n]]
        μ = _smoother(n,p; ty)
        aperi = conv(μ,aperi_pad)[2n*p+1:2n*p+nfft]
    end
    aperi
end

```

  It should be noted that at this point I have included functionality for booth smoothing (by means of convolution after `fft` and windowing by means of multiplication in the time domain) I plan on removing one of these once I figure out which one to remove.

5:28 PM - Smoothing is working well, the windowing is not correct since I am off by a factor, I need to divide something to do with the windowing function. I think I need to divide by the sum of it's `dft`.

So far my tests have consisted of checking the estimated spectral density of a AR(2) process against it's now spectral density. The ASP is showing promise as it is much cleaner (smother and estimate)  



# Tuesday, November 3, 2020

12:59 PM - Taught and grade all morning. Now I will continue testing the averaged smoothed periodogram (ASP) cross spectral density estimate. Here's the agenda:

1. ARMA(p,q) processes. These are good because they have analytically computable spectral density. Cross spectral densities will need to be tested.
2. Double-welled overdamped Langevin processes. These will be compared to the two other forms of estimation I have available (he direct method and the full, smoothed periodogram method).
3. KSE and Lorenz '63. These will again just be compared to the other estimators. I will also compare the with what I may be able to find in the literature.

All of the tests will be conducted in the file `xspect_SP_new_dev.jl` in `Tools/ToolBoxTest`.

### Experiment Notes

#### ARMA

1:20 PM - Generated the following timeseries

```julia
steps = 10^6
l = [1, -.9, .5]

X = ARMA_gen(l; steps)
```

| `nn` | `pp` | `tty` | SP    | `L`  | `Nex` | `win` | DM        | `nfft` | `n`  | `p`  | `ty`   | `win`  | ASP                                   |
| :--: | ---- | ----- | ----- | ---- | ----- | ----- | :-------- | ------ | ---- | ---- | ------ | ------ | ------------------------------------- |
|  2   | 500  | "bin" | fuzzy | 100  | 2^11  | "Par" | Very good | 10^15  | 2    | 50   | "none" | "Par"  | Very poor off by a factor.            |
|  2   | 500  | "bin" | fuzzy | 100  | 2^11  | "Par" | Very good | 10^15  | 2    | 50   | "bin"  | "none" | very good little fuzzy better than SP |

I need to fix the windowing I think I need to divide by the sum of the `dft` of the window. I fixed the scaling of the windowing procedure. It turned out I need to divide the mean square to the window function. The estimate though is still **very** noisy. If fact when I did both windowing and smoothing the it was not as good (as in it varied further from the theoretic xspectral density) than when I just smoothed it. The support of the window function is `nfft`, that is, it is nonnegative on `-nfft:nfft`. I can imagine why it would be a good idea take that any smaller since it would be eliminating data. Perhaps if there was sufficient overlap in the intervals this might be OK and yield some benefit. I don't plan on exploring that now though, and for the moment would just as soon remove the windowing functionality.

For now I will leave it in and continued testing different ARMA



# Wednesday, November 4, 2020

12:48 PM - Yesterday, I tested a lot of ARMA processes of various types. The performance of the estimators was judge qualitatively.

6:07 PM - I made it so `z_crossspect_fft`calls `z_crossspect_ASP` and ran the Wiener filter



# Thursday, November 5, 2020

12:28 PM - I am going over the plots of the various elements used in the construction of the WF. I think this is a matter of the direct method varying with in bounds in an acceptable dimension which bounds the periodogram estimate violates. For instance, the WF code has very little sensitivity, it seems, to how the xspectral is approximated at 0.



# Monday, November 9, 2020

11:56 AM - Since I finished teaching at around 9:30 AM I have been investigating differences between my implementation and the `powerspect`estimator Dr. Lin supplied over the weekend. It's performance is comparable to the ones already being used. When I feed it in to the WFing program `get_WF` it preforms similarly to the my own periodogram estimates.  



# Thursday, November 12, 2020

12:28 PM - Finished writing an exam and grading and now am ready to research. Yesterday, I met with Dr. Lin and he drew my attention to the fact that most the work I have done the past several days was using a time series that was rather lean in it's information about the underlying process. What I need to do is allow for the numerical solution of the SDE to run over a longer time. The way this can be done efficiently with out losing any important information is to skip entries in the numerical run. So, yesterday I updated the function `DataGen_DWOL` in the file `Exmples/Nonlinear Langevin/DataGen.jl` to allow for a gap in the storing of the solution entries. The crux of it is demonstrated below: (particularly the second `if` statement)

```julia
signal = zeros(d,ceil(Int,steps/gap))
    tmp = sig_init
    if scheme == "FE"
        for n = 1 : steps_tot-1
            tmp = tmp + h*V_prime(tmp) +
                        sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    ...
```

I also, change some other things about the function to suite what I hope is a refined taste in programing. I made the time step `h` an input and took out the inputs `t_start` and `t-stop`. With this intact I will now rerun all the analysis before on this richer timeseries.

1:34 PM - I installed Anaconda so I could install jupytext to reconstitute the Jupyter notebooks and test my data generator.



# Monday, November 16, 2020

11:41 AM - In Math 122B they took a test today. I just finished updating their calendar.

The problem at hand for the past few weeks has been the question: **Why does the periodogram estimate WF give a different result form the direct method estimate WF?**

Looking closely at the estimators for the cross spectral density, under analogous parameters in resolution (`nfft = Nex`) the middle of the graph (`ω ∈ [.5,6]`) the values of the cross spectral density are very small and very close. So, we concluded that the difference must be at the approximations near zero. Last week Dr. Lin pointed out that the timeseries I had been using was in fact, rather lean (short, or brief the problem here was that given the current value of `σ = 0.3` the system did not run long enough in time for the sample to process adequate in formation about it's dynamics.) So we tried it again with longer timeseries and found that (skipping every 99 samples `gap = 100`) the WF were pretty much Identical. Can I find the old results? Can I reproduce them. I would like to be able to look back at these experiments quickly. I think I need to revise what I put into my research Journal here.

12:02 PM - I am preparing to run an experiment. I am thinking about writing a module that creates a nice output after an experiment.

#### Experiment Nov 16, 2020 1
* Run on my laptop.
* timeseries data
  ```julia
  #SDE parameters
  sigma    = [.3]
  V_prime  = x -> -x.*(x.^2 .- 1)
  sig_init = [1.5]
  # Numerical estimate parameters
  scheme   = "FE"
  steps    = 10^8  # Number of time steps (not including those discarded)
  h        = .1
  discard  = steps # Number of time steps discarded
  gap      = 100     # 1 + the number of time steps between observations
  seed     = 2015

  X = @time DataGen_DWOL(;
    #SDE parameters
    sigma, V_prime, sig_init,
    # Numerical estimate parameters
    scheme, steps, h, discard, gap)
  ```
  Data was saved to "Examples/Nonlinear Langevin/data/data_11_16_2020_1.jld"
  The `h = .1` here was a mistake but I will let it play out.

  ```
  julia> X
  1×1000000 Array{Float64,2}:
 -0.656241  -0.857331  -0.879254  -1.07485  -0.920843  -0.902613  -0.948586  0.90285  1.059  0.772112  0.944862  0.813503  1.07287  0.92711  0.805445  …  0.851886  0.897801  0.585505  0.979232  0.888067  0.816635  1.11003  1.11867  0.951974  1.08343  0.963044  1.00432  1.01226  0.887688  1.05448  1.08554
  ```
 ```
* We compute the WF and compare. The code is all in the file "Examples\Nonlinear Langevin\NonLinear Langevin Double Well.jl" Here are all the Wiener Filtering parameters:
  ```julia
  Psi(x) = [x; x.^3]
  M_out = 50
  ty = "bin"
  ###       xspect_est , par    , nfft    , n    , p
  Parms = [["DM"       , 5000  , 2^17    , 2    , 5],
           ["SP"       , 5000  , 2^17    , 2    , 5]]
 ```
* Results:
  ```
  First componete of the WF by DM: Complex{Float64}[1.3297488674434907 + 2.4073936672208036e-15im, -0.4545637642213554 - 8.671031880999923e-16im]
  First componete of the WF by SP: Complex[1.325886513646724 - 1.2871558917650188e-15im, -0.45249396088237676 + 4.436097622402122e-16im]
  ```

  ```
  julia> real(h_wf_dm)
  1×2×50 Array{Float64,3}:
  [:, :, 1] =
   1.32975  -0.454564

  [:, :, 2] =
   0.145294  -0.049749

  [:, :, 3] =
   0.017089  -0.00610873

  ...

  [:, :, 48] =
   0.00235385  -0.00223777

  [:, :, 49] =
   0.00255304  -0.00291539

  [:, :, 50] =
   0.00154927  -0.00160628

   julia> real(h_wf_sp)
   1×2×50 Array{Float64,3}:
   [:, :, 1] =
    1.32589  -0.452494

   [:, :, 2] =
    0.143381  -0.0484162

   [:, :, 3] =
    0.0170552  -0.00531202

   ...

   [:, :, 48] =
    0.00126733  -0.00129325

   [:, :, 49] =
    0.0033644  -0.00290029

   [:, :, 50] =
    0.00134207  -0.00142937
  ```

#### Experiment Nov 16, 2020 2
Same set up as above the differences being **only**
* Run on my desk top.
* timeseries data
  ```julia
  h        = .01
  ```
* Data and result was saved to "Examples/Nonlinear Langevin/data/data_11_16_2020_2.jld"
  ```
  julia> X
  1×1000000 Array{Float64,2}:
  -0.862223  -0.941636  -0.988963  -0.857691  -1.09023  …  -0.75433  -1.0951  -1.03065  -1.01327  -0.713247  -0.966905
  ```

* Results:
  ```
  First componete of the WF by DM: Complex{Float64}[1.068928167997716 - 2.983446612974894e-15im, -0.30447182295937214 + 1.1961030165977653e-15im]
  First componete of the WF by SP: Complex[1.071075375767013 + 2.7847331907457456e-16im, -0.30521795675917857 - 2.5149744546223904e-17im]
  ```

  ```
  julia> real(h_wf_packs[1])
  1×2×50 Array{Float64,3}:
  [:, :, 1] =
   1.06893  -0.304472

  [:, :, 2] =
   0.193165  -0.0723479

  [:, :, 3] =
   0.0969474  -0.0359304

  ...

  [:, :, 48] =
   -0.000477763  0.00062005

  [:, :, 49] =
   -0.00230393  0.00106569

  [:, :, 50] =
   0.00169566  -0.000434639

  julia> real(h_wf_packs[8])
  1×2×50 Array{Float64,3}:
  [:, :, 1] =
   1.07108  -0.305218

  [:, :, 2] =
   0.194427  -0.0726695

  [:, :, 3] =
   0.0984909  -0.0363122

  ...

  [:, :, 48] =
   -0.000198413  0.000605992

  [:, :, 49] =
   -0.0033782  0.0014048

  [:, :, 50] =
   0.000420223  4.8856e-6
  ```


  #### Experiment Nov 16, 2020 3 (thelio job 144)
  Same set up as above the differences being **only**
  * Run on my thelio.
  * timeseries data
    ```julia
    steps    = 10^7
    h        = .1
    gap      = 1
    ```
  * Data and result was saved to "~/data/data_11_16_2020_3.jld"
    ```
    julia> X = load("data_11_16_2020_3.jld","X")
    1×10000000 Array{Float64,2}:
    0.963992  1.08766  1.03263  1.15454  1.03817  …  1.17855  1.06523  1.08276  1.10181  0.0
    ```

  * Results:
    ```
    First componete of the WF by DM: Complex{Float64}[1.0959050884568995 + 1.3295451807051166e-14im, -0.09853703933071341 - 4.945113857281478e-15im]
    First componete of the WF by SP: Complex[1.2206872066904122 - 3.948398830887976e-14im, -0.14026744444907263 + 1.1469977522780055e-14im]
    ```

    ```
    julia> real(h_wf_packs[1])
    1×2×50 Array{Float64,3}:
    [:, :, 1] =
     1.09591  -0.098537

    [:, :, 2] =
     0.00262928  -0.0010732

    [:, :, 3] =
     0.00146052  -0.000649314

    ...

    [:, :, 48] =
     -0.000233574  -6.58229e-5

    [:, :, 49] =
     0.00115693  -0.00024298

    [:, :, 50] =
     -0.000572883  0.000422376

     julia> real(h_wf_packs[8])
     1×2×50 Array{Float64,3}:
     [:, :, 1] =
      1.22069  -0.140267

     [:, :, 2] =
      0.0253809  -0.00989432

     [:, :, 3] =
      0.0197208  -0.00779881

     ...

     [:, :, 48] =
      -0.00121644  0.000298723

     [:, :, 49] =
      -0.000267469  0.000250607

     [:, :, 50] =
      -0.0017526  0.000833896
    ```

#### Experiment Nov 16, 2020 4 (thelio job 145)
Same set up as above 3, exactly save one thing `seed = 2016`.
* Data and result was saved to "~/data/data_11_16_2020_4.jld"
  ```
  julia> X = load("data_11_16_2020_4.jld","X")
  1×10000000 Array{Float64,2}:
  0.938291  0.980287  1.07918  1.07581  1.02021  …  0.330087  0.35004  0.36474  0.344513  0.0
  ```

* Results:
  ```
  First componete of the WF by DM: Complex{Float64}[1.096979669475226 + 1.404646908929906e-14im, -0.09883483887030214 - 4.739714180791289e-15im]
  First componete of the WF by SP: Complex[1.223171797163009 - 1.6490594744638936e-14im, -0.14090005527410765 + 4.523301785279795e-15im]
  ```

  ```
  julia> real(h_wf_packs[1])
  1×2×50 Array{Float64,3}:
  [:, :, 1] =
   1.09698  -0.0988348

  [:, :, 2] =
   0.00183088  -0.000800667

  [:, :, 3] =
   -0.000678346  2.57955e-5

  ...

  [:, :, 48] =
   0.000393262  8.22848e-5

  [:, :, 49] =
   -0.000153363  -2.6423e-5

  [:, :, 50] =
   -0.000700438  0.000167881

  julia> real(h_wf_packs[8])
  1×2×50 Array{Float64,3}:
  [:, :, 1] =
   1.22317  -0.1409

  [:, :, 2] =
   0.0246669  -0.0096447

  [:, :, 3] =
   0.0180439  -0.00724898

  ...

  [:, :, 48] =
   -0.000890076  0.00056583

  [:, :, 49] =
   -0.00142923  0.00045047

  [:, :, 50] =
   -0.0017453  0.00058458
  ```


# Tuesday, Nov 17, 2020

3:34 PM - Today I continue the experiments of yesterday. The goal will be to analysis the noise in the wiener filter. So, here is the plan I will repeat *Experiment Nov 16, 2020 3* 100 times with different data and save the wiener filters. Then I will find the mean and analyze the variance. I will repeat this with `steps = 10^6` and `steps = 10^8` (though maybe only running 50 or so in the ensemble).

The first thing is to write the code. I will do so in the directory `Examples/Nonlinear Langevin` under the file name `DWOL_ens_tests.jl`. These are going to be written to be run on thelio.

5:09 PM - I ran the first experiment for the day.

#### Experiment Nov 17, 2020 1 (thelio job 147)
This experiment is, as is mentioned above, an experiment in analyzing the error in the WF. Particularly, I wish to see how the error scales with the number of effective samples.
The experiment is all contained within the file `Examples/Nonlinear Langevin/DWOL_ens_tests.jl`
This is the first experiment in a series. What we do is generator an ensemble of `Nens = 100` independent runs of the timeseries with the parameters given below. For each run we compute the wiener filter using the (1) direct method (DM)for cross spectral density approximation and (2) the smoothed periodogram (SP). We then save the `2*Nens` Wiener filters for statistical analysis later.

* Data and result was saved to "~/data/data_11_16_2020_4.jld" The parameters were
  ```julia
  #SDE parameters
  sigma    = [.3]
  V_prime  = x -> -x.*(x.^2 .- 1)
  sig_init = [1.5]
  # Numerical estimate parameters
  scheme   = "FE" # Forward euler
  steps    = 10^7  # Number of time steps (not including those discarded)
  h        = .1
  discard  = steps # Number of time steps discarded
  gap      = 1     # 1 + the number of time steps between observations

  Psi(x) = [x; x.^3]

  # Model reduction Parameters
  M_out = 20
  ty = "bin"
  ### Varing parameters
  ###       xspect_est , par    , nfft    , n    , p
  #
  Parms = [["DM"       , 5000  , 2^17    , 2    , 5, ty, M_out],
           ["SP"       , 5000  , 2^17    , 2    , 5, ty, M_out]]

  Nens = 100
  ```

* Results: The data was saved in the file "~/data/data_11_17_2020_1.jld". In the file are  the above parameters and the output variable `h_wfs_ens` which is a 5-dimensional array. The first two dimension are the row in `Parms` that was used and the ensemble index, respectively the last three dimensions are for the wiener filter.

| `xspect_est`| statistic | h_wf[1,1,1] | h_wf[1,1,1] |
|---|---|---|---|
|"DM" | mean | 1.095424882549188 | -0.09837284642681761|
|"SP" | mean | 1.2212270201731081| -0.140304542799513 |
|"DM" | var | 1.1940063094075038e-6 | 1.3031776243201438e-7|
|"SP" | var | 1.2444470995957716e-5| 1.5089635615229606e-6 |

#### Experiment Nov 17, 2020 2 (thelio job 148)

This experiment was the same as the *Experiment Nov 17, 2020 1* with one exception: `steps = 10^6` the data was saved, as described above, in "~/data/data_11_17_2020_2.jld". Here are some of the results.

| `xspect_est`| statistic | h_wf[1,1,1] | h_wf[1,1,1] |
|---|---|---|---|
|"DM" | mean |1.0955532773712064 | -0.0984309622679665|
|"SP" | mean | 1.2175349889961529| -0.13915486525042936 |
|"DM" | var | 1.0964189610161555e-5 | 1.2682319503192065e-6|
|"SP" | var | 0.00014711838831749758| 1.846553738659983e-5 |



# Wednesday, Nov 18, 2020

12:18 PM - Just lost a lot of writing because Atom crashed on my as I was writing in this journal. I think there might be a problem running the markdown previewer at the same time. Anyway, I just finished grading the last midterm exam for Math 122B, and would like to do a bit more grading today, as I am a little behind. But first I want to get some experiment running. Namely, (1) on thelio I want to do an ensemble study of *Experiment Nov 16, 2020 2* with `Nens = 100`. So, that means `gap = 100`, `h = 0.01`, and `steps = 10^8`. I have that this works well for both WF estimators. Then I would like to repeat the experiment for `steps = 10^7` and `steps = 10^6` so that I can see what happens when I use too few samples. Hopefully the transition from working to not working will be instructive. I already have a lot of examples of these WF's not giving the same results. (2) I would like to look at the guts of the particular runs of the Wiener filters and see what is happening. So, I will run the *Experiment Nov 16, 2020 3* again with `info = true` in the `get_wf` call. This experiment I will run on my desktop.

**Question:** I would like to better understand the relationship between the continuous time filter for the continuous time solution of the DWOL problem and their discrete time approximations. Particularly I want to understand in what way these WF estimates may approach the theoretical continuous time filter. Is it as `steps` goes to infinity or as `nfft` or `Nex` goes to infinity or both?

#### Experiment Nov 18, 2020 1 (thelio job 150)

This experiment an ensemble study of *Experiment Nov 17, 2020 2*, in which the resulting WF matched pretty well. What I would like to do now is get an idea of the mean and variance of these filters. The experiment is running using `Langevin/DWOL_ens_tests.jl` on thelio and the data and results are saved in the file "~/data/data_11_18_2020_1.jld". Here are some of the results.

| `xspect_est`| statistic | h_wf[1,1,1] | h_wf[1,1,1] |
|---|---|---|---|
|"DM" | mean | 1.0724217934082474    | -0.30594141224404103 |
|"SP" | mean | 1.074365601988284     | -0.3066500631917417  |
|"DM" | var  | 3.966369819921137e-5  | 4.652132489933455e-6 |
|"SP" | var  | 5.3503771439725386e-5 | 7.24717987093219e-6  |

#### Experiment Nov 18, 2020 2 (thelio job 154)

This experiment is the same as *Experiment Nov 18, 2020 1*, but with `steps = 10^7` The experiment is running using `Langevin/DWOL_ens_tests.jl` on thelio and the data and results are saved in the file "~/data/data_11_18_2020_2.jld". Here are some of the results.

| `xspect_est`| statistic | h_wf[1,1,1] | h_wf[1,1,1] |
|---|---|---|---|
|"DM" | mean | 1.0685370878746394    | -0.30464623282477016|
|"SP" | mean | NaN                   | NaN                 |
|"DM" | var  | 0.0005199931437080256 | 6.218700778799889e-5|
|"SP" | var  | NaN                   | NaN                 |

#### Experiment Nov 18, 2020 3 (thelio job 155)

This experiment is the same as *Experiment Nov 18, 2020 1*, but with `steps = 10^6` The experiment is running using `Langevin/DWOL_ens_tests.jl` on thelio and the data and results are saved in the file "~/data/data_11_18_2020_3.jld". Here are some of the results.

| `xspect_est`| statistic | h_wf[1,1,1] | h_wf[1,1,1] |
|---|---|---|---|
|"DM" | mean | 1.0799339828734176   | -0.30799075003241844|
|"SP" | mean | NaN                  | NaN                 |
|"DM" | var  | 0.004598458570730727 | 0.000564803873321387|
|"SP" | var  | NaN                  | NaN                 |

It looks like the functions I used did not support times series of length less than `nfft`


#### Experiment Nov 18, 2020 4 (desktop)

This experiment is the same as *Experiment Nov 16, 2020 3*. That is it runs a DWOL model with σ = .3 (as always) and `steps = 10^7`, `h = 0.1`, and `gap = 1`. This time I am looking at the plots of the algorithm steps. `Nonlinear Langevin/NonLiner Langevin Double Well.jl` was used on my desktop and the data and results are saved in the file "Examples/Nonlinear Langevin/data/data_11_18_2020_4.jld". I did the array of plots that Dr. Lin requested last week namely:

|real | imaginary|
|---|---|
|S_yx DM and SP | S_yx DM and SP |
|S_yx/S_x^+ DM and SP | S_yx/S_x^+ DM and SP |
|{S_yx/S_x^+}_+ DM and SP |{S_yx/S_x^+}_+ DM and SP |
|{S_yx/S_x^+}_+/s_x^- DM and SP |{S_yx/S_x^+}_+/s_x^- DM and SP|

The plots were produced using the script `Nonlinear Langevin/Plots_for_DrLinNov11_2020.jl`.

In the plots it was apparent that the cross spectral density estimated using the periodogram was about twice as high as the direct method analogue. meaning the difference in power near 0 was material. That difference then get spread out by the causal-part projection operator.

#### Experiment Nov 18, 2020 5 (desktop)

This experiment is the same as *Experiment Nov 16, 2020 2*. That is it runs a DWOL model with σ = .3 (as always) and `steps = 10^8`, `h = 0.01`, and `gap = 100`. This time I am looking at the plots of the algorithm steps.  Again I used the scripts `Nonlinear Langevin/NonLiner Langevin Double Well.jl` and `Plots_for_DrLinNov11_2020.jl` on my desktop. The data and results are saved in the file "Examples/Nonlinear Langevin/data/data_11_18_2020_5.jld". The result here was that the estimator near 0 were in much closer agreement. and so there closeness survived the casual-part operator.

## UQ meeting notes

Today during UQ I reported the problems and we looked at the above plots and came to the above conclusions (I had written those after the meeting). It was also speculated that the reason for the closeness in the estimators with data of `gap = 100` over those with data of `gap = 1` (even though the time span of the run was the same) was the difference in smoothing. The `gap = 100` data (`h = 0.01`) covered a time span of 10^6 sec, but with only 10^6 points. The  `gap = 1` data (`h = 0.1`) covered a covered the same time span, but with 10^7 points.



# Thursday, Nov 19, 2020
2:50 PM - First, I want to see the effect of increasing the `gap` parameter on a parameter. So, I will have thelio (I think) generate a run `X` of `steps = 10^8` and `h=0.01`, and then run the wiener filters (one for each approximation method) on the series `X`, `X[:,1:10:end]`, and `X[:,1:100:end]` and compare the two wiener filters. I suspect the first one to be pretty different from each other but the last one should match up rather well.

#### Experiment Nov 19, 2020 1 (varying the gap for the same series)
I never got to this experiment as I think it's time to start working on KSE.
