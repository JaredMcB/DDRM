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

4:47 PM - Today the goal has been to get a big batch of scalar linear SDE's running do I can be sure the code works. And I have just barely done that. I am starting with a little thing that ran comfortably on Dell, this job is running now (job 134).

###Experiment Table: Linear SDE
| A | h   | `t_stop` |  `gap` | `M_out` | `N_eff` |
|---|---  |---       |---     |---      |---      |
|-.5|0.01 | 1e5      | 1      |20       |         |
