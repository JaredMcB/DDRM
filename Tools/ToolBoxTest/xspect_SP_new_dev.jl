"""
Here I modify the function function `z_crossspect_scalar` in `AnalysisToolBox.jl`
to allow for the xspectal estimate to have a user specified number of number
of grid points in the estimation.

The function that is used in the WF code `vector_wiener_filter_fft`in 'Model_redction_Dev.jl'
is `z_crossspect_fft` this is just a matrix-valued organizer that calls `z_crossspect_scalar`
for it's individual xsprectrum.
"""

using PyPlot
using Random
using JLD

include("../../Examples/Nonlinear Langevin/DataGen.jl")
include("../../Tools/Wiener Filtering/Scalar Wiener Filter/ARMA_Generator_DSP.jl")
include("../../Tools/AnalysisToolbox.jl")

###


### generic spectral plotter
function spect_plot(S;
        plotter = plot,
        label = "unlabled")
    plotter(2pi*(0:length(S)-1)/length(S),S,
            label = label)
    legend()
end

##########################################################
### ARMA Testing #########################################
##########################################################

### Run 1 ################################################
### Get Data
steps = 10^6
l = [1, -.9, .5]

X = ARMA_gen(l; steps)

### Plot analytic Solution
S_ana_poly = Polynomial(l)
S_ana_fun(z) = 1/(S_ana_poly(z^(-1))*S_ana_poly(z')')
Theta = 2pi*(0:1000-1)/1000
plot(Theta,S_ana_fun.(exp.(im*Theta)))

nn = 2
pp = 500
tty = "bin"
S_X_sp = z_crossspect_scalar(X,X;n = nn ,p = pp ,ty = tty)
spect_plot(S_X_sp)

X_vec = reshape(X,1,:)
L = 100
Nex = 2^11
win = "Par"
S_X_dm = z_crossspect_fft_old(X_vec,X_vec;L,Nex,win)[:]
spect_plot(S_X_dm)

nfft = 2^15
n = 2
p = 50
ty = "bin"
win = "none"
S_X_asp = z_crossspect_scalar_ASP(X,X; nfft, n, p,ty, win)
spect_plot(S_X_asp)



### Run 2 ################################################
### Get Data
steps = 10^6
l = [1, -.9, .5]
w = randn(3)

X = ARMA_gen(l,w; steps)

### Plot analytic Solution
S_ana_poly_d = Polynomial(l)
S_ana_poly_n = Polynomial(w)
S_ana_fun(z) = (S_ana_poly_n(z^(-1))*S_ana_poly_n(z')')/
        (S_ana_poly_d(z^(-1))*S_ana_poly_d(z')')
Theta = 2pi*(0:1000-1)/1000
plot(Theta,S_ana_fun.(exp.(im*Theta)))

nn = 2
pp = 500
tty = "bin"
S_X_sp = z_crossspect_scalar(X,X;n = nn ,p = pp ,ty = tty)
spect_plot(S_X_sp)

X_vec = reshape(X,1,:)
L = 100
Nex = 2^11
win = "Par"
S_X_dm = z_crossspect_fft_old(X_vec,X_vec;L,Nex,win)[:]
spect_plot(S_X_dm)

nfft = 2^15
n = 2
p = 50
ty = "bin"
win = "none"
S_X_asp = z_crossspect_scalar_ASP(X,X; nfft, n, p,ty, win)
spect_plot(S_X_asp)



### Run 3 ################################################
### Get Data
steps = 10^6
discard = 10^4

p = rand(1:9)
q = rand(1:9)

Zeros_sig = 1 .- rand(q)*2
Poles_sig = 1 .- rand(p)*2

X, P, Q = ARMA_gen(;
    Zeros = Zeros_sig,
    Poles = Poles_sig,
    steps,
    discard,
    out_poly = true)

### Plot analytic Solution
S_ana_fun(z) = (Q(z^(-1))*Q(z')')/
        (P(z^(-1))*P(z')')
Theta = 2pi*(0:1000-1)/1000
plot(Theta,
    S_ana_fun.(exp.(im*Theta)),
    label = "analytic")

nn = 2
pp = 500
tty = "bin"
S_X_sp = z_crossspect_scalar(X,X;n = nn ,p = pp ,ty = tty)
spect_plot(S_X_sp,label = "S_X_sp")

X_vec = reshape(X,1,:)
L = 1000
Nex = 2^11
win = "Par"
S_X_dm = z_crossspect_fft_old(X_vec,X_vec;L,Nex,win)[:]
spect_plot(S_X_dm,label = "S_X_dm")

nfft = 2^15
n = 2
p = 100
ty = "ave"
win = "none"
S_X_asp = z_crossspect_scalar_ASP(X,X; nfft, n, p,ty, win)
spect_plot(S_X_asp,label = "S_X_asp")


### Run 4 ################################################
### Get Data
steps = 10^6
discard = 10^4

p1 = rand(1:9)
p2 = rand(1:4)
q1 = rand(1:9)
p2 = rand(1:4)

Zeros_sig_real = 1 .- rand(q1)*2
com_zeros = rand(p2) .* exp.(2pi*im*rand(p2))
Zeros_sig = [Zeros_sig_real; com_zeros; conj(com_zeros)]
Poles_sig_real = 1 .- rand(p1)*2
com_poles = rand(p2) .* exp.(2pi*im*rand(p2))
Poles_sig = [Poles_sig_real; com_poles; conj(com_poles)]


X, P, Q = ARMA_gen(;
    Zeros = Zeros_sig,
    Poles = Poles_sig,
    steps,
    discard,
    out_poly = true)

### Plot analytic Solution
S_ana_fun(z) = (Q(z^(-1))*Q(z')')/
        (P(z^(-1))*P(z')')
Theta = 2pi*(0:1000-1)/1000
plot(Theta,
    S_ana_fun.(exp.(im*Theta)),
    label = "analytic")

nn = 2
pp = 500
tty = "bin"
S_X_sp = z_crossspect_scalar(X,X;n = nn ,p = pp ,ty = tty)
spect_plot(S_X_sp,label = "S_X_sp")

X_vec = reshape(X,1,:)
L = 500
Nex = 2^10
win = "Par"
S_X_dm = z_crossspect_fft_old(X_vec,X_vec;L,Nex,win)[:]
spect_plot(S_X_dm,label = "S_X_dm")

nfft = 2^16
n = 2
p = 50
ty = "ave"
win = "none"
S_X_asp = z_crossspect_scalar_ASP(X,X; nfft, n, p,ty, win)
spect_plot(S_X_asp,label = "S_X_asp")







## Testing
##########################################################
### DWOL Testing #########################################
##########################################################

# Model run Parameters
steps = 10^7 + 1
scheme = "FE"
t_start = 0
t_stop = 10^5
discard = 100000
sig_init = [1.5]
sigma = [.5]
V_prime = x -> -x.*(x.^2 .- 1)
SM1 = false
Obs_noise = false
d = 1
# e = randn(d,steps + discard)

dt = (t_stop - t_start)/(steps - 1)


# Get full model run
Random.seed!(2014)
X = DataGen_DWOL(
    steps;
    scheme, t_start, t_stop, discard,
    sig_init , sigma, V_prime,
    SM1, Obs_noise, d
    )
### This is how we set up the signals and
### predictors They have to me off set so
### that pred is one index behind signal
### i.e. pre(n) = psi(sig(n-1))
X_sig = X[:,2:end];

Psi(x) = [x; x.^3]
X_pred = get_pred(X,Psi) # Notice it is just
                         # X get_pred assigns
                         # psi straight across


 nn = 2
 pp = 500
 tty = "bin"
 S_X_sp = z_crossspect_scalar(X,X;n = nn ,p = pp ,ty = tty)
 spect_plot(S_X_sp,label = "S_X_sp")

 X_vec = reshape(X,1,:)
 L = 500
 Nex = 2^10
 win = "Par"
 S_X_dm = z_crossspect_fft_old(X_vec,X_vec;L,Nex,win)[:]
 spect_plot(S_X_dm,label = "S_X_dm")

 nfft = 2^16
 n = 2
 p = 50
 ty = "ave"
 win = "none"
 S_X_asp = z_crossspect_scalar_ASP(X,X; nfft, n, p,ty, win)
 spect_plot(S_X_asp,label = "S_X_asp")

### General ARMA #############################################
X = ARMA_gen(  l = [1, -5/4, 3/8],
                    w = [1];
                    r::Float64 = 1.0,
                    steps::Int64 = 10^4,
                    Zeros = [],
                    Poles = [],
                    e = [],
                    discard::Int64 = 10^3)

X_sig = X[:,2:end];

Psi(x) = [x; x.^3]
X_pred = get_pred(X,Psi) # Notice it is just
                         # X get_pred assigns
                         # psi straight across





### KSE #######################################################
gen = gen
T = 3000 # Length (in seconds) of time of run
T_disc = 1000 # Length (in seconds) of time discarded
P = 32π  # Period
N = 128  # Number of fourier modes used
h = 1e-3 # Timestep
g = x -> cos(π*x/16)*(1 + sin.(π*x/16))
obs_gap = 100

uu, vv, tt =  my_KSE_solver(T,
       T_disc  = T_disc,
       P = P,
       N = N,
       h = h,
       g = g,
       n_gap = obs_gap)








### Write the test first
