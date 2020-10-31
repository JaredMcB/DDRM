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

### Get Data
steps = 10^6
l = [1, -.9, .5]

X = ARMA_gen([1, -.9, .5]; steps)

### Plot analytic Solution
S_ana_poly = Polynomial(l)
S_ana_fun(z) = 1/(S_ana_poly(z^(-1))*S_ana_poly(z')')
Theta = 2pi*(0:1000-1)/1000
plot(Theta,S_ana_fun.(exp.(im*Theta)))

### Plot other solutions
spect_plot(S) = plot(2pi*(0:length(S)-1)/length(S),S)

n = 2
p = 50
ty = "bin"
S_X_sp = z_crossspect_scalar(X,X;n,p,ty)
spect_plot(S_X_sp)

X_vec = reshape(X,1,:)

L = 100
Nex = 2^11
win = "Par"
S_X_dm = z_crossspect_fft_old(X_vec,X_vec;L,Nex,win)[:]
spect_plot(S_X_dm)

### Work on new averaged periodogram
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

nfft = 2^15
n = 2
p = 50
ty = "none"
win = "Par"
S_X_asp = z_crossspect_scalar_ASP(X,X; nfft, n ,p ,ty ,win)

spect_plot(S_X_asp)

_window(50;win,two_sided = false)






sig = X
pred = X
n = 3
p = 10
ty = "bin"
L = nfft
win = "Par",
nfft = 2^10

# Check length of series
l_sig = length(sig)
l_pred = length(pred)
l_sig == l_pred || println("sizes must be the same, taking min and truncating")
l = min(l_sig,l_pred)

# The total nuber of subseries
R = floor(Int,l/nfft)
# The windowing function
lam = win == "none" ? ones(nfft) : _window(nfft; win, two_sided = false)
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


























plot(μ)
μ_hat = ifft(μ)
plot(μ_hat)

wind = _window(50,win = "gigi")
plot(wind)
plot(fft(wind))








function z_crossspect_scalar_ASP(sig,pred;
    nfft = 0,
    n    = 3,
    p    = 100,
    ty   = "ave",
    win  = "par")

    μ = _smoother(n,p;ty)

    l_sig = length(sig)
    l_pred = length(pred)
    l_sig == l_pred || println("sizes must be the same, taking min and truncating")
    l = min(l_sig,l_pred)

    nfft = nfft == 0 ?



    R = floor(Int,l/nfft)
    aperi = complex(zeros(nfft))
    for r = 1:R
        fftsig = fft(sig[(r-1)*nfft+1:r*nfft])
        fftpred = conj(fft(pred[(r-1)*nfft+1:r*nfft]))
        aperi .+= fftsig .* fftpred
    end

    aperi_pad = [aperi[end - p*n + 1 : end]; aperi; aperi[1:p*n]]

    μ = _smoother(n,p; ty)
    z_crsspect_smoothed = conv(μ,aperi_pad)[2n*p+1:2n*p+nfft]





    # nfft == l || println("adjusted size from $l to $nfft")
    sig_pad = l < nfft ? [sig[1:l]; zeros(nfft - l)] : sig[1:nfft]
    pred_pad = l < nfft ? [pred[1:l]; zeros(nfft - l)] : pred[1:nfft]

    fftsig = fft(sig_pad)
    fftpred = conj(fft(pred_pad))

    peri = fftsig .* fftpred / nfft
    peri_pad = [peri[end - p*n + 1 : end]; peri; peri[1:p*n]]
    z_crsspect_smoothed = conv(μ,peri_pad)[2n*p+1:2n*p+nfft]
end














## Testing

### DWOL ###########################################
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
