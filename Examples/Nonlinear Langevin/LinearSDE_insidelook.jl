using PyPlot
using Random
using JLD


# Get software to generate model
include("../LinearSDE/modgen_LSDE.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")

# Model run Parameters
t_start = 0
t_stop  = 1e4
h       = 1e-2


A       = -[-0.5 1; 0 -0.2]*[-0.5 1; 0 -0.2]'/1.5 #reshape([-0.5],1,1)
σ       = I + 0*A
Xo      = [1; 1.5]
t_disc  = 1e3
gap     = 1

# Get full model run
Random.seed!(2014)
X = modgen_LSDE(t_start,t_stop,h;
    A, σ, Xo, t_disc, gap)
steps = floor(Int,(t_stop - t_start)/(h*gap)) +1

data = Dict("steps" => steps,
            "t_stop" => t_start,
            "sigma" => σ,
            "X" => X)
save("Examples/Nonlinear Langevin/data/data_10_23_2020Linear.jld",data)

data = load("Examples/Nonlinear Langevin/data/data_10_23_2020Linear.jld")
X = data["X"]

auto_times(X[1,:])

# Put in Psi functions
Psi(x) = x

# Model reduction Parameters
M_out = 50
n = 3
p = 500
par = 55
ty = "bin"
xspec_est = "DM"
rl = true
Preds = false
PI = false
rtol = 1e-6


### Varing parameters
###            xspect_est, par    , nfft    , n    , p
#
Parms = [["DM"       , 500  , 2^10    , 2    , 5],
         ["SP"       , 500  , 2^10    , 2    , 5],
         ["DM"       , 100    , 2^10    , 3    , 500],
         ["DM"       , 500    , 2^16    , 3    , 500],
         ["DM"       , 1000   , 2^16    , 3    , 500],
         ["DM"       , 5000   , 2^16    , 3    , 500],
         ["DM"       , 10000  , 2^16    , 3    , 500],
         ["DM"       , 5000   , 2^20    , 3    , 500],
         ["DM"       , 10000  , 2^20    , 3    , 500],
         ["SP"       , 10000  , 2^10    , 2    , 5],
         ["SP"       , 10000  , 2^16    , 2    , 5],
         ["SP"       , 10000  , 2^17    , 2    , 5],
         ["SP"       , 10000  , 2^17    , 3    , 10]]

nfft      = Parms[1][3]

P = 2#length(Parms)

h_wf_packs  = []
times = zeros(P)
for i = 1:P
    Out = @timed get_wf(X, Psi;
        M_out, ty, rl, Preds, PI, rtol, info = true,
        xspec_est = Parms[i][1],
        par       = Parms[i][2],
        nfft      = Parms[i][3],
        n         = Parms[i][4],
        p         = Parms[i][5]);

    append!(h_wf_packs,Out.value)
    times[i]      = Out.time
end

h_wf_dm = h_wf_packs[1]
h_wf_sp = h_wf_packs[8]

h_wf_dm[:,:,1]
h_wf_sp[:,:,1]

I +h*A

z_crossspect_sigpred_num_fft_dm = h_wf_packs[2]
z_crossspect_sigpred_num_fft_sp = h_wf_packs[9]

S_sigpred_overS_plus_fft_num_dm = h_wf_packs[5]
S_sigpred_overS_plus_fft_num_sp = h_wf_packs[12]

S_sigpred_overS_plus_plus_num_fft_dm = h_wf_packs[6]
S_sigpred_overS_plus_plus_num_fft_sp = h_wf_packs[13]

H_num_dm = h_wf_packs[7]
H_num_sp = h_wf_packs[14]

fig, axs = subplots(4,2,sharey = "row",sharex = true)

axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm1")
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[2,2,:]),label = "dm1")

axs[1].set_title("Direct Method")
axs[1].set_ylabel("S_XY/S_X^+ (real)")
axs[1].axis([0, 6.28,-5e-2,2e-1])
axs[1].grid("on")
axs[1].legend()

axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm1")
axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[2,2,:]),label = "dm1")
axs[2].set_ylabel("S_XY/S_X^+  (lmag)")
axs[2].grid("on")
axs[2].axis([0, 6.28,-2e-1,2e-1])

axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm1")
axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[2,2,:]),label = "dm1")
axs[3].set_ylabel("{S_XY/S_X^+}_+  (real)")
axs[3].grid("on")
axs[3].axis([0, 6.28,-1e-2,2e-1])

axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm1")
axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[2,2,:]),label = "dm1")
axs[4].set_ylabel("{S_XY/S_X^+}_+  (imag)")
axs[4].grid("on")
axs[4].axis([0, 6.28,-2e-1,2e-1])
axs[4].set_xlabel("frequencies")


axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp1")
axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[2,2,:]),label = "sp1")
axs[5].set_title("Periodogram")
axs[5].grid("on")
axs[5].legend()

axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp1")
axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[2,2,:]),label = "sp1")
axs[6].grid("on")

axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp1")
axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[2,2,:]),label = "sp1")
axs[7].grid("on")

axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp1")
axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[2,2,:]),label = "sp1")
axs[8].grid("on")
axs[8].set_xlabel("frequencies")
fig.suptitle("Causal projection of Estimated spectral densities over S_pred_plus")



























## Extimated Spectral Densities
z_crossspect_sigpred_num_fft_dm = h_wf_packs[2]
z_crossspect_sigpred_num_fft_sp = h_wf_packs[9]

semilogy(2pi*(0:nfft-1)/nfft,(real(z_crossspect_sigpred_num_fft_dm[1,1,:])),label = "dm1")
semilogy(2pi*(0:nfft-1)/nfft,(real(z_crossspect_sigpred_num_fft_sp[1,1,:])),label = "sp1")
legend()
axis([.5, 5.5,-1e-1,1e1])

nfft = 2^17
plot(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_sp[1,1,:]),label = "sp1")
legend()
axis([.5, 5.5,-1e-2,2.5e-2])

plot(2pi*(0:nfft-1)/nfft,imag(z_crossspect_sigpred_num_fft_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,imag(z_crossspect_sigpred_num_fft_dm[1,2,:]),label = "dm2")
plot(2pi*(0:nfft-1)/nfft,imag(z_crossspect_sigpred_num_fft_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,imag(z_crossspect_sigpred_num_fft_sp[1,2,:]),label = "sp2")
legend()
axis([.5, 5.5,-2.5e-2,2.5e-2])

fig, axs = subplots(1,2,sharey = true)
axs[1].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_dm[1,1,:]),label = "dm1")
axs[1].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_dm[1,2,:]),label = "dm2")
axs[1].set_title("Direct Method")
axs[1].set_xlabel("frequencies")
axs[1].set_ylabel("Spectral density")
axs[1].grid("on")
axs[1].legend()

axs[2].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_sp[1,1,:]),label = "sp1")
axs[2].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_sp[1,2,:]),label = "sp2")
axs[2].set_title("Periodogram")
axs[2].set_xlabel("frequencies")
axs[2].set_ylabel("Spectral density")
axs[2].grid("on")
axs[2].legend()
fig.suptitle("Cross spectral density estimates")

SVD1s = [svd(S_pred_plus[:,:,n]).S[1] for n = 1:nfft]
SVD2s = [svd(S_pred_plus[:,:,n]).S[2] for n = 1:nfft]



## Estimated spectral densities over S_pred_plus

S_sigpred_overS_plus_fft_num_dm = h_wf_packs[5]
S_sigpred_overS_plus_fft_num_sp = h_wf_packs[12]

semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_sp[1,1,:])),label = "sp1")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_sp[1,2,:])),label = "sp2")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_dm[1,1,:])),label = "dm1")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_dm[1,2,:])),label = "dm2")
legend()
axis([.5, 5.5,-1e-1,1e1])

plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
legend()
title("Estimated spectral densities over S_pred_plus")
axis([.25, 6,-2.5e-2,2.5e-2])

plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
legend()
axis([.5, 5.5,-2.5e-2,2.5e-2])

# Dr. Lin
fig, axs = subplots(1,2,sharey = true)
axs[1].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[1].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
axs[1].set_title("Direct Method")
axs[1].set_xlabel("frequencies")
axs[1].set_ylabel("Spectral density")
axs[1].grid("on")
axs[1].legend()

axs[2].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[2].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "sp2")
axs[2].set_title("Periodogram")
axs[2].set_xlabel("frequencies")
axs[2].set_ylabel("Spectral density")
axs[2].grid("on")
axs[2].legend()
fig.suptitle(" Estimated spectral densities over S_pred_plus")


## Causaul projection


S_sigpred_overS_plus_plus_num_fft_dm = h_wf_packs[6]
S_sigpred_overS_plus_plus_num_fft_sp = h_wf_packs[13]

nfft      = Parms[1][3]
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
legend()
title("Causal projection of Estimated spectral densities over S_pred_plus")
axis([.25, 6,-2.5e-2,2.5e-2])

fig, axs = subplots(1,2,sharey = true)
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
axs[1].set_title("Direct Method")
axs[1].set_xlabel("frequencies")
axs[1].set_ylabel("Spectral density")
axs[1].grid("on")
axs[1].axis([.25, 6,-2.5e-2,2.5e-2])
axs[1].legend()

axs[2].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[2].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
axs[2].set_title("Periodogram")
axs[2].set_xlabel("frequencies")
axs[2].set_ylabel("Spectral density")
axs[2].grid("on")
axs[2].axis([.25, 6,-2.5e-2,2.5e-2])
axs[2].legend()
fig.suptitle("Causal projection of Estimated spectral densities over S_pred_plus")



## Where things go wrong

fig, axs = subplots(2,2,sharey = true)

axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
axs[1].set_title("Direct Method")
axs[1].set_ylabel("Spectral density over S_pred_plus")
axs[1].grid("on")
axs[1].legend()

axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
axs[3].set_title("Periodogram")
axs[3].grid("on")
axs[3].legend()
#fig.suptitle(" Estimated spectral densities over S_pred_plus")

axs[2].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[2].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
axs[2].set_xlabel("frequencies")
axs[2].set_ylabel("Causal Projection of Spectral density over S_pred_plus")
axs[2].grid("on")
axs[2].axis([.25, 6,-2.5e-2,2.5e-2])
axs[2].legend()

axs[4].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[4].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
axs[4].set_xlabel("frequencies")
axs[4].set_ylabel("Spectral density")
axs[4].grid("on")
axs[4].axis([.25, 6,-2.5e-2,2.5e-2])
axs[4].legend()
fig.suptitle("Causal projection of Estimated spectral densities over S_pred_plus")


## Where things go wrong

fig, axs = subplots(4,2,sharey = "row",sharex = true)

axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
axs[1].set_title("Direct Method")
axs[1].set_ylabel("S_XY/S_X^+ (real)")
axs[1].axis([0, 6.28,-5e-2,2e-1])
axs[1].grid("on")
axs[1].legend()

axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
axs[2].set_ylabel("S_XY/S_X^+  (lmag)")
axs[2].grid("on")
axs[2].axis([0, 6.28,-2e-1,2e-1])

axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
axs[3].set_ylabel("{S_XY/S_X^+}_+  (real)")
axs[3].grid("on")
axs[3].axis([0, 6.28,-1e-2,2e-1])

axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
axs[4].set_ylabel("{S_XY/S_X^+}_+  (imag)")
axs[4].grid("on")
axs[4].axis([0, 6.28,-2e-1,2e-1])
axs[4].set_xlabel("frequencies")


axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
axs[5].set_title("Periodogram")
axs[5].grid("on")
axs[5].legend()

axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
axs[6].grid("on")

axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
axs[7].grid("on")

axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
axs[8].grid("on")
axs[8].set_xlabel("frequencies")
fig.suptitle("Causal projection of Estimated spectral densities over S_pred_plus")

## For Dr. Lin

fig, axs = subplots(4,2,sharey = true,sharex = true)

axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[1].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
axs[1].set_title("Direct Method")
axs[1].set_ylabel("Spectral density over S_pred_plus (real)")
axs[1].axis([.1, 6.2,-1e-1,1e-1])
axs[1].grid("on")
axs[1].legend()

axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[2].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
axs[2].set_ylabel("Spectral density over S_pred_plus (lmag)")
axs[2].grid("on")
axs[2].legend()

axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[3].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
axs[3].set_ylabel("Causal Projection of Spectral density over S_pred_plus (real)")
axs[3].grid("on")
axs[3].legend()

axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[4].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
axs[4].set_ylabel("Causal Projection of Spectral density over S_pred_plus (imag)")
axs[4].grid("on")
axs[4].set_xlabel("frequencies")
axs[4].legend()


axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[5].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
axs[5].set_title("Periodogram")
axs[5].grid("on")
axs[5].legend()

axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[6].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
axs[6].grid("on")
axs[6].legend()

axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[7].plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
axs[7].grid("on")
axs[7].legend()

axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[8].plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
axs[8].grid("on")
axs[8].set_xlabel("frequencies")
axs[8].legend()




















N=10
H_num_from_wf_dm = [fill(real(h_wf_dm[1,1,1]),N) fill(real(h_wf_dm[1,2,1]),N)]
plot(2pi*(0:N-1)/(N-1),H_num_from_wf_dm)
H_num_from_wf_sp = [fill(real(h_wf_sp[1,1,1]),N) fill(real(h_wf_sp[1,2,1]),N)]
plot(2pi*(0:N-1)/(N-1),H_num_from_wf_sp)




function spect_plot(S;
        plotter = plot,
        label = "unlabled")
    plotter(2pi*(0:length(S)-1)/length(S),S,
            label = label)
    legend()
end

H_num_dm = h_wf_packs[7]
H_num_sp = h_wf_packs[14]

nfft      = Parms[1][3]
plot(2pi*(0:nfft-1)/nfft,H_num_dm[1,1,:],label = "dm1")
plot(2pi*(0:nfft-1)/nfft,H_num_dm[1,2,:],label = "dm2")
plot(2pi*(0:nfft-1)/nfft,H_num_sp[1,1,:],label = "sp1")
plot(2pi*(0:nfft-1)/nfft,H_num_sp[1,2,:],label = "sp2")
legend()






i = 3

Parms[i]
Out = @timed get_wf(X, Psi;
    M_out, ty, rl, Preds, PI, rtol, info = true,
    xspec_est = Parms[i][1],
    par       = Parms[i][2],
    nfft      = Parms[i][3],
    n         = Parms[i][4],
    p         = Parms[i][5]);
Out.value[1][1,:,1]
z_crossspect_sigpred_num_fft_dm = Out.value[2]
H_num_dm = Out.value[7]
nfft = 2^10
plot(2pi*(0:nfft-1)/nfft,(real(z_crossspect_sigpred_num_fft_dm[1,1,:])),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,(real(z_crossspect_sigpred_num_fft_dm[1,2,:])),label = "dm2")

H_num_dm = Out.value[7]
nfft      = Parms[3][3]
plot(2pi*(0:nfft-1)/nfft,H_num_dm[1,1,:],label = "dm1")
plot(2pi*(0:nfft-1)/nfft,H_num_dm[1,2,:],label = "dm2")

S_sigpred_overS_plus_plus_num_fft_dm = Out.value[6]

nfft      = Parms[3][3]
plot(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:])),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:])),label = "dm2")

c = [real(h_wf_dm[1,1,1]) - real(h_wf_sp[1,1,1]);
     real(h_wf_dm[1,2,1]) - real(h_wf_sp[1,2,1])]


plot(2pi*(0:nfft-1)/nfft,real(S_C[1,1,:]),label = "S_C1")
plot(2pi*(0:nfft-1)/nfft,real(S_C[1,2,:]),label = "S_C2")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
legend()
axis([.5, 5.5,0,2.5e-2])

S_pred_plus = h_wf_packs[3]

S_C = complex(zeros(size(S_sigpred_overS_plus_plus_num_fft_sp)))
for o = 1:2^17
    S_C[:,:,o] = S_sigpred_overS_plus_plus_num_fft_sp[:,:,o] - c'*S_pred_minus[:,:,o]
end

H_num = complex(zeros(1,2,nfft))
for i = 1: nfft
    # matlog2[:,i] = svd(z_spect_pred_minus_num_fft[:,:,i]).S ###
    H_num[:,:,i] = S_C[:,:,i]/
                   S_pred_minus[:,:,i]
end
h_num_raw = ifft(H_num, 3)
