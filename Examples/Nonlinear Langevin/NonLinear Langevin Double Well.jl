using PyPlot
using Random


include("DataGen.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")

steps = 10^6 + 1
scheme = "FE"
t_start = 0
t_stop = 10^4
discard = 100000
sig_init = [1.5]
sigma = [.3]
V_prime = x -> -x.*(x.^2 .- 1)
SM1 = false
Obs_noise = false
d = 1
#e = randn(d,steps + discard)

# Get full model run
Random.seed!(2014)
X = DataGen_DWOL(
    steps;
    scheme, t_start, t_stop, discard,
    sig_init , sigma, V_prime,
    SM1, Obs_noise, d
    )

T = range(t_start,stop = t_stop, length = steps)
X

data = Dict("steps" => steps,
            "t_stop" => t_start,
            "sigma" => sigma,
            "X" => X)
save("Examples/Nonlinear Langevin/data/data_10_23_2020.jld",data)

data = load("Examples/Nonlinear Langevin/data/data_10_23_2020.jld")
X = data["X"]

auto_times(X[1,:])

# Put in Psi functions
Psi(x) = [x; x.^3]

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
# Parms = [["DM"       , 100    , 2^10    , 3    , 500],
#          ["DM"       , 500    , 2^16    , 3    , 500],
#          ["DM"       , 1000   , 2^16    , 3    , 500],
#          ["DM"       , 5000   , 2^16    , 3    , 500],
#          ["DM"       , 10000  , 2^16    , 3    , 500],
#          ["DM"       , 5000   , 2^20    , 3    , 500],
#          ["DM"       , 10000  , 2^20    , 3    , 500],
#          ["SP"       , 10000  , 2^10    , 2    , 5],
#          ["SP"       , 10000  , 2^16    , 2    , 5],
#          ["SP"       , 10000  , 2^17    , 2    , 5],
#          ["SP"       , 10000  , 2^17    , 3    , 10]]

Parms = [["DM"       , 5000  , 2^17    , 2    , 5],
         ["SP"       , 5000  , 2^17    , 2    , 5]]

P = length(Parms)

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
h_wf_sp = h_wf_packs[9]

h_wf_dm[1,:,1]
h_wf_sp[1,:,1]

function spect_plot(S;
        plotter = plot,
        label = "unlabled")
    plotter(2pi*(0:length(S)-1)/length(S),S,
            label = label)
    legend()
end

H_num_dm = h_wf_packs[8]
H_num_sp = h_wf_packs[16]

nfft      = Parms[1][3]
plot(2pi*(0:nfft-1)/nfft,H_num_dm[1,1,:],label = "dm1")
plot(2pi*(0:nfft-1)/nfft,H_num_dm[1,2,:],label = "dm2")
plot(2pi*(0:nfft-1)/nfft,H_num_sp[1,1,:],label = "sp1")
plot(2pi*(0:nfft-1)/nfft,H_num_sp[1,2,:],label = "sp2")
legend()


S_sigpred_overS_plus_plus_num_fft_dm = h_wf_packs[7]
S_sigpred_overS_plus_plus_num_fft_sp = h_wf_packs[15]

semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:])),label = "dm1")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:])),label = "dm2")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:])),label = "sp1")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:])),label = "sp2")
legend()
axis([.5, 5.5,-1e-1,1e1])


plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
legend()
axis([.5, 5.5,0,2.5e-2])

plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2")
legend()
axis([.5, 5.5,0,2.5e-2])



S_sigpred_overS_plus_fft_num_dm = h_wf_packs[6]
S_sigpred_overS_plus_fft_num_sp = h_wf_packs[14]

semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_dm[1,1,:])),label = "dm1")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_dm[1,2,:])),label = "dm2")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_sp[1,1,:])),label = "sp1")
semilogy(2pi*(0:nfft-1)/nfft,(real(S_sigpred_overS_plus_fft_num_sp[1,2,:])),label = "sp2")
legend()
axis([.5, 5.5,-1e-1,1e1])


plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
legend()
axis([.5, 5.5,0,2.5e-2])

plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
plot(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2")
legend()
axis([.5, 5.5,0,2.5e-2])






## h = 0.1

h_wf[1,:,1,:]'

steps = 10^6 + 1
scheme = "FE"
t_start = 0
t_stop = 10^5
discard = 100000
sig_init = [1.5]
sigma = [.3]
V_prime = x -> -x.*(x.^2 .- 1)
SM1 = false
Obs_noise = false
d = 1
#e = randn(d,steps + discard)

# Get full model run
Random.seed!(2014)
X2 = DataGen_DWOL(
    steps;
    scheme, t_start, t_stop, discard,
    sig_init , sigma, V_prime,
    SM1, Obs_noise, d
    )

T = range(t_start,stop = t_stop, length = steps)

data = Dict("steps" => steps,
            "t_stop" => t_start,
            "sigma" => sigma,
            "X2" => X2)
save("data/data_10_28_2020.jld",data)

data = load("data/data_10_28_2020.jld")
X2 = data["X2"]

auto_times(X2[1,:])

###         Varing parameters
###            xspect_est, par    , nfft    , n    , p

# Parms = [["DM"       , 100    , 2^10    , 3    , 500],
#          ["DM"       , 500    , 2^16    , 3    , 500],
#          ["DM"       , 1000   , 2^16    , 3    , 500],
#          ["DM"       , 5000   , 2^16    , 3    , 500],
#          ["DM"       , 10000  , 2^16    , 3    , 500],
#          ["DM"       , 5000   , 2^20    , 3    , 500],
#          ["DM"       , 10000  , 2^20    , 3    , 500],
#          ["SP"       , 10000  , 0       , 3    , 500],
#          ["SP"       , 10000  , 10^6    , 2    , 5000],
#          ["SP"       , 10000  , 2^16    , 3    , 500]]

P = length(Parms)

h_wf  = zeros(1,2,M_out,P)
times = zeros(P)
for i = 1:P
    Out = @timed get_wf(X2, Psi;
        M_out, ty, rl, Preds, PI, rtol,
        xspec_est = Parms[i][1],
        par       = Parms[i][2],
        nfft      = Parms[i][3],
        n         = Parms[i][4],
        p         = Parms[i][5]);

    h_wf[:,:,:,i] = Out.value
    times[i]      = Out.time
end

h_wf[1,:,1,:]'

Parms_now = [["DM"       , 10000  , 2^20    , 3    , 500],
             ["SP"       , 10000  , 10^6    , 2    , 2000]]

P = length(Parms_now)

h_wf  = zeros(1,2,M_out,P)
times = zeros(P)
for i = 1:P
    Out = @timed get_wf(X2, Psi;
        M_out, ty, rl, Preds, PI, rtol,
        xspec_est = Parms_now[i][1],
        par       = Parms_now[i][2],
        nfft      = Parms_now[i][3],
        n         = Parms_now[i][4],
        p         = Parms_now[i][5]);

    h_wf[:,:,:,i] = Out.value
    times[i]      = Out.time
end

h_wf[1,:,1,:]'

Parms_now = [["DM"       , 10000  , 10^6    , 3    , 500],
             ["SP"       , 10000  , 10^6    , 2    , 2000]]

P = length(Parms_now)

h_wf  = zeros(1,2,M_out,P)
times = zeros(P)
for i = 1:P
    Out = @timed get_wf(X2, Psi;
        M_out, ty, rl, Preds, PI, rtol,
        xspec_est = Parms_now[i][1],
        par       = Parms_now[i][2],
        nfft      = Parms_now[i][3],
        n         = Parms_now[i][4],
        p         = Parms_now[i][5]);

    Out_info = Out.value
    times[i]      = Out.time
end

h_wf[1,:,1,:]'

X_sig = X2[:,2:end];

Psi(x) = [x; x.^3]
X_pred = get_pred(X,Psi) # Notice it is just
                         # X get_pred assigns
                         # psi straight across

# Model reduction Parameters
L = Parms_now[1][2]
Nex = Parms_now[1][3]

S_dm = z_crossspect_fft_old(X_sig, X_pred; L, Nex);

n = Parms_now[2][4]
p = Parms_now[2][5]
ty = "bin"
nfft = Parms_now[2][3]

S_sp = z_crossspect_fft(X_sig, X_pred; nfft, n, p, ty);

N = size(S_sp,3)
semilogx(2π*(0:N-1)/N,abs.(S_sp[1,:,:]'),label = "sp")
N = size(S_dm,3)
semilogx(2π*(0:N-1)/N,abs.(S_dm[1,:,:]'),":",label = "dm")
legend()

Parms_now = [["DM"       , 10000  , 10^6    , 3    , 500],
             ["SP"       , 10000  , 10^6    , 2    , 2000]]

P = length(Parms_now)
i=1
Out_DM = get_wf(X2, Psi;
        M_out, ty, rl, Preds, PI, rtol,
        xspec_est = Parms_now[i][1],
        par       = Parms_now[i][2],
        nfft      = Parms_now[i][3],
        n         = Parms_now[i][4],
        p         = Parms_now[i][5],info = true);

i = 2
Out_sp = get_wf(X2, Psi;
        M_out, ty, rl, Preds, PI, rtol,
        xspec_est = Parms_now[i][1],
        par       = Parms_now[i][2],
        nfft      = Parms_now[i][3],
        n         = Parms_now[i][4],
        p         = Parms_now[i][5],info = true);

H_sp = Out_sp[8]

H_dm = Out_DM[8]

plot(2π*(0:nfft-1)/nfft,H_dm[1,:,:]',label = "dm")
legend()

plot(2π*(0:nfft-1)/nfft,H_sp[1,:,:]',label = "sp")
legend()

semilogx(2π*(0:nfft-1)/nfft,H_dm[1,1,:],".",label = "dm")
semilogx(2π*(0:nfft-1)/nfft,H_sp[1,1,:],label = "sp")
legend()


semilogx(2π*(0:nfft-1)/nfft,H_dm[1,2,:],".",label = "dm")
semilogx(2π*(0:nfft-1)/nfft,H_sp[1,2,:],label = "sp")
legend()













M_h = size(h_ana,3)

Y_hat = zeros(size(Y));
Y_hat[:,1:M_h] = Y[:,1:M_h]
for i=M_h:steps-1
    Y_hat[:,i] = sum(h_ana[:,:,k+1]*pred[:,i-k]
                    for k = 0:M_h-1)
end

err = Y - Y_hat

Lags = -100:10
C1 = my_crosscov(pred[1,:],err[:],Lags)
C2 = my_crosscov(pred[2,:],err[:],Lags)


plot(Lags,real([C1 C2]))
