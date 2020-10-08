include("modgen_LSDE.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")

using JLD
using PyPlot
using DSP: nextfastfft
using Distributions

A = reshape([-0.5],1,1)
σ = reshape([1],1,1)
Xo = [1]
t_disc = 1000
gap = 1
scheme = "EM"

t_start = 0
t_stop  = 1e4
h       = 1e-2

Δt      = h*gap
M_out   = 20

n = 2
p = 150

X = modgen_LSDE(t_start,t_stop,h,
    A = A,
    σ = σ,
    Xo = Xo,
    t_disc = t_disc,
    gap = gap)

d, N = size(X)

nfft = nextfastfft(N)
X = [X zeros(d,nfft-N)]


τ_exp, τ_int    = auto_times(X[:])*Δt
N_eff           = N*Δt/τ_int

Psi(x) = x

n = 2
p = 150
ty = "bin"
par = 2000

@time h_wf = get_wf(X, Psi; par,n,p,ty);
X = X[:,1:N]

h_wf

signal = X # Vector valued process
Psi # column vector valued function

M_out = 20
rl = true
Preds = false
PI = false
rtol = 1e-6

# We would like a presample since we want the
# times series to be offset by one.

sig = signal[:,2:end] # sig is now one a head of signal
d, steps = size(sig)
nu = size(Psi(zeros(d,1)),1)

pred = complex(zeros(nu, steps))
for n = 1:steps
    pred[:,n] = Psi(signal[:,n])
end # pred is now even with signal and therefore one step
    # step behind sig. I.e. pred[:,n] = Psi(sig[:,n-1])
    # which is what we want so as to ensure the reduced
    # model can run explicitly.

win = "Par"

d, stepsy = size(sig)
nu, stepsx = size(pred)

stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])
nfft = nextfastfft(steps)
nffth = Int(floor(nfft/2))
L = par

R_pred_smoothed = matrix_autocov_seq(pred,
   L = L,
   steps = steps,
   nu = nu,
   win = win
   )

plot(real(R_pred_smoothed[1,1,:]))

# Compute coefficients of spectral factorization of z-spect-pred
l = PI ? spectfact_matrix_CKMS_pinv(R_pred_smoothed,rtol = rtol) :
         spectfact_matrix_CKMS(R_pred_smoothed)

l_pad_minus = nfft >= L+1 ? cat(dims = 3,l,zeros(nu,nu,nfft - L - 1)) :
                           l[:,:,1:nfft]

z_spect_pred_minus_num_fft = fft(l_pad_minus,3)
z_spect_pred_plus_num_fft = complex(zeros(nu,nu,nfft))
for i = 1 : nfft
    z_spect_pred_plus_num_fft[:,:,i] = z_spect_pred_minus_num_fft[:,:,i]'
end

#recovered spectral density of pred
z_spect_pred_recoverd = complex(zeros(nu,nu,nfft))
for i = 1 : nfft
    z_spect_pred_recoverd[:,:,i] =z_spect_pred_minus_num_fft[:,:,i]*
                                  z_spect_pred_plus_num_fft[:,:,i]
end

z_spect_pred = z_crossspect_fft(pred, pred,
                    nfft = nfft, n = 2, p = 150, ty = "bin");

semilogy(2π*(1:nfft)/nfft,[real(z_spect_pred[1,1,
                :]) real(z_spect_pred_recoverd[1,1,
                :])]) # z_spect_pred[1,2,
                #:] z_spect_pred[2,2,:]]))
# axis([0,2000,-20,300])

# Compute z-cross-spectrum of sigpred
z_crossspect_sigpred_num_fft = z_crossspect_fft(sig, pred;
                    nfft, n, p, ty);

Theta = 2π*(1:nfft)/nfft



uz_spect_pred = real([exp(-Theta[i]*im)*z_spect_pred[1,1,i] for i in 1:nfft])

semilogy(2π*(1:nfft)/nfft,[real(z_crossspect_sigpred_num_fft[1,1,:]) uz_spect_pred])



(τ_exp, N_eff)

1 .+ h*A

h_wf

# d, N  = size(X)
# nu    = size(Psi(X[:,1]),1)
# M_out = size(h_wf,3)

# X_rm = zeros(d,N); X_rm[:,1:M_out] = X[:,1:M_out]

# PSI = zeros(nu,N);
# for i = 1:M_out
#     PSI[:,i] = Psi(X_rm[:,i])
# end

# for i = M_out + 1 : N
#     X_rm[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out, dims = 2)
#     PSI[:,i] = Psi(X_rm[:,i])
# end



d, N  = size(X)
nu    = size(Psi(X[:,1]),1)
M_out = 2 #size(h_wf,3)

X_rm = zeros(d,N); X_rm[:,1:M_out] = X[:,1:M_out]

PSI = zeros(nu,N);
for i = 1:M_out
    PSI[:,i] = Psi(X_rm[:,i])
end

for i = M_out + 1 : N
    X_rm[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out, dims = 2) + sqrt(h)*σ*randn(d)
    PSI[:,i] = Psi(X_rm[:,i])
end

X_rm

data = Dict(
        "h_wf" => h_wf,
        "A" => A,
        "σ" => σ,
        "Xo" => Xo,
        "t_disc" => t_disc,
        "gap" => gap,
        "scheme" => scheme,
        "t_start" => t_start,
        "t_stop" => t_start,
        "h" => h,
        "X_55" => X,
        "X__rm_55_h_4" => X_rm)

save("LSDE_Data\\LSDE_wfs_M$M_out.jld",data)

data = load("LSDE_Data\\LSDE_wfs_M10.jld")

blup = findall(isnan,X_rm[1,:])[1]

findall(x -> x > 10^1,X_rm[1,:])[1]

plot([X[1:1000:end] X_rm[1:1000:end]])

m = mean(X[1:9900001])
m_rm = mean(X_rm[1:9900001])
v = var(X[1:9900001])
v_rm = var(X_rm[1:9900001])

dist_rm = Normal(m_rm,v_rm)

f_rm(x) = pdf(dist_rm, x)

NN = 1000
XX = -5:1/NN:5

plot(XX,f.(XX),":",label = "orig normal")
plot(XX,f_rm.(XX),":",label = "rm normal")

emp_pdf(X[1:9900001]);
emp_pdf(X_rm[1:9900001]);

legend()

lags = 0:1000

A_rm = my_autocor(X_rm[:],lags)
A    = my_autocor(X[:],lags)

semilogy(lags*h,[A A_rm])

z_spect = z_spect_scalar(X[:], n = 3, p=100, ty = "ave")
z_spect_rm = z_spect_scalar(X_rm[:], n = 3, p=100, ty = "ave")

Θ = 2*π*(1:1000:nfft)/nfft
Z = exp.(im*Θ)

a = 1 + h*A[1,1]
z_spect_ana_fun(z) = h*σ/( (1-a*z^(-1))*(1-a*z) )
z_spect_ana = real(z_spect_ana_fun.(Z))

semilogy(Θ,[z_spect[1:1000:nfft] z_spect_rm[1:1000:nfft] z_spect_ana])
