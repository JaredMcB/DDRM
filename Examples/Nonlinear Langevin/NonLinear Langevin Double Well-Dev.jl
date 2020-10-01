using PyPlot


include("DataGen.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")

steps = 10^4 + 1
scheme = "FE"
t_start = 0
t_stop = 10^3
discard = 100000
sig_init = [1.5]
sigma = [.3]
V_prime = x -> -x.*(x.^2 .- 1)
SM1 = false
Obs_noise = false
d = 1
e = randn(d,steps + discard)


Y = DataGen_DWOL(
    steps;
    scheme, t_start, t_stop, discard,
    sig_init , sigma, V_prime,
    SM1, Obs_noise, d, e
    )

T = range(t_start,stop = t_stop, length = steps)

Psi(x) = [x; x.^3]

# h_wf, pred = get_wf(Y,Psi, Preds = true);

# get_wf(Y,Psi, Preds = true);

# get_wf(
signal = Y  # Vector valued process
Psi; # column vector valued function
M_out = 20
n = 3
p = 1500
par = 1500
rl = true
Preds = false
PI = false
rtol = 1e-6
# )

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

# h_wf = vector_wiener_filter_fft(sig, pred, M_out,
#         n = n, p = p, par = par, PI = PI, rtol = rtol)

# vector_wiener_filter_fft(
sig
pred
M_out
par
win = "Par"
n
p
PI
rtol
# )

d, stepsy = size(sig)
nu, stepsx = size(pred)

stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])
nfft = nextfastfft(steps)
nffth = Int(floor(nfft/2))
L = 1500 #par

R_pred_smoothed = matrix_autocov_seq(pred,
   L = L,
   steps = steps,
   nu = nu,
   win = win
   )

lags = -L:L

# Smoothed viewing window
lam = _window(L, win = win, two_sided = false)

R_pred_smoothed = zeros(Complex,nu,nu,length(0:L))

### Break in
i = 1; j = 1
X = pred[i,1:steps]
C = my_crosscov(X,X,lags)



temp = my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
### break out
for i = 1 : nu
    for j = 1 : nu
        temp = my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
        temp = .5*(temp[L+1:end] + conj(reverse(temp[1:L+1])))
        R_pred_smoothed[i,j,:] = lam .* temp
    end
end
R_pred_smoothed


### Break in

plot(real([R_pred_smoothed[1,1,
                :] R_pred_smoothed[1,2,
                :] R_pred_smoothed[2,2,:]]))

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

### Break in
#recovered spectral density of pred
z_spect_pred_recoverd = complex(zeros(nu,nu,nfft))
for i = 1 : nfft
    z_spect_pred_recoverd[:,:,i] =z_spect_pred_minus_num_fft[:,:,i]'*
                                  z_spect_pred_plus_num_fft[:,:,i]
end

z_spect_pred = z_crossspect_fft(pred, pred,
                    nfft = nfft, n = n, p = p, win = "Par");

z_spect_pred_recoverd

semilogx(real(z_spect_pred_recoverd[1,1,
                :])) # z_spect_pred_recoverd[1,2,
                #:] z_spect_pred_recoverd[2,2,:]]))
# axis([0,2000,-20,100])

semilogx(real(z_spect_pred[1,1,
                :])) # z_spect_pred[1,2,
                #:] z_spect_pred[2,2,:]]))
# axis([0,2000,-20,300])

# Compute z-cross-spectrum of sigpred
z_crossspect_sigpred_num_fft = z_crossspect_fft(sig, pred,
                    nfft = nfft, n = n, p = p, win = "Par");

# This computes the impule response (coefficeints of z) for S_{yx}{S_x^+}^{-1}
S_sigpred_overS_plus_fft_num = complex(zeros(d,nu,nfft))

for i = 1 : nfft
    S_sigpred_overS_plus_fft_num[:,:,i] = z_crossspect_sigpred_num_fft[:,:,i]/
                                          z_spect_pred_plus_num_fft[:,:,i]
end

S_sigpred_overS_plus_fft_num_fft = ifft(S_sigpred_overS_plus_fft_num,3)

# Extracts causal part coefficinets of S_{yx}{S_x^+}^{-1}, {S_{yx}{S_x^+}^{-1}}_+
S_sigpred_overS_plus_fft_plus_num_fft = cat(dims = 3,
                S_sigpred_overS_plus_fft_num_fft[:,:,1: nffth],
                zeros(d,nu,nfft - nffth))

# Computes causal part of S_{yx}/S_x^+, {S_{yx}/S_x^+}_+
S_sigpred_overS_plus_plus_num_fft = fft(S_sigpred_overS_plus_fft_plus_num_fft,3);

# Obtain transfer function H by dividing {S_{yx}/S_x^+}_+ by S_x^-

H_num = complex(zeros(d,nu,nfft))
for i = 1: nfft
    H_num[:,:,i] = S_sigpred_overS_plus_plus_num_fft[:,:,i]/
                   z_spect_pred_minus_num_fft[:,:,i]
end

# Extrct tranferfunction coeffifcients (impulse responce of Weiner filter)
h_num_raw = ifft(H_num, 3)

# Truncate
M_out > nfft && println("M_out > nfft, taking min")
M = min(M_out, nfft)
h_num_fft = h_num_raw[:,:,1:M]

h_wf = rl ? real(h_wf) : h_wf
Preds ? [h_wf, pred] : h_wf

Δt = (t_stop - t_start)/(steps - 1)
h_ana = zeros(1,2,10)
h_ana[:,:,1] = [1+Δt -Δt]

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
