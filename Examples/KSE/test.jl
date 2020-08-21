using JLD, Dates

include("Model_Reduction_Dev.jl")
pwd()
Dat = load("..\\Server_scripts\\Data\\KSE_sol3.jld")

uu = Dat["dat_uu"]
vv = Dat["dat_vv"]
tt = Dat["dat_tt"]
h = Dat["h"]
d = 5
h = Dat["h"]
N = Dat["N"]
P = Dat["P"]
obs_gap = Dat["obs_gap"]
Δt = h*obs_gap

V_obs = vv[2:d+1,1:end]


# Build PSI
function InvBurgRK4_1step(x)
 lx = length(x)
 function F(x)
     𝑥 = [conj(reverse(x, dims = 1));0; x]
     -im/2*(2π/P*(1:lx)/N) .* conv(𝑥,𝑥)[2*lx+2:3*lx+1]
 end

 Δt = h*obs_gap

 k1 = F(x)
 k2 = F(x .+ Δt*k1/2)
 k3 = F(x .+ Δt*k2/2)
 k4 = F(x .+ Δt*k3)
 A =  @. x + Δt/6*(k1 + 2k2 + 2k3 + k4)
end

function Inertialman_part(x)
  lx = length(x)
  𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

  L = complex(zeros(lx^2))
  for j = 1:lx
     for k = 1:lx
        L[ (j-1)*lx+k] = 𝑥(j+lx)*𝑥(j+lx-k)
     end
  end
  L
end

function Inertialman_part_short(x)
  lx = length(x)
  𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

  L = complex(zeros(binomial(lx+1,2)))
  i = 1
  for j = 1:lx
     for k = j:lx # k should normaly go from 1 to lx but i changed it to go from j to lx.
        L[i] = 𝑥(j+lx)*𝑥(j+lx-k)
        i += 1
     end
  end
  L
end

Psi(x) = [x; InvBurgRK4_1step(x); Inertialman_part_short(x)]

M_out = 100

@time h_wf = get_wf(V_obs,Psi,
     M_out = M_out,
     rl = false,
     PI = false)

## Expanded
signal = V_obs
sig = V_obs[:,2:end]
d, steps = size(sig)
nu = size(Psi(zeros(d,1)),1)

pred = complex(zeros(nu, steps))
for n = 1:steps
   pred[:,n] = Psi(signal[:,n])
end

# vector_wiener_filter_fft(pred, sig, M_out,
                 # PI = PI, rtol = rtol)
# Expand

# function vector_wiener_filter_fft(pred::Array{Complex{Float64},2}, sig,
#     M_out = 20;
#     par::Int64 = 1500,
#     win = "Par",
#     n = 3,
#     p = 2500,
#     PI = true,
#     rtol = 1e-6
#     )

M_out = 20
par = 1500
win = "Par"
n = 3
p = 2500
PI = true
rtol = 1e-6




d, stepsy = size(sig)
nu, stepsx = size(pred)

stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])
nfft = nextfastfft(steps)
nffth = Int(floor(nfft/2))

L = par
lags = 0:L;

# Smoothed viewing window
lam = _window(L, win = win, two_sided = false)

R_pred_smoothed = zeros(Complex,nu,nu,length(lags))
for i = 1 : nu
    for j = 1 : nu
        R_pred_smoothed[i,j,:] = lam .* my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
    end
end

# Compute coefficients of spectral factorization of z-spect-pred
l = spectfact_matrix_CKMS_pinv(R_pred_smoothed)

l_pad_minus = nfft >= L+1 ? cat(dims = 3,l,zeros(nu,nu,nfft - L - 1)) :
                           l[:,:,1:nfft]

z_spect_pred_minus_num_fft = fft(l_pad_minus,3)
z_spect_pred_plus_num_fft = zeros(Complex,nu,nu,nfft)
for i = 1 : nfft
    z_spect_pred_plus_num_fft[:,:,i] = z_spect_pred_minus_num_fft[:,:,i]'
end

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
h_wf = h_num_fft
dat = Dict("dat_h_wf" => h_wf)
save(("..\\Server_scripts\\Data\\KSE_wf_New.jld"),dat)
