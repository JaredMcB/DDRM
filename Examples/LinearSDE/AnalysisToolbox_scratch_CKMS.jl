using PyPlot
using JLD
using Dates
using FFTW

include("c:\\Users\\JaredMcBride\\Desktop\\"*
            "DDMR\\Tools\\AnalysisToolbox.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")
include("AnalysisToolbox_scratch.jl")

Data = load("c:\\Users\\JaredMcBride\\Desktop\\"*
            "DDMR\\Examples\\KSE\\Data\\KSE_sol_lin.jld")

uu = Data["dat_uu"]
vv = Data["dat_vv"]
tt = Data["dat_tt"]

h = Data["h"]
N = Data["N"]
P = Data["P"]
obs_gap = Data["obs_gap"]
Δt = h*obs_gap

t_start = 0
t_stop = 150
ind_start = floor(Int,t_start/Δt)+1
ind_stop =floor(Int,t_stop/Δt)

H1 = imshow(uu[:,ind_start:ind_stop]', extent=[0,21.55,0,150], aspect="auto")

pred = vv[5:8,:]
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

LL, Err = spectfact_matrix_CKMS_SC(R_pred_smoothed,
    ϵ = 1e-8)

Err
semilogy(Err)
title("2x2 - Matrix function factorization")
xlabel("No. of interations")
ylabel("∞-norm of one-step difference")





plot((0:L)*Δt,real(R_pred_smoothed[1,1,:]))

expons = 1:.25:4.5
Num_e = length(expons)
NN_ckms = map(x-> floor(Int, 10^x),expons)
LL = complex(zeros(nu,nu,L+1,Num_e))
for i in 1:Num_e
    LL[:,:,:,i] = spectfact_matrix_CKMS_SC(R_pred_smoothed,
        N_ckms = NN_ckms[i])
end

Norm = zeros(Num_e,0)
for i = 1:nu
    for j = 1:nu
        global Norm
        Norm = [Norm map(k -> norm(LL[i,j,:,Num_e] .- LL[i,j,:,k],Inf),1:Num_e)]
    end
end
Norm

semilogy(NN_ckms, Norm)



## Test a random function
nu = 2  # mension of the square matrix valued function
nfft = 10^4  # No. of fft grid points
dL_ana = randn(nu,nu,nfft-1)  # White noise
L_ana = cumsum(cat(dims = 3, zeros(nu,nu), dL_ana),dims = 3) # Brownian moition

BB_ana = zeros(nu,nu,nfft) # Declare Brownian Bridge
S_ana = zeros(nu,nu,nfft) # Declasr the function
B = rand(nu,nu) # some shift to get the function off zero.

for i = 1:nfft
    # Form braonian bridge
    BB_ana[:,:,i] = L_ana[:,:,i] - (i-1)/(nfft-1) * L_ana[:,:,nfft] + B
    # Craft posdef valued function
    S_ana[:,:,i] = BB_ana[:,:,i]'*BB_ana[:,:,i]
end

plot(2π*(0:nfft-1)/(nfft-1),S_ana[1,1,:])

# Low pass filter
L = 199
S_ana_cos = fft(S_ana,3)

S_ana_coef_plus = S_ana_cos[:,:,1:L+1]
S_ana_coef_minus = complex(zeros(nu,nu,L))
for i = 1:L
    S_ana_coef_minus[:,:,i] = S_ana_coef_plus[:,:,L+2-i]'
end


S_ana_trunc = ifft(cat(S_ana_coef_plus,
                        zeros(nu,nu,nfft - 2L-1),
                        S_ana_coef_minus, dims = 3),3)
plot(2π*(0:nfft-1)/(nfft-1),S_ana_trunc[1,1,:])

# Spectral factor
LL, Err = spectfact_matrix_CKMS_SC(S_ana_cos[:,:,1:L+1],
    ϵ = 1e-8)

semilogy(Err)

LL
LL_pad_minus = cat(dims = 3,LL,zeros(nu,nu,nfft - L - 1))
S_num_minus = ifft(LL_pad_minus,3)
S_num_plus = zeros(Complex,nu,nu,nfft)
for i = 1 : nfft
    S_num_plus[:,:,i] = S_num_minus[:,:,i]'
end

S_num = complex(zeros(nu,nu,nfft))
for i = 1:nfft
    S_num[:,:,i] = S_num_minus[:,:,i]*S_num_plus[:,:,i]
end
S_num
plot(S_num[1,1,:])



BB_num = ifft(cat(LL,zeros(nu,nu,nfft - 1501), dims = 3))

plot([BB_num[1,1,:] BB_ana[1,1,:]])

S_num = complex(zeros(nu,nu,nfft))
for i = 1:nfft
    S_num[:,:,i] = BB_num[:,:,i]'*BB_num[:,:,i]
end

plot([S_num[1,1,:] S_ana[1,1,:]])



## Test a normal function
nu = 2  # mension of the square matrix valued function
nfft = 10^4  # No. of fft grid points

S_ana_fun(z) = [2z^(-1)+6+2z 11z^(-1)+22+7z; 7z^(-1)+22+11z 38z^(-1)+84+38z]
S_ana_fun(1)

Θ = 2π*(0:nfft-1)/nfft
Z = exp.(im*Θ)

S_ana = complex(zeros(nu,nu,nfft))
for i = 1:nfft
    S_ana[:,:,i] = S_ana_fun(Z[i])
end

plot(2π*(0:nfft-1)/nfft,S_ana[1,1,:])

# Low pass filter
L = 199
S_ana_fft = fft(S_ana,3)

S_ana_coef_plus = S_ana_fft[:,:,1:L+1]
S_ana_coef_minus = complex(zeros(nu,nu,L))
for i = 1:L
    S_ana_coef_minus[:,:,i] = S_ana_coef_plus[:,:,L+2-i]'
end


S_ana_trunc = ifft(cat(S_ana_coef_plus,
                        zeros(nu,nu,nfft - 2L-1),
                        S_ana_coef_minus, dims = 3),3)
plot(2π*(0:nfft-1)/(nfft-1),S_ana_trunc[1,1,:])

S_ana_cos = fft(S_ana,3)/nfft
S_ana_coef_input = complex(zeros(nu,nu,L+1))
S_ana_coef_input[:,:,1] = S_ana_cos[:,:,1]
for l = 1:L
    S_ana_coef_input[:,:,l+1] = S_ana_cos[:,:,end + 1 - l]
end
S_ana_coef_input
# Spectral factor
LL, Err = spectfact_matrix_CKMS_SC(S_ana_coef_input[:,:,1:L+1],
    ϵ = 1e-8)

semilogy(Err)
loglog(Err)

LL
LL_pad_minus = cat(dims = 3,LL,zeros(nu,nu,nfft - L - 1))
S_num_minus = ifft(LL_pad_minus,3)
S_num_plus = zeros(Complex,nu,nu,nfft)
for i = 1 : nfft
    S_num_plus[:,:,i] = S_num_minus[:,:,i]'
end

S_num = complex(zeros(nu,nu,nfft))
for i = 1:nfft
    S_num[:,:,i] = S_num_minus[:,:,i]*S_num_plus[:,:,i]
end
S_num
plot(S_num[1,1,:])





S = zeros(2,2,10)
S[:,:,1] = [6 22; 22 84]
S[:,:,2] = [2 11; 7 38]

l = spectfact_matrix_CKMS_pinv(S)
