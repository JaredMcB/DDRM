include("..\\Model_Reduction_Dev.jl")
include("AnalysisToolbox_scratch_ckms.jl")


using PyPlot
using JLD

function visual_test_ckms(P,l,nfft;semilog = false)
    d  = size(P,1)
    lp = size(P,3)
    ll = size(l,3)
    S_fun(z)    = P[:,:,1] + sum(P[:,:,i]*z^(-i+1) + P[:,:,i]'*z^(i-1) for i = 2:lp)
    S_fun_minus(z) = sum(l[:,:,i]*z^(-i+1) for i = 1:ll)
    S_fun_plus(z) = sum(l[:,:,i]'*z^(i-1) for i = 1:ll)

    Θ = 2π*(0:nfft-1)/nfft
    Z = exp.(im*Θ)
    S = complex(zeros(d,d,nfft))
    S_l = complex(zeros(d,d,nfft))
    for i = 1:nfft
        S[:,:,i] = S_fun(Z[i])
        S_l[:,:,i] = S_fun_minus(Z[i])*S_fun_plus(Z[i])
    end


    for i = 1:d
        for j = i:d
            semilog ? semilogy(Θ,real(S[i,j,:]), label = "S ($i,$j)") :
                      plot(Θ,real(S[i,j,:]), label = "S ($i,$j)")

            semilog ? semilogy(Θ,real(S_l[i,j,:]), label = "S_l ($i,$j)") :
                      plot(Θ,real(S_l[i,j,:]), label = "S_l ($i,$j)")
        end
    end
    legend()
end

P = zeros(1,1,2)
P[1,1,1] = 10
P[1,1,2] = 3

L = spectfact_matrix_CKMS(P)

d, d1, m = size(L)

steps = 10^6

u = randn(d, 2*steps)
X = zeros(d, 2*steps)
for n = m:2*steps
    X[:,n] = sum(L[:,:,j]*u[:,n-j+1] for j = 1:m)
end
X = X[:,steps+1:end]

P_autocov = matrix_autocov_seq(X,
       L = 1500,
       steps = steps,
       nu = d,
       win = "Par"
       )

z_spect_X = z_spect_scalar(X[1,:], n = 3, p = 1500, ty = "ave")
gap = 100
Θ = 2π*(1:gap:steps)/steps
plot(Θ,z_spect_X[1:gap:steps])
nfft = 1000
visual_test_ckms(P,L,nfft)

P_fft_full = fft(z_spect_X)/steps;

P_fft = reshape(P_fft_full[1:2],1,1,2)

L_fft = spectfact_matrix_CKMS(P_fft)

L_autocov = spectfact_matrix_CKMS(P_autocov)

visual_test_ckms(P_autocov,L_autocov,nfft)



## 2×2 Matrix
P = zeros(2,2,2)
P[:,:,1] = [174 17 ; 17 42]
P[:,:,2] = [-2 -4; -79 -8]

l_ana = zeros(2,2,2)
l_ana[:,:,1] = [1 13; 2 1]
l_ana[:,:,2] = [-2 0; -1 -6]

l_num = spectfact_matrix_CKMS(P)

d, d1, m = size(l_ana)

steps = 10^6

u = randn(d, 2*steps)
X = zeros(d, 2*steps)
for n = m:2*steps
    X[:,n] = sum(l_ana[:,:,j]*u[:,n-j+1] for j = 1:m)
end
X = X[:,steps+1:end]

P_autocov = matrix_autocov_seq(X,
       L = 1500,
       steps = steps,
       nu = d,
       win = "Par"
       )

L_autocov = spectfact_matrix_CKMS(P_autocov)
visual_test_ckms(P_autocov,L_autocov,nfft)
visual_test_ckms(P,l_ana,nfft)

d = 4; m = 10
P = zeros(d,d,m)
for i = 2 : m
    b = randn(d,d)
    P[:,:,i] = b
end
b = randn(d,d)
P[:,:,1] = b*b' + sum(P[:,:,i]*P[:,:,i]' for i = 2:m) + 1e-10*I;

d  = size(P,1)
lp = size(P,3)
S_fun(z)    = P[:,:,1] + sum(P[:,:,i]*z^(-i+1) + P[:,:,i]'*z^(i-1) for i = 2:lp)

Θ = 2π*(0:nfft-1)/nfft
Z = exp.(im*Θ)
detS = zeros(nfft)
for i = 1:nfft
    detS[i] = real(det(S_fun(Z[i])))
end
plot(Θ,detS)

l_num, Err = spectfact_matrix_CKMS_SC(P, N_ckms = 1000, update = 100);

d, d1, m = size(l_num)

steps = 10^6

u = randn(d, 2*steps)
X = zeros(d, 2*steps)
for n = m:2*steps
    X[:,n] = sum(l_num[:,:,j]*u[:,n-j+1] for j = 1:m)
end
X = X[:,steps+1:end]

P_autocov = matrix_autocov_seq(X,
       L = 55,
       steps = steps,
       nu = d,
       win = "Par"
       )

L_autocov = spectfact_matrix_CKMS(P_autocov)

visual_test_ckms(P_autocov,L_autocov,nfft)
visual_test_ckms(P,l_num,nfft)
legend(bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)



## Kse Data
data_loc = "C:\\Users\\jared\\Desktop\\Github Repos"*
    "\\DDMR\\Examples\\KSE\\Data\\KSE_sol_Lin.jld"

Data = load(data_loc)

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

X = vv[2:4,:]
d, steps = size(X)

L = 2000
d = 2
P_autocov = matrix_autocov_seq(X,
       L = L,
       steps = steps,
       nu = d,
       win = "Par"
       )

v2 = vv[3,:]
lags = -L:L

A_v2 = my_crosscov(v2,v2,lags)
plot(lags*Δt,A_v2)
plot(Δt*(0:L),P_autocov[2,2,:])

L_autocov = spectfact_matrix_CKMS(P_autocov)

nfft = 1000
visual_test_ckms(P_autocov[1:1,1:1,:],L_autocov[1:1,1:1,:],nfft, semilog = true)

z_spect_fft = z_crossspect_fft(
    X,
    X,
    nfft = 0,
    n = 3,
    p = 2500,
    win = "Par")

Nfft = size(z_spect_fft,3)
Θ_N = 2π*(0:Nfft-1)/Nfft
d=1
for i = 1:d
    for j = i:d
        semilogy(Θ_N,z_spect_fft[i,j,:], label = "S_N ($i,$j)")
    end
end

legend(bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)

z_spect_fft
