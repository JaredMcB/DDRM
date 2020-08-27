using PyPlot
using JLD
using Dates
using FFTW

include("c:\\Users\\JaredMcBride\\Desktop\\"*
            "DDMR\\Tools\\AnalysisToolbox.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")
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

NN_ckms = map(x-> floor(Int, 10^x),2:.25:4)
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


semilogy(NN_ckms, Norm)
