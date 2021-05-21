# # Lasso Playspace
using JLD
using DSP # For conv function in Psi
using Dates
using GLMNet

using PyPlot

mrl = include("../../Tools/WFMR_lasso.jl")
mr  = include("../../Tools/WFMR.jl")
at  = include("../../Tools/AnalysisToolbox.jl")

# Load Old Data

# +
gen = "lin1e3"     # this is just a reference designation it shows up in the
                # output file. I think of generatrion.

server = startswith(pwd(), "/u5/jaredm") ? true : false
println("on server = $server")

sol_file = server ? "../../../data/KSE_Data/ks_sol_$gen.jld" :
   "C:/Users/JaredMcBride/Desktop/DDMR/Examples/KSE/Data/ks_sol_$gen.jld"
println("Sol load location: " * sol_file)

@time vv = load(sol_file,"dat_vv");
# -

# # Get Reduced model 
#
# Model reductrion parameters

d = 5
h = 0.1
# collect observations
obs_gap = 1
V_obs = vv[2:d+1,1:obs_gap:end]
vv = []

Psi = mrl.get_Psi_2017(h)

# Get Wiener filter
#@time h_wf = get_wf(V_obs,Psi, M_out = M_out,PI = true)
signal = V_obs
V_obs = []
M_out = 20
n = 3; p = 1500; par = 1500
ty = "bin"
xspec_est = "old"
nfft = 0
rl = true
Preds = false
N_ckms = 3000
PI = false
rtol = 1e-6
Verb = false
tm = now()

paramaters = Dict(
    "M_out" => M_out,
    "n" => n,
    "p" => p,
    "par" => par,
    "ty" => ty,
    "xspec_est" => xspec_est,
    "nfft" => nfft,
    "rl" => rl,
    "Preds" => Preds,
    "N_ckms" => N_ckms,
    "rtol" => rtol,
    "tm" => tm
)

# +
Len = 5000

d, steps = size(signal[:,1:Len])
nu = size(Psi(zeros(d,1)),1)

sig = @view signal[:,2:steps]   # sig is now one a head of signal
steps -= 1                      # this makes steps the length of sig

pred = zeros(ComplexF64, nu, steps)
for n = 1:steps
    pred[:,n] = Psi(@view signal[:,n])
end # pred is now even with signal and therefore one step
# behind sig. I.e. pred[:,n] = Psi(sig[:,n-1])
# which is what we want so as to ensure the reduced
# model can run explicitly.

d, stepsy = size(sig)
nu, stepsx = size(pred)
stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])

sig = Array(transpose(sig))
pred = Array(transpose(pred))

Pred = zeros(ComplexF64,steps-M_out+1,M_out*nu)
for m = 1:M_out
    Pred[:,nu*(M_out-m)+1:nu*(M_out-m+1)] = @view pred[m:(steps + m - M_out),:]
end
# -


PRED = [real(Pred) -imag(Pred); imag(Pred) real(Pred)]
SIG = [real(sig[M_out:steps,:]); imag(sig[M_out:steps,:])];

PRED

h = zeros(ComplexF64,M_out*nu,d)

h_temp = glmnet(PRED,SIG[:,2])


h_temp = glmnet(PRED,SIG[:,1],lambda = reverse(Array(0.00001:0.00001:0.0008)))

h_temp = glmnet(PRED,SIG[:,2],lambda = reverse(Array(0.00001:0.00001:0.0008)))




B = h_temp.betas
i=67
h = B[:,i]
h_c = [complex(h[k],h[k+700]) for k = 1:700]
b = findall(x -> !iszero(x),h)
foreach(x -> println(floor(Int,x/35)+1,"   ", x % 35),b)




X = randn(1000,10) .+ 40
Bet_T = rand(10)*54
y = X*Bet_T

mx = mean_and_std(X)
my = mean_and_std(y)


Xc = zscore(X)
yc = zscore(y)



Out = glmnet(Xc,yc)

Bc = Out.betas

cv = glmnetcv(X,y)

Betc = coef(cv)



for i = 1:d
    cv = glmnetcv(PRED,SIG[:,i])
    locbetamin = argmin(cv.meanloss)
    h_temp = cv.path.betas[:,locbetamin]
    for j = 1:M_out*nu
        h[j,i] = complex(h_temp[j],h_temp[j+M_out*nu])
    end
end

h_wfls = zeros(ComplexF64,d,nu,M_out)
for m = 1:M_out
    h_wfls[:,:,m] = (@view h[(nu*(m-1) + 1):nu*m,:])'
end
h_wfls




# Save Wienerfilter
dat = Dict("dat_h_wf" => h_wf)
Data = merge(paramaters,dat)
# save("../data/KSE_Data/KSE_sol_wienerfilter.jld",Data)

wf_file = server ? "../../../data/KSE_Data/ks_wf_$gen-Len$Len.jld" :
   "C:/Users/jared/Desktop/DDMR/Examples/KSE/Data/ks_wf_$gen-Len$Len.jld"
save(wf_file,Data)
println("Wiener filter saved")

# # Load old Wiener Filter
# LD = load("Data\\KSE_sol_wienerfilter.jld")
# h_wf = LD["dat_h_wf"]
# println("Wiener filter load


# #############################################################################
# #############################################################################


function get_pred(signal, Psi)
    d, steps = size(signal)
    nu = size(Psi(zeros(d,1)),1)

    pred = zeros(ComplexF64, nu, steps)
    for n = 1:steps
        pred[:,n] = Psi(@view signal[:,n])
    end
    pred
end

d, steps = size(signal)
nu = size(Psi(zeros(d,1)),1)

sig = @view signal[:,2:steps]   # sig is now one a head of signal
steps -= 1                      # this makes steps the length of sig

pred = zeros(ComplexF64, nu, steps)
for n = 1:steps
    pred[:,n] = Psi(@view signal[:,n])
end

d, stepsy = size(sig)
nu, stepsx = size(pred)
stepsx == stepsy || print("X and Y are not the same length. Taking min.")
steps = minimum([stepsx stepsy])

Sig = Array(transpose(sig))
Pred = Array(transpose(pred))

PRed = zeros(ComplexF64,steps-M_out+1,M_out*nu)
for m = 1:M_out
    PRed[:,nu*(M_out-m)+1:nu*(M_out-m+1)] = @view Pred[m:(steps + m - M_out),:]
end

using StatsBase

PRED = [real(PRed) -imag(PRed); imag(PRed) real(PRed)]
SIG = [real(Sig[M_out:steps,:]); imag(Sig[M_out:steps,:])]


dt = fit(ZScoreTransform, PRED; dims = 1)

dt












using LinearAlgebra

lambda = 1

A = PRED'*PRED + lambda * I
B = PRED'*sig[M_out:steps,:]

h = A \ B

h_wfrr = zeros(ComplexF64,d,nu,M_out)
for m = 1:M_out
    h_wfrr[:,:,m] = (@view h[(nu*(m-1) + 1):nu*m,:])'
end
h_wfrr

h = PRED \ sig[M_out:steps,:]
h_wfbs = zeros(ComplexF64,d,nu,M_out)
for m = 1:M_out
    h_wfbs[:,:,m] = (@view h[(nu*(m-1) + 1):nu*m,:])'
end
h_wfbs
