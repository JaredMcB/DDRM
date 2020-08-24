
include("modgen_LSDE.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")
include("..\\..\\Tools\\AnalysisToolbox.jl")

using JLD
using PyPlot
using DSP: nextfastfft

A = reshape([-0.5],1,1)
σ = reshape([1],1,1)
Xo = [1]
t_disc = 1000
gap = 10
scheme = "EM"

d = size(A,1)

t_start = 0
t_stop  = 1e6
h       = 1e-2

Δt      = h*gap

@time X = modgen_LSDE(t_start,t_stop,h,
    A = A,
    σ = σ,
    Xo = Xo,
    t_disc = t_disc,
    gap = gap,
    scheme = scheme)

N = size(X,2)

nfft = nextfastfft(N)
X = [X zeros(d,nfft - N)]

X2 = copy(X)
X1 = copy(X)

lags = -10000:10000
A1 = my_crosscor(X1[:],X1[:],lags)
A2 = my_crosscor(X2[:],X2[:],lags)

plot(lags*h*gap,A2)

τ_exp1, τ_int1 = auto_times(X1[:])*Δt
τ_exp2, τ_int2 = auto_times(X2[:])*Δt
N_ef1 = N*Δt/τ_int1
N_ef2 = N*Δt/τ_int2

Psi(x) = x

@time h_wf = get_wf(X,Psi)

h_wf

1 .+ h*A

h_wf_55 = h_wf

h_wf = h_wf_4 # _5, _55

d, N  = size(X)
nu    = size(Psi(X[:,1]),1)
M_out = size(h_wf,3)

X_rm = zeros(d,N); X_rm[:,1:M_out] = X[:,1:M_out]

PSI = zeros(nu,N);
for i = 1:M_out
    PSI[:,i] = Psi(X_rm[:,i])
end

for i = M_out + 1 : N
    X_rm[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out, dims = 2)
    PSI[:,i] = Psi(X_rm[:,i])
end



data = Dict(
        "h_wf_4" => h_wf_4,
        "h_wf_5" => h_wf_5,
        "h_wf_55"=> h_wf_55,
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

save("Data\\LSDE_wfs.jld",data)

data = load("Data\\LSDE_wfs.jld")

h_wf_4 = data["h_wf_4"]

abs.(eigen(h_wf_4[:,:,1]).values)

X_rm = data["X__rm_55_h_4"]

blup = findall(isnan,X_rm[1,:])[1]

findall(x -> x>10^6,X_rm[1,:])[1]

semilogy([abs.(X_rm[:,1:150]') C*bas.^(1:150)])

bas = (abs.(X_rm[1,80]/X_rm[1,60]))^(0.05)
C = abs(X_rm[1,20])/bas^20

bas

Psi(x) = x
X = data["X_55"]
h_wf_55 = data["h_wf_55"]

h_wf = h_wf_55 # _5, _55

d, N  = size(X)
nu    = size(Psi(X[:,1]),1)
M_out = size(h_wf,3)

X_rm = zeros(d,N); X_rm[:,1:M_out] = X[:,1:M_out]

PSI = zeros(nu,N);
for i = 1:M_out
    PSI[:,i] = Psi(X_rm[:,i])
end

for i = M_out + 1 : N
    X_rm[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out, dims = 2)
    PSI[:,i] = Psi(X_rm[:,i])
end


blup = findall(isnan,X_rm[1,:])[1]

bas = (abs.(X_rm[1,80]/X_rm[1,60]))^(0.05)
C = abs(X_rm[1,20])/bas^20

bas

E = eigen(h_wf[:,:,1]).values

abs.(E)
