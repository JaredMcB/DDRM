include("modgen_LSDE.jl")
include("..\\..\\Tools\\Model_Reduction_Dev.jl")

using JLD
using PyPlot
using DSP: nextfastfft

A = reshape([-0.5],1,1)
σ = reshape([1],1,1)
Xo = [1]
t_disc = 1000
gap = 1
scheme = "EM"

t_start = 0
t_stop  = 1e6
h       = 1e-2

Δt      = h*gap
M_out   = 20

X = modgen_LSDE(t_start,t_stop,h,
    A = A,
    σ = σ,
    Xo = Xo,
    t_disc = t_disc,
    gap = gap,
    scheme = scheme)

d, N = size(X)

nfft = nextfastfft(N)
X = [X zeros(d,nfft-N)]


τ_exp, τ_int    = auto_times(X[:])*Δt
N_eff           = N*Δt/τ_int

Psi(x) = x

@time h_wf = get_wf(X, Psi, par = 2000)

(τ_exp, N_eff)

1 .+ h*A

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
M_out = size(h_wf,3)

X_rmn = zeros(d,N); X_rmn[:,1:M_out] = X[:,1:M_out]

PSI = zeros(nu,N);
for i = 1:M_out
    PSI[:,i] = Psi(X_rmn[:,i])
end

for i = M_out + 1 : N
    X_rmn[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out, dims = 2) + sqrt(h)*σ*randn(d)
    PSI[:,i] = Psi(X_rmn[:,i])
end



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

save("Data\\LSDE_wfs.jld",data)

data = load("Data\\LSDE_wfs.jld")

blup = findall(isnan,X_rm[1,:])[1]

findall(x -> x>10^1,X_rm[1,:])[1]

plot([X_rmn[:] X[:]])

hist(X[:],100,alpha = .5)
hist(X_rmn[:],100,alpha = .5);

lags = 0:1000

A_rm = my_autocor(X_rm[:],lags)
A    = my_autocor(X[:],lags)

plot(lags*h,[A A_rm])

z_spect = z_spect_scalar(X[:], n = 3, p=100, ty = "ave")
z_spect_rm = z_spect_scalar(X_rm[:], n = 3, p=100, ty = "ave")

Θ = 2*π*(1:1000:nfft)/nfft
Z = exp.(im*\Theta)

a = 1 + h*A[1,1]
z_spect_ana_fun(z) = h*\sigma/( (1-a*z^(-1))*(1-a*z) )
z_spect_ana = real(z_spect_ana_fun.(Z))

semilogy(Θ,[z_spect[1:1000:nfft] z_spect_rm[1:1000:nfft] z_spect_ana])
