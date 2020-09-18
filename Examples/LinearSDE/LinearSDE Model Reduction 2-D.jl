using PyPlot
using Statistics
using Distributions

include("..\\..\\Tools\\Model_Reduction_Dev.jl")
include("modgen_LSDE.jl")

A       = - .5*[1 1; 1 1.1]
σ       = I + zeros(2,2)
Xo      = [1; 10]
t_disc  = 1000
gap     = 10
d       = size(A,1)
t_start = 0
t_stop  = 1e5
h       = 1e-2
Δt      = h*gap
M_out   = 100

F = I + h*A
eigen(F)
cond(F)

X = modgen_LSDE(t_start,t_stop,h;
    A, σ, Xo, t_disc, gap)
X_full = X
X = X[1:1,:]
N       = size(X,2)
nfft    = nextfastfft(N)
d_rm = 1
X = [X zeros(d_rm,nfft - N)]

Psi(x) = x

@time h_wf = get_wf(X[1:1,:],Psi; M_out)

X = X[:,1:N]
M_out = 20
nu    = size(Psi(X[:,1]),1)

X_rm = zeros(d_rm,N); X_rm[:,1:M_out] = X[:,1:M_out]

PSI = zeros(nu,N);
for i = 1:M_out
    PSI[:,i] = Psi(X_rm[:,i])
end

for i = M_out + 1 : N
    X_rm[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k] for k = 1:M_out,
                dims = 2) + sqrt(Δt)*rand(D,d_rm)
    PSI[:,i] = Psi(X_rm[:,i])
    isnan(X_rm[1,i]) && break
end
maximum(X_rm)

plot(X[1,:],X[2,:])
plot(X_rm[1,:],X_rm[2,:])
size(X_rm)

plot([X[1,:] X_rm[1,:]])

onesteperrors = zeros(N-M_out)
for i = M_out + 1 : N
    onesteperrors[i-M_out] = (X[:,i] - sum(h_wf[:,:,k]*X[:,i-k] for k = 1:M_out, dims = 2))[1]
end

onesteperrors

plot(onesteperrors)

emp_pdf(onesteperrors)
μ = mean(onesteperrors)
v = var(onesteperrors)

D = Normal(μ,sqrt(v))

xx = -7.5:1/100:7.5
pdf_D(x) = pdf(D,x)

plot(xx,pdf_D.(xx))
