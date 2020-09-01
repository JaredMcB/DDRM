include("..\\Model_Reduction_Dev.jl")
include("..\\AnalysisToolbox_scratch_ckms.jl")

using LinearAlgebra
using PyPlot
pygui(true)

function visual_test_ckms(P,l,nfft;diff = false)
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

    if diff
        for i = 1:d
            for j = i:d
                plot(abs.(S[i,j,:] - S_l[i,j,:]), label = "S $i $j")
            end
        end
    else
        for i = 1:d
            for j = i:d
                plot(S[i,j,:], label = "S $i $j")
                plot(S_l[i,j,:], label = "S_l $i $j")
            end
        end
    end
    legend()
end

d = 3; m = 2
P = zeros(3,3,m)
b = randn(3,3)
P[:,:,i] = b*b' + 1e-10*I
for i = 2 : m
    b = randn(3,3)
    P[:,:,i] = b
end

# T = [isposdef(P[:,:,i]) for i = 1:m]

nfft = 10^4
l, Err = spectfact_matrix_CKMS_SC(P);
Err
semilogy(Err)
visual_test_ckms(P,l,nfft,diff = true)

 
