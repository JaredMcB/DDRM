include("..\\Model_Reduction_Dev.jl")
include("AnalysisToolbox_scratch_ckms.jl")

using PyPlot
pygui(true)
### Example 1###
P = zeros(1,1,2)
P[1,1,1] = 10
P[1,1,2] = 3

L = spectfact_matrix_CKMS(P)

l, Err = spectfact_matrix_CKMS_SC(P);
visual_test_ckms(P,l,nfft)


### Example 2 ###
P = zeros(2,2,2)
P[:,:,1] = [10 0 ; 0 84]
P[:,:,2] = [3 0; 0 38]

l = spectfact_matrix_CKMS(P)

ρ = (-42 + sqrt(42^2 - 38^2))/38

(sqrt(-38/ρ),-ρ*sqrt(-38/ρ))


### Example 3 ###
P = zeros(2,2,10)
P[:,:,1] = [100 30 ; 12 84]
P[:,:,2] = [3 1; -1 38]
P[:,:,3] = [3 1; 0 9]
P[:,:,4] = [3 0; 1 9]


l, Err = spectfact_matrix_CKMS_SC(P);
l

nfft = 10^3
visual_test_ckms(P,l,nfft)



function visual_test_ckms(P,l,nfft)
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
            plot(S[i,j,:], label = "S $i $j")
            plot(S_l[i,j,:], label = "S_l $i $j")
        end
    end
    legend()
end


### Example 4 from JLE ###
P = zeros(2,2,2)
P[:,:,1] = [6 22; 22 84]
P[:,:,2] = [2 7; 11 38]

l1 = spectfact_matrix_CKMS_SC(P);

l = zeros(2,2,2)
l[:,:,1] = [2 1; 7 3]
l[:,:,2] = [1 0; 5 1]

l2 = l
l1 = l1[1]
ll = size(l1,3)

S1_fun_minus(z) = sum(l1[:,:,i]*z^(-i+1) for i = 1:ll)
S2_fun_minus(z) = sum(l2[:,:,i]*z^(-i+1) for i = 1:ll)

res(z) = S1_fun_minus(z)*S1_fun_minus(z^(-1))' -
            S2_fun_minus(z)*S2_fun_minus(z^(-1))'

d= 2; nfft = 10^3
Res = complex(zeros(d,d,nfft))
for i = 1:nfft
    Res[:,:,i] = res(exp(im*2π*i/nfft))
end
M = zeros(size(Res[:,:,1]))
for i = 1:d
  for j = i:d
    plot(abs.(Res[i,j,:]), label = "($i,$j)")
    M[i,j] = maximum(abs.(Res[i,j,:]))
  end
end
legend()
M1 = M
M2 = M

M2


M1
