include("../AnalysisToolbox.jl")

X = randn(1000000) + im*randn(1000000)
lags = -100:10
T = my_crosscov(X,X,-100:10)
x = y = X
L = min(length(x),length(y))
m = length(lags)

zx = x .- mean(x)
zy = y .- mean(y)

C = zeros(Complex,m)
for k = 1:m
    l = lags[k]
    C[k] = ( l >= 0 ? dot(zx[1+l : L],zy[1 : L-l]) : dot(zx[1 : L+l],zy[1-l : L]))/L
end
C
dot(x,y)
l = 1


using LinearAlgebra

M = 2*10^4
dot(randn(M) + im*randn(M),randn(M) + im*randn(M))
