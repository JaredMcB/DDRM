mr = include("../WFMR.jl")
at = include("../AnalysisToolbox.jl")


using PyPlot
using DSP
using Polynomials


 ## 
steps = 10^4
x = randn(steps)

Zeros = [.5 -.6]

h = coeffs(prod([Polynomial(1) [Polynomial([1, -z]) for z in Zeros]]))

y = fftfilt(h,x)

# y1 = zeros(steps)
# for i = length(h) : steps
#     y1[i] = sum(h[j]*x[i-j+1] for j = 1:length(h))
# end
#
# my_diff(x,y) = maximum(abs.(x - y))
#
# my_diff(y[3:end],y1[3:end])

@time h_whf = mr.get_whf(y)

h_whf = real(reshape(h_whf,:))

x_rec = fftfilt(h_whf,y)

x
x_rec

hist(x)
figure()
hist(x_rec)

lags = -1000:1000
A_rec = at.my_autocov(x_rec,lags)
A = at.my_autocov(x,lags)

plot(lags,[A A_rec])
