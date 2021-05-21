mr = include("../WFMR.jl")
at = include("../AnalysisToolbox.jl")


using PyPlot
using DSP
using Polynomials


 ##
steps = 10^4
x = randn(steps)

Zeros = [.5 -.6 .9]

P = prod([Polynomial(1) [Polynomial([1, -z]) for z in Zeros]])
h = coeffs(P)

y = fftfilt(h,x)

y1 = zeros(steps)
for i = length(h) : steps
    y1[i] = sum(h[j]*x[i-j+1] for j = 1:length(h))
end

my_diff(x,y) = maximum(abs.(x - y))

my_diff(y[4:end],y1[4:end])

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


nfft = 10^4
Theta = 2π*(0:nfft-1)/nfft

s_ana(z) = P(z)*P((z^-1)')'
S_ana = real(map(s_ana,exp.(im*Theta)))

plot(Theta,S_ana)

S_num = at.z_crossspect_dm(y,y, L = 55, Nex = nfft)

plot(Theta,S_num);plot(0,0)
plot(Theta,S_num_rec);plot(0,0)
S_num = at.z_crossspect_dm(x,x, L = 500, Nex = nfft)
S_num_rec = at.z_crossspect_dm(x_rec,x_rec, L = 500, Nex = nfft)
S_num
