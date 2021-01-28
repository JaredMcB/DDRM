using PyPlot
using FFTW

ks = include("ks.jl")
at = include("../Tools/AnalysisToolbox.jl")


P = 32π
steps = 100
N = 128

# Trefethen
vv = ks.run_ks(steps,1.5;
             nsubsteps = 6,
             verbose = false, L = P)

# get solution
uu = real(N*ifft([zeros(1, steps+1); vv; reverse(conj(vv),dims = 1)],1))
H1 = imshow(reverse(reverse(uu',dims = 1),dims = 2), extent=[0,P,0,150], aspect="auto")


VV = [vv_lin, vv_myNA, vv_myA]
UU = [uu_my,uu_lin, uu]

sty = ["-","--",":"]
labels = ["vv_lin", "vv_myNA", "vv_myA"]
for i = 1:5
    for j = 1 : length(VV)
        subplot(5,1,i)
        plot(1.5*(1:300),real(VV[j][i,1:300]),linestyle = sty[j])
    end
    title("model $i")
end
suptitle("Comparison of L's nonlinear part (solid) with M's nonlinear part (dealiased, dashed;aliased, dotted) in L's ODE solver")






v = ks.kt_init()

vc = real(N*ifft([0; v[:];0;reverse(conj(v[:]))]))

vi = ks.kt_init(real = true)

maximum(abs.(vi - vc))

EE = abs2.(vv)

semilogy(mean(EE,dims = 2))

# 2017
steps = 10000
vv = ks.run_ks(ks.fei_init(), steps,.1;
             nsubsteps = 100,
             verbose = false)

A = at.my_autocor(vv[2,:],-150:150)

plot(A)
