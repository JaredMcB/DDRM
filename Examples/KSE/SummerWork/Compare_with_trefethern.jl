using PyPlot
using Statistics: mean, var

at = include("../../../Tools/AnalysisToolbox.jl")
kse = include("../Model_KSE.jl")
ksed = include("../Model_KSE_Dev.jl")

T        = 1000 # Length (in seconds) of time of run
T_disc   = T ÷ 2    # Length (in seconds) of time discarded
P        = 21.55  # Period
N        = 96  # Number of fourier modes used
h        = 0.001  # Timestep
g        = x -> cos(x/16)*(1 + sin.(x/16))
obs_gap  = 100

Δt = h*obs_gap
uu, vv, tt    =  @time kse.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)
uud, vvd, ttd =  @time ksed.my_KSE_solver(T; T_disc, P, N, h, g, n_gap = obs_gap)

sum(abs.(uu - uud).^2)

t_start = 0
t_stop = 150
ind_start = floor(Int,t_start/Δt)+1
ind_stop =floor(Int,t_stop/Δt)

H1 = imshow(reverse(uu[:,1:40],dims = 2)', extent=[0,P,0,10], aspect="auto")
figure()
H2 = imshow(reverse(uud[:,1:40],dims = 2)', extent=[0,P,0,10], aspect="auto")


mean(uu, dims = 2)

u = rand(128)
using FFTW

x = P*(1:N)/N
u = g.(x)
v = fft(u)

function NonLin(v)
    v_pad = [v; zeros(N)]
    nv = fft(real(ifft(v_pad)).^2)
    nv[1:N]
end

vd = ifft(u)

function NonLind(v)
    v_pad = [v; zeros(N)]
    nv = ifft(real(fft(v_pad)).^2)
    nv[1:N]
end

Nv = NonLin(v)
Nvd = NonLind(vd)


EE = (abs.(vv[:,40:end]).^2)
EEd =(abs.(vvd[:,40:end]).^2)

H1 = imshow(reverse(EE,dims = 2)', extent=[0,P,0,10], aspect="auto")
figure()
H2 = imshow(reverse(EEd,dims = 2)', extent=[0,P,0,10], aspect="auto")

E = mean(EE,dims = 2)
Ed = mean(EEd,dims = 2)

semilogy([E 16000*Ed])

v2 =vv[3,:]
vd2 =vvd[3,:]
lags = 1:300
A = complex(zeros(5,length(lags)))
Ad = complex(zeros(5,length(lags)))

for i = 1:5
    A[i,:] = at.my_autocor(vv[i+1,:],lags)
    Ad[i,:] = at.my_autocor(vvd[i+1,:],lags)
end
for i = 1:5
    figure()
    plot(real(A[i,:]), label = "old $i")
    plot(real(Ad[i,:]),":",label = "new $i")
    legend()
end

EE = mean(abs2.(vv),dims = 2)
EEd = mean(abs2.(vvd),dims = 2)
for i = 2:6
    plot(EE[i,:])
end

semilogy([EE EEd])
