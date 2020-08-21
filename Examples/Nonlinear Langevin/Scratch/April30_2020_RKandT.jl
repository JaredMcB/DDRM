
using Plots
ENV["MPLBACKEND"]="qt5agg" # Change backend to allow output plots
pyplot()


include("..\\DataGen.jl")
include("..\\RedModRun.jl")
include("..\\..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")


t_start = 0
t_stop = 10^4
steps = 2*10^4 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [.5]
sigma_v = sigma
d = 1

Nen = 500

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)
Psi(x) = [x; x.^3]

M_out = 20

# Derivitive parameters
steps_tot = steps + discard
nu = size(Psi(0),1)
dt = (t_stop - t_start)/(steps)
tim = range(t_start,t_stop,length = steps)
e = randn(d,10^4+10^3)

##True run Euler-Maruyamma dt = 10^-4
t_start = 0
t_stop = 100
dt = 2^-10
steps = Int64(round( (t_stop - t_start)/dt )) + 1
sig_init = [1.5]
sigma = [.5]
tim = range(t_start,t_stop,length = steps)

signal, e = DataGen_DWOL(steps,
    scheme = "FE",
    t_start = t_start,
    t_stop = t_stop,
    discard = 0,
    sig_init = sig_init,
    sigma = sigma,
    d = 1,
    V_prime = dVdx,
    SM1 = false,
    Obs_noise = true)

P = timeseries_plot(tim,signal[1,:],title = "Truth")

## test the other schemes with deltat = 10^-2

deltat = 2^-2

steps_N = Int64(round( (t_stop - t_start)/deltat )) + 1
tim_N = range(t_start,t_stop,length = steps_N)

scale = Int(round(deltat/dt))
e_N = zeros(1,steps_N)
for i = 1 : steps_N-1
    e_N[:,i] = sum(e[:,scale*(i-1)+1:scale*i],dims = 2)
end
signal_et = DataGen_DWOL(steps_N,
    scheme = "FE",
    t_start = t_start,
    t_stop = t_stop,
    discard = 0,
    sig_init = sig_init,
    sigma = sigma,
    d = 1,
    V_prime = dVdx,
    SM1 = false,
    Obs_noise = false,
    e = e_N)

plot(tim[1:1000],signal[1,1:1000])
plot!(tim_N[1:5],signal_et[1,1:5])











P_et = timeseries_plot(tim,signal_et[1,:],title = "Expl trap")
P_em = timeseries_plot(tim,signal_em[1,:],title = "Expl Midp")
plot(P_et,P_em,layout = (2,1))
wind = 1:200
plot(tim[wind],[signal_fe[1,wind] signal_t2[1,wind]])

Psi(x) = [x; x.^3; x.^5]
h_wf = get_wf(signal_fe[:,1:5:end],Psi)

e_hat = zeros(1,2*10^2 + 2*10^3)
for i = 1 : 2*10^2 + 2*10^3
    e_hat[1,i] = sum(e[:,5*(i-1)+1 : 5*i])
end
sig_hat = redmodrun(
    h_wf, Psi,
    sigma_v = [1],
    sig_init = [1.5],
    steps = 2*10^3,
    t_start = 0,
    t_stop = 10^3,
    discard = 2*10^2,
    obs_noise = e_hat
    )

signal_t2 = DataGen_DWOL(2*10^3,
    scheme = "T2",
    t_start = 0,
    t_stop = 10^3,
    discard = 2*10^2,
    sig_init = sig_init,
    sigma = sigma,
    d = 1,
    V_prime = dVdx,
    SM1 = false,
    Obs_noise = false,
    e = e_hat)
tim_hat = range(0,10^3,length = 2*10^3)

plot(tim[1:200],signal_fe[1,1:200])
plot!(tim_hat[1:40],[sig_hat[1,1:40] signal_t2[1,1:40]])
