using Plots
pyplot()

include("..\\DataGen.jl")

t_start = 0
t_stop = 10^3
steps = 10^5 + 1 # (this includes the t_stop)
discard = 10^3
steps_tot = steps + discard

sig_init = [1.5]
sigma = [1]
sigma_v = sigma
d = 1

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)

Psi(x) = [x; x.^3]
M_out = 20


steps_FE = 10^6
steps_RK = 10^5

sig_FE = DataGen_Langevin_FE(steps_FE,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx)

sig_RK = DataGen_Langevin_FE(steps_RK,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx)

t_FE = range(t_start,stop = t_stop,length = 10^6)
t_RK = range(t_start,stop = t_stop,length = 10^5)

l = @layout([a b])

p1 = plot(t_FE, sig_FE[1,:]

p1 = plot([sig_FE[1,:] sig_FE[1,:]], sig_RK[1,:] sig_RK[1,:] ],
    layout=l,
    t=[:line :histogram :line :histogram])
