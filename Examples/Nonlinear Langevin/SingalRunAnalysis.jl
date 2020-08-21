# ##########################################################################
#  Title : SingalRunAnalysis.jl
#  Author: Jared McBride (04/10/2020)
#
#  Description: This script serves as a platform for
#  analysing the preformance of the reduced model.
#
# ##########################################################################

# Packages used
using Statistics
using Plots
pyplot()

include("..\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
include("DataGen.jl")
include("RedModRun.jl")

## Preference parameters
t_start = 0
t_stop = 10^3
steps = 10^5 + 1 # (this includes the t_stop)
discard = 10^4
steps_tot = steps + discard

sig_init = [1.5]
sigma = [1]
sigma_v = sigma
d = 1

dVdx(x) = -x.*(x.^2 .- 1) # Symetric Double well
# dVdx(x) = -(6x.^5 - 8x.^3 + 2x)
Psi(x) = [x; x.^3]

M_out = 20


## Derivitive parameters
steps_tot = steps + discard

nu = size(Psi(0),1)
h = (t_stop - t_start)/(steps)


## Generate signal and predictors
signal = DataGen_Langevin_FE(steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    d = d,
    V_prime = dVdx,
    SM1 = true)
# Notice since SM! is true the output series will include one
# presample. Thus effecting an offset by one in the index.
# Signal is one time step behind, the desired sig.

sig = signal[:,2:end]

plot(sig[1,:])

pred = zeros(nu, steps)
for n = 1:steps
    pred[:,n] = Psi(signal[n])
end

## Compute Wiener Filter
h_wf = vector_wiener_filter_fft(pred, sig, M_out)
plot(real.(h_wf[1,1,:]))

## Run reduced model
sig_rm = redmodrun(
    h_wf, Psi,
    sigma_v = sigma_v,
    sig_init = sig_init,
    steps = steps,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard
    )

win = 1:10500

l = @layout([a b; c d])
plot([sig[1,win], sig[1,:], sig_rm[win], sig_rm], layout=l,
    t=[:line :histogram :line :histogram],
    label = ["sig" "sig" "sig_rm" "sig_rm"],
    leg = true,
    border=:none)

plot(sig_rm)

plot(real.(h_wf[1,2,:]))
