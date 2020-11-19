# using PyPlot
using JLD


include("DataGen.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")

function run_get_wf(;
    Psi      = x -> [x; x.^3],
    # Model run parameters
    sigma    = [.3],
    V_prime  = x -> -x.*(x.^2 .- 1),
    sig_init = [1.5],
    # Numerical estimate parameters
    scheme   = "FE",
    steps    = 10^7, # Number of time steps (not including those discarded)
    h        = .1,
    discard  = steps ,# Number of time steps discarded
    gap      = 1,    # 1 + the number of time steps between observations
    # WF parameters
    Parms,
    rl = true,
    Nens = 100)

    d = size(sig_init,1)
    nu = size(Psi(sig_init),1)
    P = length(Parms)

    h_wfs_ens = rl ? zeros(P,Nens,d,nu,M_out) :
                     complex(zeros(P,Nens,d,nu,M_out))
    for n = 1 : Nens
        X = @time DataGen_DWOL(;
            #SDE parameters
            sigma, V_prime, sig_init,
            # Numerical estimate parameters
            scheme, steps, h, discard, gap)

        for i = 1:P
            h_wfs_ens[i,n,:,:,:] = get_wf(X, Psi;
                info = false, rl,
                xspec_est = Parms[i][1],
                par       = Parms[i][2],
                nfft      = Parms[i][3],
                n         = Parms[i][4],
                p         = Parms[i][5],
                ty        = Parms[i][6],
                M_out     = Parms[i][7]);
        end
    end
    h_wfs_ens
end

#SDE parameters
sigma    = [.3]
V_prime  = x -> -x.*(x.^2 .- 1)
sig_init = [1.5]
# Numerical estimate parameters
scheme   = "FE"
steps    = 10^7  # Number of time steps (not including those discarded)
h        = .01
discard  = steps # Number of time steps discarded
gap      = 100     # 1 + the number of time steps between observations

Psi(x) = [x; x.^3]

# Model reduction Parameters
M_out = 20
ty = "bin"
### Varing parameters
###       xspect_est , par    , nfft    , n    , p
#
Parms = [["DM"       , 5000  , 2^17    , 2    , 5, ty, M_out],
         ["SP"       , 5000  , 2^17    , 2    , 5, ty, M_out]]

Nens = 1

h_wfs_ens = run_get_wf(;Psi,
    sigma, V_prime, sig_init,
    scheme, steps, h, discard, gap,
    Parms,
    Nens)

data = Dict("sigma"         => sigma,
            "V_prime_str"   => "x -> -x.*(x.^2 .- 1)",
            "sig_init"      => sig_init,
            # Numerical estimate parameters
            "scheme"        => scheme,
            "steps"         => steps,
            "h"             => h,
            "discard"       => discard,
            "gap"           => gap,
            "h_wfs_ens"      => h_wfs_ens)

# save("../../../data/DWOL_Data/data_11_18_2020_2.jld",data)

v_sp = var(h_wfs_ens[2,:,1,1,1])
m_sp = mean(h_wfs_ens[2,:,1,1,1])

v_dm = var(h_wfs_ens[1,:,1,1,1])
m_dm = mean(h_wfs_ens[1,:,1,1,1])

println("===================================================")
println("mean and varieance of DM WF: $m_dm and $v_dm")
println("mean and varieance of Sp WF: $m_sp and $v_sp")
