# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Julia 1.5.2
#     language: julia
#     name: julia-1.5
# ---

# # Varying `gap` with on long timeseries

using PyPlot
using Random
using JLD
include("DataGen.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")

# SDE Parameters
sigma    = [.3]
V_prime  = x -> -x.*(x.^2 .- 1)
sig_init = [1.5]
# NSDE Parameters
scheme   = "FE"
steps    = 10^6
h        = .01
discard  = steps
gap      = 1
seed     = 2016

## Get full model run
Random.seed!(seed)
X = @time DataGen_DWOL(;
    sigma, V_prime, sig_init,
    scheme, steps, h, discard, gap)

## Check the run
plot(X[1,1:end÷8000:end],ms=1,".")

# +
## Save the data
data = Dict("sigma"         => sigma,
            "V_prime_str"   => "x -> -x.*(x.^2 .- 1)",
            "sig_init"      => sig_init,
            # Numerical estimate parameters
            "scheme"        => scheme,
            "steps"         => steps,
            "h"             => h,
            "discard"       => discard,
            "gap"           => gap,
            "X"             => X)
            
save("../../../data/DWOL_Data/data_11_19_2020_1X.jld",data)
# -

## Load Data
X = load("../../../data/DWOL_Data/data_11_19_2020_1X.jld","X")

# ## Collect WF information

## Put in Psi functions
Psi(x) = [x; x.^3]

## Model reduction Parameters
M_out = 20
ty = "bin"

## Varing parameters
##       xspect_est , par    , nfft    , n    , p
Parms = [["DM"       , 5000  , 2^17    , 2    , 5],
         ["SP"       , 5000  , 2^17    , 2    , 5]]
Gaps = [100,10,1]

h_wf_packs  = []
for g in Gaps
    for i = 1 : length(Parms)
        Out = get_wf(X[:,1:g:end], Psi;
            M_out, ty, info = true,
            xspec_est = Parms[i][1],
            par       = Parms[i][2],
            nfft      = Parms[i][3],
            n         = Parms[i][4],
            p         = Parms[i][5]);
        append!(h_wf_packs,Out)
    end
end




