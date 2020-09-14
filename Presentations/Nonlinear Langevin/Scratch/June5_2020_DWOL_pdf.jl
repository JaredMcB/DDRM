using JLD
using Printf
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
# ENV["MPLBACKEND"]="agg" # Change backend to allow output plots
pyplot()

cd("C:\\Users\\JaredMcBride"*
    "\\Desktop\\Git Stuff\\ModelReduction\\"*
    "Spectral Methods\\Model Reduction Julia Scripts\\Scratch")

include("..\\DataGen.jl")
include("..\\DWOL_eqidist_sampler.jl")

## Compute emperical pdf
data_dict = load("data\\full_model_run_pdf.jld")
pdf_N = data_dict["dat_pdf_N"]
b_mpts = data_dict["dat_b_mpts"]

function p(z)
    i = argmin(abs.(b_mpts .- z))
    if (i == 1) | (i == length(b_mpts))
        x = [1; 2; 3]
        y = zeros(3)
    else
        x = b_mpts[i-1:i+1]
        y = pdf_N[i-1:i+1]
    end
    a = [x.^2 x ones(3)] \ y
    P = a[1]*z^2 + a[2]*z + a[3]
end

# Numerically compute exact pdf by Fokker-Planck equations
σ = data_dict["sigma"][1]

function RHS!(dU,U,μ,t)
    dU[1] = U[2]
    dU[2] = -1/μ*(3t^2 - 1)*U[1] - 1/μ*(t^3 -t)*U[2]
end

U0 = [p(0);0]
tspan = (0.0,2.0)
μ = σ^2/2
prob = ODEProblem(RHS!,U0,tspan,μ)
sol = solve(prob)

p_num(x) = x>=0 ? sol(x)[1] : sol(-x)[1]

## Compare results

x = -2:0.01:2
plot(x,p.(x))
plot!(x,p_num.(x),line = :dash)
