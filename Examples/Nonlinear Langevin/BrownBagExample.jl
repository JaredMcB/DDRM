include("C:\\Users\\jared\\Desktop\\Github Repos\\DDMR\\Tools\\Wiener Filtering\\Matrix Wiener Filter\\wiener_filter_Matrix_fft.jl")
# include("C:\\Users\\jared\\Desktop\\Github Repos\\DDMR\\Tools\\Model_Reduction_Dev.jl")

include("DataGen.jl")
include("RedModRun.jl")


t_start = 0
t_stop  = 150000

sig_init = [1.5]
sigma = [.35]
sigma_v = sigma
d = 1

gap = 10^2

dt = 10^-3
Δt = gap*dt

Time = range(t_start,t_stop, step = dt)
N_grid = range(t_start,t_stop, step = Δt)

discard_N = 10^5
discard_T = gap*discard_N

T = length(Time)
N = length(N_grid)



nu = 2
M_out = 20
Nen = 100
scheme = "EM"

dW = randn(1,T + discard_T)
W = cumsum(dW,dims = 2)
ΔW = zeros(d,N + discard_N)
for i = 1:(N + discard_N - 1)
    ΔW[:,i+1] = sqrt(dt/Δt)*(W[:,gap*i] - W[:,gap*(i-1)+1])
end

@time signal_T = DataGen_DWOL(T,
    scheme = scheme,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard_T,
    sig_init = sig_init,
    sigma = sigma,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    e = dW)

signal_N = signal_T[1,1:gap:end]


parameters = Dict(
    "t_start" => t_start,
    "t_stop" => t_stop,
    "sig_init" => t_stop,
    "sigma" => sigma,
    "sigma_v" => sigma_v,
    "d" => d,
    "gap" => gap,
    "dt" => dt,
    "Δt" => Δt,
    "discard_N" => discard_N,
    "discard_T" => discard_T,
    "Psi" => "Double Well potential",
    "nu" => nu,
    "M_out" => M_out,
    "Nen" => Nen,
    "scheme" => scheme
)


dat = Dict(
    "dat_signal_N" => signal_N)



save("Examples\\Nonlinear Langevin\\data\\full_model_run.jld",
    merge(parameters,dat))

data_dict = load("Examples\\Nonlinear Langevin\\data\\full_model_run.jld")
signal_N = data_dict["dat_signal_N"]

plot(N_grid,signal_N,
    title = "Double-welled Langevin Run",
    xaxis = "time",
    yaxis = "position",
    leg = :none)


Psi1(x) = [x ; x.^3]

signal_N_vec = reshape(signal_N,1,:)
h_wf = get_wf(signal_N, Psi1)
h_wf_new = h_wf
h_wf_old = zeros(1,2,20)
for i = 1:M_out
    h_wf_old[:,:,i] = h_wf_new[:,:,i]'
end
sig_rm = redmodrun(
    h_wf_old, Psi1;
    sigma_v,
    sig_init,
    steps = N,
    t_start,
    t_stop,
    discard = discard_N
    )
