
using StatsBase
using Statistics
using DSP
using FFTW
using Distributions
using JLD

function DataGen_DWOL(;
    #SDE parameters
    sigma    = [1],
    V_prime  = x -> -x.*(x.^2 .- 1),
    sig_init = [1.5],
    # Numerical estimate parameters
    scheme   = "FE",
    steps    = 10^6, # Nomber of time steps (not including those discarded)
    h        = .01,
    discard  = steps, # Number of time steps discarded
    gap      = 100, #1 + the number of time steps between observations
    ObsNoise = false, # if true keeps track of an returns noise e
    e        = ObsNoise ? randn(size(sigma,1),steps + discard) : 0
    )

    d = size(sigma,1)

    if discard == 0
        sig_init = [DWOL_dist_samp(1,σ = sigma)]
    end

    size(sig_init,1) == size(sigma,1) || print("Dimension of initial value "*
                                               "and sigma do not agree")

    steps_tot = steps + discard

    # So the multiplication is defined when d=1
    d == 1 && (sigma = reshape(sigma,1,1))

    # Here we genereate the signal process.
    signal = zeros(d,ceil(Int,steps/gap))
    tmp = sig_init
    if scheme == "FE"
        for n = 1 : steps_tot-1
            tmp = tmp + h*V_prime(tmp) +
                        sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "T2"
        for n = 1 : steps_tot-1
            tmp = tmp .+ h*(-tmp.^3 .+ tmp) .+
                            h^2/2*( tmp.*(tmp.^2 .- 1).*(tmp.^2 .- 1) ) .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "EM"
        for n = 1 : steps_tot-1
            k1 = tmp .+ h/2*V_prime(tmp)
            tmp = tmp .+ h*V_prime(k1) .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "ET"
        for n = 1 : steps_tot-1
            k0 = V_prime(tmp)
            k1 = tmp .+ h*k0
            tmp = tmp .+ h/2*V_prime(k1) .+ h/2*K0 .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "RK4"
        for n = 1 : steps_tot-1
            k1 = tmp
            k2 = tmp .+ h/2*V_prime(k1)
            k3 = tmp .+ h/2*V_prime(k2)
            k4 = tmp .+ h*V_prime(k3)
            tmp = tmp .+ h/6*( V_prime(k1) .+ 2*V_prime(k1) .+
                            2*V_prime(k1) .+ V_prime(k1) ) .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    end

    result = ObsNoise ? [signal, e] : signal
end

function DWOL_dist_samp(
    N::Int64 = 1; # N is the number of samples you want
    σ = [.35] # The σ for the noise of the the DWOL
    )
    # sampling
    μ_pos = 1
    μ_neg = -1
    v = 1/(sqrt(2π)*1.05)

    D_neg = Normal(μ_pos,v)
    D_pos = Normal(μ_neg,v)

    q(x) = exp(-(x-μ_pos)^2/(2v^2))/(2v*sqrt(2\pi)) +
            exp(-(x-μ_neg)^2/(2v^2))/(2v*sqrt(2\pi))

    μ = σ[1]^2/2
    p(x) =  exp(-(x^2-1)^2/μ)
    c = 1.05*p(0)/q(0)
    cq(x) = c*q(x)
    Z = zeros(N)
    for n = 1:N
        z = (rand() < .5 ? rand(D_neg) : rand(D_pos))
        while rand()*cq(z) > p(z)
            z = (rand() < .5 ? rand(D_neg) : rand(D_pos))
        end
        Z[n] = z
    end
    N == 1 ? Z[1] : Z
end


function get_pred(sig, Psi)
    d, steps = size(sig)
    nu = size(Psi(zeros(d,1)),1)

    pred = zeros(nu, steps)
    for n = 1:steps
        pred[:,n] = Psi(sig[:,n])
    end
    pred
end

function Run_and_get_WF(
    Psi;
    simulator = DataGen_Langevin_FE,
    steps = 10^5 + 1,
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [1],
    d = 1,
    V_prime = (x -> -x.*(x.^2 .- 1)),
    Nen = 100,
    M_out = 20,
    rl = true # if ture output made REAL
    )

    ##

    nu = size(Psi(zeros(d,1)),1)
    h = (t_stop - t_start)/(steps)

    # The number of series in the ensamble
    h_wf_ens = complex(zeros(d,nu,M_out,Nen))

    for k = 1:Nen
        if simulator ==  DataGen_Lorenz_63
            signal = simulator(steps,
                t_start = t_start,
                t_stop = t_stop,
                discard = discard,
                sig_init = sig_init,
                sigma = sigma,
                SM1 = true)
        else
            signal = simulator(steps,
                t_start = t_start,
                t_stop = t_stop,
                discard = discard,
                sig_init = sig_init,
                sigma = sigma,
                d = d,
                V_prime = V_prime,
                SM1 = true)
        end
        # Notice since SM! is true the output series will include one
        # presample. Thus effecting an offset by one in the index.
        # Signal is one time step behind, the desired sig.

        sig = signal[:,2:end]

        pred = zeros(nu, steps)
        for n = 1:steps
            pred[:,n] = Psi(signal[:,n])
        end
        h_wf = vector_wiener_filter_fft(pred, sig, M_out)
        h_wf_ens[:,:,:,k] = h_wf
    end
    rl == true ? real.(h_wf_ens) : h_wf_ens
end


function Run_and_get_WF_DWOL(
    Psi;
    scheme = "FE",
    steps = 10^5 + 1,
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [1],
    d = 1,
    V_prime = (x -> -x.*(x.^2 .- 1)),
    Nen = 100,
    M_out = 20,
    rl = true # if ture output made REAL
    )

    ##

    nu = size(Psi(zeros(d,1)),1)
    h = (t_stop - t_start)/(steps)

    # The number of series in the ensamble
    h_wf_ens = complex(zeros(d,nu,M_out,Nen))
    for k = 1:Nen
        signal = DataGen_DWOL(steps,
            scheme = scheme,
            t_start = t_start,
            t_stop = t_stop,
            discard = discard,
            sig_init = sig_init,
            sigma = sigma,
            V_prime = x -> -x.*(x.^2 .- 1),
            SM1 = true,
            Obs_noise = false)

        # Notice since SM! is true the output series will include one
        # presample. Thus effecting an offset by one in the index.
        # Signal is one time step behind, the desired sig.

        h_wf_ens[:,:,:,k] = get_wf(signal,Psi)
    end
    rl == true ? real.(h_wf_ens) : h_wf_ens
end

function VarTrunc(h_wf)
    d, nu, M_out = size(h_wf)
    h_trunc = zeros(d, nu, M_out,M_out)
    for i = 1 : M_out
        h_trunc[:,:,1:i,i] = h_wf[:,:,1:i]
    end
    h_trunc
end

function DataGen_Lorenz_63(
    steps = 10^5 + 1;
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = randn(3,1),
    sigma = I,
    SM1 = false,
    Obs_noise = false
    )

    d = 3

    steps_tot = steps + discard
    h = (t_stop - t_start)/(steps)
    t = range(t_start,stop = t_stop,length = steps)

    sigm = 10;
    bet = 8/3;
    rho = 28;

    F(x) = [sigm*(x[2] - x[1]);
        x[1]*(rho - x[3]) - x[2];
        x[1]*x[2] - bet*x[3]]

    # Here we genereate the signal process.
    e = randn(d,steps_tot)
    signal = zeros(d,steps_tot)
    signal[:,1] = sig_init
    for n = 2 : steps_tot
        signal[:,n] = signal[:,n-1] .+ h*F(signal[:,n-1]) + sqrt(h)*sigma*e[:,n]
    end

    # Here we discard the transient part
    sig = signal[:,discard+1:steps_tot]

    # sig_minus_1 is returned as it is needed in produceing the predictor series
    # which is offset by one.
    result = SM1 ? [signal[:,discard] sig] : sig
    result = Obs_noise ? [result, e] : result
end

function Autocov(steps::Int64;
    lag::UnitRange{Int64} = 0:999,
    bN = 100,
    scheme = "FE",
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    V_prime = x -> -x.*(x.^2 .- 1),
    vari = false)

    A = zeros(length(lag),bN)
    for n = 1:bN
        signal = DataGen_DWOL(steps,
            scheme = "FE",
            t_start = t_start,
            t_stop = t_stop,
            discard = discard,
            sig_init = sig_init,
            sigma = sigma,
            V_prime = x -> -x.*(x.^2 .- 1),
            SM1 = true,
            Obs_noise = false)
        A[:,n] = autocov_con(signal[1,:],lag)
    end

    Autoc = mean(A,dims = 2)
    Autoc_var = var(A,dims = 2)
    vari ? [Autoc, Autoc_var] : Autoc
end

function Autocov(h_wf::Array{Float64,3}, Psi,
    steps::Int64;
    lag::UnitRange{Int64} = 0:999,
    bN = 100,
    t_start = t_start,
    t_stop = t_stop,
    discard = discard,
    sig_init = sig_init,
    sigma = sigma,
    vari = false,
    ALL = false)

    A = zeros(length(lag),bN)
    for n = 1:bN
        signal = redmodrun(h_wf, Psi,
            sigma_v = sigma_v,
            sig_init = sig_init,
            steps = steps,
            t_start = t_start,
            t_stop = t_stop,
            discard = discard
            )
        A[:,n] = autocov_con(signal[1,:],lag)
    end

    Autoc = mean(A,dims = 2)
    Autoc_var = var(A,dims = 2)
    out = vari ? [Autoc, Autoc_var] : Autoc
    ALL ? A : out
end

function redmodrun(
    h_wf, Psi;
    sigma_v = [1],
    sig_init = [1.5],
    steps = 10^5 + 1,
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    obs_noise::Array{Float64,2} = [1.0 0; 0 0]
    )

    d, nu, M_h = size(h_wf)
    d == size(sigma_v,1) || print("Dimension of signal and sigma
                                    do not agree")

    steps_tot = steps + discard
    h = (t_stop - t_start)/(steps)
    t = range(t_start,stop = t_stop,length = steps)

    # So the multiplication is defined when d=1
    d == 1 && (sigma_v = reshape(sigma_v,1,1))

    # Here we genereate the signal process.
    v = (obs_noise == [1.0 0; 0 0] ? randn(d,steps_tot) : obs_noise)
    signal_rm = complex(zeros(d,steps_tot))
    signal_rm[:,1] = sig_init
    PSI = complex(zeros(nu,steps_tot))
    for i = 1: steps_tot-1
        PSI[:,i] = Psi(signal_rm[:,i])
        signal_rm[:,i+1] = sum([real.(h_wf)[:,:,k+1]*PSI[:,i-k]
            for k = 0:min(i - 1,M_h - 1)]) + sqrt(h)*sigma_v*v[:,i+1]
    end

    sig_rm = signal_rm[:,discard+1:steps_tot]
end


function timeseries_plot(T,TS; title = "original plot")

    P1 = plot(T,TS,
        xlabel = "Time",
        ylabel = "Space",
        title = title * " time series" )

    P2 =  histogram(TS,
        xlabel = "Space",
        title = title * " histogram" )

    P = plot(P1,P2,
        layout = (1,2))
end

function No_of_trans(series,ϵ)
    S = (series .> ϵ) - (series .< -ϵ)
    well = (series[1] .> 0) - (series[1] .< 0)
    counter = 0
    for i = 2:length(S)
        if S[i] == -well
            counter += 1
            well   *= -1
        end
    end
    counter
end

mean_trans_time(series,Time;ϵ = .5) = No_of_trans(series,ϵ)/Time[end]



function my_hist(series, bin_num)
    series = reshape(series,length(series))
    sort!(series)
    bins = (series[end] - series[1])/bin_num.*(0:bin_num) .+ series[1]
    cu_hist = zeros(bin_num+1)
    for i = 1:bin_num
        cu_hist[i+1] = sum(series .< bins[i+1])
    end
    hist = zeros(bin_num)
    bin = zeros(bin_num)
    for i=1:bin_num
        hist[i] = cu_hist[i+1] - cu_hist[i]
        bin[i] = (bins[i] + bins[i+1])/2
    end
    plot(bin,hist,t=[:line],leg = :none)
end

function Model_filter(series; par = 55)
    L = par
    R_series = autocov_con(series,0:L)

    # Smoothing for z-spect-pred
    LL = Int(floor(L/2))
    lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
    lam2 = 2*(1 .- (LL+1:L)/L).^3
    lam = [lam1; lam2]

    R_series = R_series.*lam

    l = Scalar_CKMS_c(R_series);
end
