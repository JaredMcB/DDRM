
using StatsPlots
using DataFrames
using StatsBase
using Statistics
using DSP
using FFTW

include("DWOL_eqidist_sampler.jl")

function DataGen_Langevin_FE(
    steps = 10^5 + 1;
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [1],
    d = 1,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false
    )

    d == size(sigma,1) || print("Dimension of signal and sigma do not agree")

    steps_tot = steps + discard
    h = (t_stop - t_start)/(steps)
    t = range(t_start,stop = t_stop,length = steps)

    # So the multiplication is defined when d=1
    d == 1 && (sigma = reshape(sigma,1,1))

    # Here we genereate the signal process.
    e = randn(d,steps_tot)
    signal = zeros(d,steps_tot)
    signal[:,1] = sig_init
    for n = 2 : steps_tot
        signal[:,n] = signal[:,n-1] + h*V_prime(signal[:,n-1]) + sqrt(h)*sigma*e[:,n]
    end

    # Here we discard the transient part
    sig = signal[:,discard+1:steps_tot]

    # sig_minus_1 is returned as it is needed in produceing the predictor series
    # which is offset by one.
    result = SM1 ? [signal[:,discard] sig] : sig
    result = Obs_noise ? [result, e] : result
end




function DataGen_DWOL(
    steps = 10^5 + 1;
    scheme = "FE",
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [1],
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false,
    d = 1,
    e = randn(d,steps + discard)
    )

    if discard == 0
        sig_init = [DWOL_dist_samp(1,sigma[1])]
    end

    size(sig_init,1) == size(sigma,1) || print("Dimension of initial value "*
                                               "and sigma do not agree")

    steps_tot = steps + discard
    Δt = (t_stop - t_start)/(steps-1)

    # So the multiplication is defined when d=1
    d == 1 && (sigma = reshape(sigma,1,1))

    # Here we genereate the signal process.

    signal = zeros(d,steps_tot)
    signal[:,1] = sig_init
    if scheme == "FE"
        for n = 1 : steps_tot-1
            signal[:,n+1] = signal[:,n] + Δt*V_prime(signal[:,n]) + sqrt(Δt)*sigma*e[:,n+1]
        end
    elseif scheme == "T2"
        for n = 1 : steps_tot-1
            signal[:,n+1] = signal[:,n] .+ Δt*(-signal[:,n].^3 .+ signal[:,n]) .+
                            Δt^2/2*( signal[:,n].*(signal[:,n].^2 .- 1).*(signal[:,n].^2 .- 1) ) .+
                            sqrt(Δt)*sigma*e[:,n+1]
        end
    elseif scheme == "EM"
        for n = 1 : steps_tot-1
            k1 = signal[:,n] .+ Δt/2*V_prime(signal[:,n])
            signal[:,n+1] = signal[:,n] .+ Δt*V_prime(k1) .+
                            sqrt(Δt)*sigma*e[:,n+1]
        end
    elseif scheme == "ET"
        for n = 1 : steps_tot-1
            k0 = V_prime(signal[:,n])
            k1 = signal[:,n] .+ Δt*k0
            signal[:,n+1] = signal[:,n] .+ Δt/2*V_prime(k1) .+ Δt/2*K0 .+
                            sqrt(Δt)*sigma*e[:,n+1]
        end
    elseif scheme == "RK4"
        for n = 1 : steps_tot-1
            k1 = signal[:,n]
            k2 = signal[:,n] .+ Δt/2*V_prime(k1)
            k3 = signal[:,n] .+ Δt/2*V_prime(k2)
            k4 = signal[:,n] .+ Δt*V_prime(k3)
            signal[:,n+1] = signal[:,n] .+ Δt/6*( V_prime(k1) .+ 2*V_prime(k1) .+
                            2*V_prime(k1) .+ V_prime(k1) ) .+
                            sqrt(Δt)*sigma*e[:,n+1]
        end
    end

    # Here we discard the transient part
    sig = signal[:,discard+1:steps_tot]

    # sig_minus_1 is returned as it is needed in produceing the predictor series
    # which is offset by one.
    result = SM1 ? [signal[:,discard] sig] : sig
    result = Obs_noise ? [result, e] : result
end

function DWOL_dist_samp(N::Int64 = 1, σ = [.35])
    # sampling

    μ = σ[1]^2/2
    p(z) =  exp(-(x^2-1)^2/μ)
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

function DataGen_Langevin_RK4(
    steps = 10^5 + 1;
    t_start = 0,
    t_stop = 10^3,
    discard = 10^4,
    sig_init = [1.5],
    sigma = [1],
    d = 1,
    V_prime = x -> -x.*(x.^2 .- 1),
    SM1 = false,
    Obs_noise = false
    )

    d == size(sigma,1) || print("Dimension of signal and sigma do not agree")

    steps_tot = steps + discard
    h = (t_stop - t_start)/(steps)
    t = range(t_start,stop = t_stop,length = steps)

    # So the multiplication is defined when d=1
    d == 1 && (sigma = reshape(sigma,1,1))

    # Here we genereate the signal process.
    e = randn(d,steps_tot)
    signal = zeros(d,steps_tot)
    signal[:,1] = sig_init
    for n = 2 : steps_tot

        k1 =  V_prime(signal[:,n])
        k2 =  V_prime(signal[:,n] + h*k1/2)
        k3 =  V_prime(signal[:,n] + h*k2/2)
        k4 =  V_prime(signal[:,n] + h*k3)

        signal[:,n] = signal[:,n-1] +
                    h/6*(k1 + 2*k2 + 2*k3 + k4) +
                    sqrt(h)*sigma*e[:,n]
    end

    # Here we discard the transient part
    sig = signal[:,discard+1:steps_tot]

    # sig_minus_1 is returned as it is needed in produceing the predictor series
    # which is offset by one.
    result = SM1 ? [signal[:,discard] sig] : sig
    result = Obs_noise ? [result, e] : result
end



function get_wf(signal,Psi;
    M_out = 20,
    rl = true)
    # We would like a presample since we want the
    # times series to be offset by one.

    sig = signal[:,2:end]
    d, steps = size(sig)
    nu = size(Psi(zeros(d,1)),1)

    pred = complex(zeros(nu, steps))
    for n = 1:steps
        pred[:,n] = Psi(signal[:,n])
    end

    h_wf = rl ? real(vector_wiener_filter_fft(pred, sig, M_out)
                    ) : vector_wiener_filter_fft(pred, sig, M_out)
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

function analyse_h_ens(h_wf_ens; plt = true)

    h_wf_mean = mean(h_wf_ens,dims = 4)[:,:,:,1]
    h_wf_var = var(h_wf_ens,dims = 4)[:,:,:,1]

    if plt
        T = h_wf_ens[1,1,:,:]'
        T[:,1] .-= 1

        S = h_wf_ens[1,2,:,:]'

        dfT, dfS = DataFrame(T), DataFrame(S)
        p1 = @df dfT boxplot(T,
            marker=(0.3,:orange,stroke(.5)),
            alpha=0.75,
            leg = :none,
            title = "First coefficients")
        p2 = @df dfS boxplot(S,
            marker=(0.3,:orange,stroke(.5)),
            alpha=0.75,
            leg = :none,
            title = "Second coefficients")

        P = plot(p1,p2, layout = (2,1))
    end

    plt == true ?  [h_wf_mean, h_wf_var, P] :  [h_wf_mean, h_wf_var]
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

function autocov_con(x::AbstractVector{<:Number},lags::UnitRange{Int})
    lx = size(x,1)
    if maximum(lags) > lx
        print("lag cannot be greater than lenght of series")
        lags = filter(x -> abs(x) < lx, lags)
    end

    x .-= mean(x)
    A = conv(x,conj(reverse(x)))/lx
    A = [A[k + lx] for k in lags]
end

function auto_times(x::AbstractVector{<:Real};plt = false)
    lx = size(x,1)
    L = minimum([lx, 10^6])

    lags = 0:L
    A = autocov_con(x,lags)

    end_int = try
                findall(A.<0)[1] - 1
            catch e
                if isa(e, BoundsError)
                    L
                end
            end
    end_exp = Int64(round(end_int/3))

    A_mat = [ones(end_exp,1) reshape(1:end_exp,end_exp,1)]
    b = inv(A_mat'*A_mat)*A_mat'*log.(A[1:end_exp])
    τ_exp = -1/b[2]

    A ./= A[1]
    τ_int = .5 + sum(A[2:end_int])
    if plt
        P = plot(0:(end_int-1),log.(A[1:end_int]),
            ylabel = "Log of Autocov",
            xlabel = "Lags",
            label = "log(A)")
        P = plot!(P,0:end_exp-1,A_mat*b,
            label = "linear approx of log(A)")
        end
    plt ? [τ_exp, τ_int, P] : [τ_exp, τ_int]
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

function emp_cdf(series;
    plt_only = true)
    l = length(series)
    series = reshape(series,l)
    sort!(series)

    bw = 2*iqr(series)/l^(1/3)
    bn = Int64(ceil((series[end] - series[1])/bw))
    b_pts = bw*(0:bn) .+ series[1]

    cdf = zeros(bn+1)
    for i = 1:bn
        cdf[i+1] = sum(series .< b_pts[i+1])
    end
    cdf /= l

    P = plot(b_pts,cdf)
    plt_only ? P : [cdf,b_pts,P]
end

function emp_pdf(series;
    plt_only = true)

    P = emp_cdf(series,plt_only = false)
    cdf = P[1]
    b_pts = P[2]
    bn = length(cdf) - 1
    bw = b_pts[2] - b_pts[1]

    pdf = zeros(bn)
    b_midpts = zeros(bn)
    for i=1:bn
        pdf[i] = cdf[i+1] - cdf[i]
        b_midpts[i] = (b_pts[i] + b_pts[i+1])/2 # Midpoint as bin location
    end
    pdf /= bw

    P = plot(b_midpts,pdf)
    plt_only ? P : [pdf,b_midpts,P]
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

function Scalar_CKMS_f(s_ana; N = 100)
    # N is the number of point used in the fourier transorm it is the number of
    # terms of the approximating lauren polynomial.

    S = [s_ana(exp(im*2*pi*j/N)) for j = 0 : N-1];

    S_fft = fft(S)/N;

    Ne = N - Int64(floor(N/2)); # This is to get only the casual coefficients

    S_fft_1 = S_fft[1:Ne]; # These are only the casual coefficeints

    m = Ne - 1
    F = [[zeros(1,m-1); I] zeros(m)]
    h = [zeros(1,m-1) 1]
    NN = reverse(S_fft_1[2:end]);

    K = reshape(NN,m,1)
    L = reshape(NN,m,1)
    Re = reshape([S_fft_1[1]],1,1)
    Rr = reshape([S_fft_1[1]],1,1)
    for i = 1:100
        hL = h*L
        FL = F*L

        K_new = K - FL/Rr*hL'
        L_new = FL - K/Re*hL
        Re_new = Re - hL/Rr*hL'
        Rr_new = Rr - hL'/Re*hL

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end
    k = K/Re
    re = Re[1]

    l = sqrt(re)*[k[n] for n = m:-1:1]

    s_plus_num(z) = sqrt(re) + sum([l[n]*z^(-n) for n = 1:m]);

    return s_plus_num, l
end

function Scalar_CKMS_c(R)

    m = length(R) - 1
    F = [[zeros(1,m-1); I] zeros(m)]
    h = [zeros(1,m-1) 1]
    NN = reverse(R[2:end]);

    K = reshape(NN,m,1)
    L = reshape(NN,m,1)
    Re = reshape([R[1]],1,1)
    Rr = reshape([R[1]],1,1)
    for i = 1:100
        hL = h*L
        FL = F*L

        K_new = K - FL/Rr*hL'
        L_new = FL - K/Re*hL
        Re_new = Re - hL/Rr*hL'
        Rr_new = Rr - hL'/Re*hL

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end
    k = K/Re
    re = Re[1]

    l = [sqrt(abs(re))]
    [l; sqrt(abs(re))*[k[n] for n = m:-1:1]]
end

function crosscov_con(x,y,lags)
    lx = size(x,1); ly = size(y,1)
    lx == ly || throw(DimensionMismatch("lengths of inputs must match"))

    if maximum(lags) >= lx
        print("absolute lag must be less than lenght of series")
        lags = filter(x -> abs(x) < lx, lags)
    end

    x = x .- mean(x)
    y = y .- mean(y)

    A = conv(x,conj(reverse(y)))/lx
    A = [A[k + lx]  for k in lags]
end
