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
    signal_rm = zeros(d,steps_tot)
    signal_rm[:,1] = sig_init
    PSI = zeros(nu,steps_tot)
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
