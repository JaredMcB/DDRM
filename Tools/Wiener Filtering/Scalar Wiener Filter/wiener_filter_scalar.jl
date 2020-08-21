################################################
#
# Title: wiener_filter_scalar.jl
# Author: Jared McBride (11-15-2019)
#
# Here we find the Wiener filter for two scalar
# processes, the signal and the predictors (or
# observations).
# The inputs are
#   sig : the timeseries we are approximating
#   pred : the timeseries with which we will
#       be approximating the signal.
#   M : the degree of the approximating Laurent
#       polynomials used in construction of the
#       spectra.
# The output is
#   h : the sequence of length 2M+1.
#
##############################################
# Spectral Method
using StatsBase

# Time Domain
using JuMP
using OSQP

include("SFbyCKMS.jl")

function z_crossspect(sig,pred,L; win = "Bart")
    lags = -L:L;
    C = crosscov(sig,pred,lags)
    if win == "Bar"
        lam = 1 .- (0:L)/L
    elseif win == "Tuk"
        lam = .5*(1 .+ cos.(pi/L*(0:L)))
    elseif win == "Par"
        LL = Int(floor(L/2))
        lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
        lam2 = 2*(1 .- (LL+1:L)/L).^3
        lam = [lam1; lam2]
    else
        lam = ones(L+1)
    end

    z_crossspect_num(z) = sum([lam[abs(i) + 1]*C[L+1+i]*z^(-i) for i = -L : L])
end

function wiener_filter_Nu(pred, sig, M_out = 513; par::Int64 = 55, Nex::Int64 = 2^10)
    Theta = 2pi*(0:Nex-1)/(Nex-1);
    Z = map(th -> exp(im*th),Theta);

    M = par
    R_pred = autocov(pred,0:M)

    #     Smoothing
    LL = Int(floor(M/2))
    lam1 = 1 .- 6*((0:LL)/M).^2 .+ 6*((0:LL)/M).^3
    lam2 = 2*(1 .- (LL+1:M)/M).^3
    lam = [lam1; lam2]

    R_pred = R_pred.*lam
    l = Scalar_CKMS_c(R_pred)

    z_crossspect_sigpred_num = z_crossspect(sig, pred, par, win = "Par")
    z_spect_pred_minus_num(z) = sum([l[i+1]z^(-i) for i = 0:M])
    z_spect_pred_plus_num(z) = z_spect_pred_minus_num(z)'

    S_sigpred_invS_plus_fft_num = fft(z_crossspect_sigpred_num.(Z)./z_spect_pred_plus_num.(Z))/Nex

    Nexh = Int(floor(Nex/2))
    S_sigpred_invS_plus_fft_plus_num = [S_sigpred_invS_plus_fft_num[1]; zeros(Nexh - 1); S_sigpred_invS_plus_fft_num[Nexh + 1:end]];

    S_sigpred_invS_plus_plus_num = ifft(S_sigpred_invS_plus_fft_plus_num)*Nex

    H_num = S_sigpred_invS_plus_plus_num./z_spect_pred_minus_num.(Z);
    h_num_raw = fft(H_num)/Nex
    h_num = [h_num_raw[1]; reverse(h_num_raw)[1:Nexh]]
    h_num[1:M_out]
end


function wiener_filter_TD(pred,sig,M)
    # pred are the predictors
    # sig is the sequence to be approximated
    # M is the degree of the approximating Laurent polynomials
    N = length(sig)

    # Approximate autocorr and crosscor (Laurent degree M)
    function Loss(h)
        S = 0
        m = length(h)
        for n = 1:N

            nvm = minimum([n m])
            con = dot(reverse(h[1:nvm]),pred[(n - nvm +1) : n])
            S = S + (sig[n] - con)^2/N
        end
        S
    end

    model = Model(with_optimizer(OSQP.Optimizer))

    @variable(model, h[1:M])

    @objective(model, Min, Loss(h))

    optimize!(model)

    h_TD = value.(h)
end

function wiener_filter(pred,sig,M_in,N,M_out)
    # pred are the predictors
    # sig is the sequence to be approximated
    # M is the degree of the approximating Laurent polynomials

    # Approximate autocorr and crosscor (Laurent degree M)
    M = M_in
    R_pred    = autocov(pred,0:M)
    R_predsig = crosscov(pred,sig, -M:M)

    #     Smoothing
    LL = Int(floor(M/2))
    lam1 = 1 .- 6*((0:LL)/M).^2 .+ 6*((0:LL)/M).^3
    lam2 = 2*(1 .- (LL+1:M)/M).^3
    lam = [lam1; lam2]

    R_pred = R_pred.*lam
    R_predsig = R_predsig.*[reverse(lam[2:M+1]); lam]

    l = Scalar_CKMS_c(R_pred);

    z_spect_minus_num(z) = sum([l[i+1]z^(-i) for i = 0:M])
    z_spect_plus_num(z) = z_spect_minus_num(z)'

    z_spect_sigpred(z) = sum([R_predsig[i]*z^(M+1 - i) for i = 1: 2*M + 1])

    N
    Nh = Int(floor(N/2))
    Z_N = exp.(2*pi*im/N*(0:N-1))
    Theta_N = 2*pi/N*(0:N-1)
    S_plus_invS_sigpred_fft = fft(z_spect_sigpred.(Z_N)./
                                  z_spect_plus_num.(Z_N))/N
    S_plus_invS_sigpred_coef = S_plus_invS_sigpred_fft[1:Nh];

    S_plus_invS_sigpred_coef_N = [S_plus_invS_sigpred_coef; zeros(N-Nh)]
    z_spect_invS_sigpred_plus = ifft(S_plus_invS_sigpred_coef_N)*N

    H =  fft(z_spect_invS_sigpred_plus./
             z_spect_minus_num.(Z_N))/N

    M = minimum([M_out Nh])
    h = H[1:M]
end
