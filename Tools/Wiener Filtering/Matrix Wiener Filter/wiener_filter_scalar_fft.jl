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

include("SFbyCKMS.jl")

function z_crossspect_fft(sig,pred,L,Nex; win = "Par")
    Nexh = Int(floor(Nex/2))
    lags = -L:L;
    c = size(sig) == size(pred) ? crosscov(sig,pred,lags) :
        error("sig and pred must be same size")

    # Smoothing window
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
    Lam = [lam[L+1:-1:2]; lam]

    C_smoothed = Lam .* c;
    C_padded = [zeros(Nex - Nexh - L); C_smoothed; zeros(Nexh - L - 1)];
    C = fftshift(C_padded);

    z_crossspect_num_fft = fft(C)
end

function wiener_filter_fft(pred, sig, M_out = 100; par::Int64 = 55, Nex::Int64 = 2^10)

    Nexh = Int(floor(Nex/2))

    L = par
    R_pred = autocov(pred,0:L)

    # Smoothing for z-spect-pred
    LL = Int(floor(L/2))
    lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
    lam2 = 2*(1 .- (LL+1:L)/L).^3
    lam = [lam1; lam2]

    R_pred = R_pred.*lam

    # Compute coefficients of spectral factorization of z-spect-pred
    l = Scalar_CKMS_c(R_pred);
    l_pad = [l; zeros(Nex - (L+1))]
    z_spect_pred_minus_num_fft = fft(l_pad)
    z_spect_pred_plus_num_fft = conj.(z_spect_pred_minus_num_fft)

    # Compute z-cross-spectrum of sigpred
    z_crossspect_sigpred_num_fft = z_crossspect_fft(sig, pred, par, Nex, win = "Par");

    # This computes the impule response (coefficeints of z) for S_{yx}/S_x^+
    S_sigpred_overS_plus_fft_num_fft = ifft(z_crossspect_sigpred_num_fft./z_spect_pred_plus_num_fft)

    # Extracts causal part coefficinets of S_{yx}/S_x^+, {S_{yx}/S_x^+}_+
    S_sigpred_overS_plus_fft_plus_num_fft = [S_sigpred_overS_plus_fft_num_fft[1: Nexh]; zeros(Nex - Nexh)]

    # Computes causal part of S_{yx}/S_x^+, {S_{yx}/S_x^+}_+
    S_sigpred_overS_plus_plus_num_fft = fft(S_sigpred_overS_plus_fft_plus_num_fft);

    # Obtain transfer function H by dividing {S_{yx}/S_x^+}_+ by S_x^-
    H_num_fft = S_sigpred_overS_plus_plus_num_fft ./ z_spect_pred_minus_num_fft

    # Extrct tranferfunction coeffifcients (impulse responce of Weiner filter)
    h_num_raw_fft = ifft(H_num_fft)

    # Truncate
    h_num_fft = h_num_raw_fft[1:M_out]
end
