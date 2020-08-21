using StatsBase

include("SFbyCKMS.jl")

function wiener_filter_NT(pred,sig,M_in,N,M_out)
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

    #   Now we us FFT to complete the calculations
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

    M = M_out
    h = H[1:M]
end
