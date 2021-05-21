module WhiteningFiltersScalar

using LinearAlgebra
using DSP: nextfastfft
using FFTW

at = include("AnalysisToolbox.jl")
mr = include("WFMR.jl")

function get_whf(X::Vector{T};
    subrout = CKMS,
    nfft = 0,
    M = 1000,               # number of output lags
    par = 1500,             # Order of approximating Lauenat Poly
    win = "Par",            # Type of smoother for autocov
    flags...) where T<: Number;

    steps = length(X)
    if subrout == CKMS
        nfft = nfft == 0 ? nextfastfft(steps) : nfft
    else
        nfft = max(nfft,par,2*M)
    end
    
    M > steps - 1 && println("M is bigger than length of time series")
    M   = min(M, steps - 1)
    par = max(M,par)
    par = min(par,steps - 1)

    R_pred_smoothed = mr.matrix_autocov_seq(reshape(X,1,:); L = par, win)

    h_m, h_w = subrout(R_pred_smoothed; nfft, M, flags...)
end


function CKMS(R;          # fitst three are common to both subroutines
    nfft,                 # resolution of z-spectrum on unit cirlce
    M,                    # number of out put lags.
    tol_ckms = 1e-10,     # the rest are subroutine specific parameters
    N_ckms = 10^4)

    par = size(R,3) 
    
    # Compute coefficients of spectral factorization of z-spect-pred
    S⁻ = mr.spectfact_matrix_CKMS(R; ϵ = tol_ckms, N_ckms)

    Err  = S⁻[2] ###
    S⁻ = S⁻[1][:]                            # the model filter ###

    h_m = S⁻[1:M]

    S⁻ = nfft >= par ? [S⁻; zeros(eltype(S⁻), nfft - par)] : S⁻

    fft!(S⁻)   # z-spectrum of model filter
    Sinv = (S⁻).^-1

    h_w = ifft(Sinv)[1:M];
    h_m, h_w
end


function whf_cholesky(R;  # fitst three are common to both subroutines
    nfft,                 # resolution of z-spectrum on unit cirlce
    M,                    # number of out put lags.
    eps = 1e-4)
    
    par = size(R,3)
    
    R = R[:]              # R will come in as 3D array but this is a scallar context so we flatten it
    
    R = nfft >= par ? [R; zeros(eltype(R), nfft - par)] : R
                            
    σ = sqrt(R[1])*eps/sqrt(2)
    
    T = Array(Toeplitz(conj(R),R)) + σ^2*I
    T = conj(cholesky(T).L)
    
    h_m = T[pad - M + 1:end,pad - M + 1]
    
    T = inv(T);

    h_w = T[pad - M + 1:end,pad - M + 1]
    h_m, h_w
end

end # module







# ## Old stuff just in case

# ### Original



# function get_whf(X::Array{T,2};
#     par = 1500,
#     nfft = nextfastfft(size(X,2)),
#     win = "Par",
#     tol_ckms = 1e-10,
#     N_ckms = 10^4,
#     verb = false,
#     model = false) where T <: Number

#     d, steps = size(X)

#     nfft = nfft == 0 ? nextfastfft(steps) : nfft
#     nffth = nfft ÷ 2
#     L = min(par,steps-1)

#     R_pred_smoothed = @timed mr.matrix_autocov_seq(X; L, win)
#     if verb
#         println("Time taken for autocov: ", R_pred_smoothed.time)
#         println("Bytes Allocated: ", R_pred_smoothed.bytes)
#     end

#     # Compute coefficients of spectral factorization of z-spect-pred
#     S⁻ = @timed mr.spectfact_matrix_CKMS(R_pred_smoothed.value; ϵ = tol_ckms, N_ckms, verb)
#     if verb
#         println("Time taken for spectfact: ",S⁻.time)
#         println("Bytes Allocated: ",S⁻.bytes)
#     end

#     Err  = S⁻.value[2] ###
#     S⁻ = S⁻.value[1]                            # the model filter ###

#     S⁻ = nfft >= L+1 ? cat(dims = 3,S⁻,zeros(d,d,nfft - L - 1)) :
#                             (@view S⁻[:,:,1:nfft])

#     fft!(S⁻, 3)                                 # z-spectrum of model filter

#     S⁻inv = complex(zeros(d,d,nfft))
#     for i = 1 : nfft
#         S⁻inv[:,:,i] = inv(@view S⁻[:,:,i])
#     end                                             # the final S_pred⁺

#     h_whf = ifft(S⁻inv, 3)
# #     model ? [h_whf, h_mf] : h_whf
#     verb ? [h_whf, Err] : h_whf
# end

