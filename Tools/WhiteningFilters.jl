module WhiteningFiltersScalar

using LinearAlgebra
using DSP
using FFTW
using ToeplitzMatrices: Toeplitz
using StatsBase: var

at = include("AnalysisToolbox.jl")
mr = include("WFMR.jl")

function get_whf(X::Vector{T};
    subrout = CKMS,
    nfft = 0,    
    par = 1500,             # Order of approximating Lauenat Poly
    M = par,               # number of output lags
    win = "Par",            # Type of smoother for autocov
    flags...) where T<: Number;

    steps = length(X)
    if subrout == CKMS
        nfft = nfft == 0 ? nextfastfft(steps) : nfft
    else
        nfft = max(nfft,par,2*M)
    end
    
    ## We need M ≤ par < steps
    M > steps - 1 && println("M is bigger than length of time series")
    M   = min(M, steps - 1)
    par = max(M,par)
    par = min(par,steps - 1)

    R_pred_smoothed = at.my_autocov(reshape(X,1,:); L = par, win)
 
    h_m, h_w = subrout(R_pred_smoothed; nfft, M, flags...)
end


function CKMS(R;          # fitst three are common to both subroutines
    nfft,                 # resolution of z-spectrum on unit cirlce
    M,                    # number of out put lags.
    tol_ckms = 0,     # the rest are subroutine specific parameters
    N_ckms = 10^4,
    verb = false)

    par = size(R,3) 
    
    # Compute coefficients of spectral factorization of z-spect-pred
    S⁻ = mr.spectfact_matrix_CKMS(R; ϵ = tol_ckms, N_ckms)

    Err  = S⁻[2] ###
    S⁻ = S⁻[1][:]                            # the model filter ###

    h_m = S⁻[1:M]

    S⁻ = nfft >= par ? [S⁻; zeros(eltype(S⁻), nfft - par)] : S⁻[1:nfft]

    fft!(S⁻)   # z-spectrum of model filter
    Sinv = (S⁻).^-1

    h_w = ifft(Sinv)[1:M];
    verb ? (h_m, h_w, Err) : (h_m, h_w)
end


function whf_cholesky(R;  # fitst three are common to both subroutine
    M = 100,                    # number of out put lags.
    nfft = 2M,            # The demesion of the covariance matrix we factor
    eps = 1e-4)
    
    println("You are using Cholesky")
    
    par = size(R,3)
    
    R = R[:]              # R will come in as 3D array but this is a scallar context so we flatten it
    
    R = nfft >= par ? [R; zeros(eltype(R), nfft - par)] : R[1:nfft]
                            
    σ = sqrt(R[1])*eps/sqrt(2)
    
    T = Array(Toeplitz(conj(R),R)) + σ^2*I
    T = conj(cholesky(T).L)
    
    h_m = reverse(T[end,end-M-1:end]) #T[nfft - M + 1:end,nfft - M + 1]
    
    T = inv(T);

    h_w = reverse(T[end,end-M-1:end]) #T[nfft - M + 1:end,nfft - M + 1]
    h_m, h_w
end

"""
    get_itr_whf is a function that computes whitening filters in an iteratatively

Examples:
~~~~~
x = pred[:]

h_m, h_w = whf.get_itr_whf(x; maxit = 3, par = 100, getter = whf.get_whf_B);

wx = filter_and_plot(h_w,x)
~~~~~
~~~~~
x = pred[:]

h_m, h_w = whf.get_itr_whf(x; maxit = 5, par = 1000, subrout = whf.whf_cholesky)

wx = filter_and_plot(h_w,x)
~~~~~
"""

function get_itr_whf(X::Vector{T};
    getter = get_whf,
    maxit = 5,
    par = 1500,
    flags_itr...) where T<: Number;
    
    h_w = [1]
    h_m = [1]
    println("Starting iterations")
    for i = 1 : maxit
        wx       = filt(h_w, X)
        Out      = getter(wx[maxit*par*(i-1)+1:end]; par, flags_itr...)[1:2];
        h_m      = conv(h_m, Out[1])
        h_w      = conv(h_w, Out[2])
    end
    println("Ending iterations")
    h_m, h_w
end


function Get_PREDhash(X::Vector{T};M_out = 30) where T<: Number;
    steps = length(X)
    
    Xh = zeros(ComplexF64,steps-M_out+1,M_out)
    for m = 1:M_out
        Xh[:,m] = @view X[m:(steps + m - M_out)]
    end
    
    Xh
end

function get_whf_B(X::Vector{<:Number}; par = 50, eps = 1e-6)
    # 1. Form X-hash
    M_out = par
    
    Xh = Get_PREDhash(X; M_out)

    # 2. Factor
    # a small amount of regularization is needed so that the matrix we factor is PD.
    σ = eps^2*var(X)/2
    W = Xh'Xh/size(Xh,1) + σ*I
    println(cond(W))

    W = conj(cholesky(W).L)
    println(cond(W))
    
    h_m = reverse(W[end,:])
    
    W = inv(W)
    
    h_w = reverse(W[end,:])
    
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

# function get_filters(x;par = 5000, win = "Par")
#     X = size(x,2) > 1 ? x : X = reshape(x,1,:)
#     tol_ckms = 0; N_ckms = 10000; # par = 
#     nfft = nextfastfft(size(X,2))
    

#     d, steps = size(X)

#     nfft = nfft == 0 ? nextfastfft(steps) : nfft
#     nffth = nfft ÷ 2
#     L = min(par,steps-1)

#     R_pred_smoothed = mr.matrix_autocov_seq(X; L, win);

#     # Compute coefficients of spectral factorization of z-spect-pred
#     S⁻ = mr.spectfact_matrix_CKMS(R_pred_smoothed; ϵ = tol_ckms, N_ckms)

#     Err  = S⁻[2] ###
#     S⁻ = S⁻[1]                            # the model filter ###
    
#     h_m = S⁻[1:par]
    
#     S⁻ = nfft >= L+1 ? cat(dims = 3,S⁻,zeros(d,d,nfft - L - 1)) :
#                             (@view S⁻[:,:,1:nfft])

#     fft!(S⁻, 3)   # z-spectrum of model filter
#     Sinv = (S⁻[:]).^-1

#     h_w = ifft(Sinv)[1:par];
#     h_m, h_w
# end

