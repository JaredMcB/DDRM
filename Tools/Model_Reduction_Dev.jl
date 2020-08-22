using FFTW
using LinearAlgebra
using DSP: conv, nextfastfft
using Polynomials
using StatsBase
using SparseArrays



"""
    `get_wf` provides the causal wiener filter that best approximates `signal`
in terms of `Psi(signal[:,i-1])`. it has built in the one step prediction. First,
it generates the predictor series. and then employs the `vector_wiener_filter_fft`
function to obtain the Wiener filter. Inputs are the signal, a 2d array, d-by-N
where d is the dimension of the state space an N is the number of point. There
are two key word arguments, `M_out` specifies the number of coefficents desired
from the Weiner filter. And `rl` which if true just take the real part of the
Wiener filter as the output.
"""

function get_wf(signal, Psi;
    M_out = 20,
    n = 3, p = 1500, par = 1500,
    rl = true,
    PI = false,
    rtol = 1e-6)
    # We would like a presample since we want the
    # times series to be offset by one.

    sig = signal[:,2:end]
    d, steps = size(sig)
    nu = size(Psi(zeros(d,1)),1)

    pred = complex(zeros(nu, steps))
    for n = 1:steps
        pred[:,n] = Psi(signal[:,n])
    end

    h_wf = vector_wiener_filter_fft(pred, signal[:,1:end-1], M_out,
            n = n, p = p, par = par, PI = PI, rtol = rtol)

    h_wf = rl ? real(h_wf) : h_wf
end



"""
    spectfact_matrix_CKMS take in m+1 coefficents of a (d x d)-matrix-values Laurent
polynomial of the form
        M(z) = P[:,:,1] + sum_{i=1}^{m} P[:,:,i+1]*z^{-i} + P[:,:,i+1]'z^i
and produce the first m+1 coefficients of the (dxd)-matrix valued polynomial
(in z^{-1})
        L(z) = l[:,:,1] + sum_{i=1}^{m} l[:,:,i+1]*z^{-i}
which satisfies
        M(z) = L(z)*L(z^{-1}}')'
    For this function the input is P and the output is l.
"""

function spectfact_matrix_CKMS(P; N_ckms = 200)
    d = size(P)[1];
    m = size(P)[3] - 1

    NN = reverse(P[:,:,2:end],dims = 3)
    Re = Rr = p0 = P[:,:,1]

    F = sparse([[zeros(d,d*(m-1)); I] zeros(d*m,d)])
    h = sparse([zeros(d,d*(m-1)) I])

    K = complex(zeros(d*m,d))
    for i = 0 : m-1
        K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
    end
    L = K

    # spectfactLog = zeros(4,N_ckms)

    for i = 1:N_ckms
        hL = h*L; FL = F*L

        K_new = K - FL/Rr*hL'
        L_new = FL - K/Re*hL
        Re_new = Re - hL/Rr*hL'
        Rr_new = Rr - hL'/Re*hL

        # spectfactLog[:,i] = [cond(Rr),
        #                      cond(Re),
        #                      norm(K - K_new),
        #                      norm(L - L_new)]

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end

    k = K/Re
    re = Re

    sqrt_re = sqrt(re)

    l = complex(zeros(d,d,m+1))
    l[:,:,1] = sqrt_re;
    for i = m-1:-1:0
        l[:,:,m-i+1] = k[d*i + 1: d*(i+1),:]*sqrt_re
    end

    # save("Data\\CKMS_dat.jld",
    #     "spectfactLog",
    #     spectfactLog)

    l
end

function spectfact_matrix_CKMS_pinv(P; N_ckms = 200, rtol = 1e-6)
    d = size(P)[1];
    m = size(P)[3] - 1

    NN = reverse(P[:,:,2:end],dims = 3)
    Re = Rr = p0 = P[:,:,1]

    F = sparse([[zeros(d,d*(m-1)); I] zeros(d*m,d)])
    h = sparse([zeros(d,d*(m-1)) I])

    K = complex(zeros(d*m,d))
    for i = 0 : m-1
        K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
    end
    L = K

    for i = 1:N_ckms
        hL = h*L; FL = F*L

        K_new = K - FL*pinv(Rr,rtol = rtol)*hL'
        L_new = FL - K*pinv(Re,rtol = rtol)*hL
        Re_new = Re - hL*pinv(Rr,rtol = rtol)*hL'
        Rr_new = Rr - hL'*pinv(Re,rtol = rtol)*hL

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end

    k = K*pinv(Re,rtol = rtol)
    re = Re

    sqrt_re = sqrt(re)

    l = complex(zeros(d,d,m+1))
    l[:,:,1] = sqrt_re;
    for i = m-1:-1:0
        l[:,:,m-i+1] = k[d*i + 1: d*(i+1),:]*sqrt_re
    end
    l
end


function _smoother(n=4,p=5; ty = "bin")
    if ty == "bin"
        μ = Polynomial(ones(p+1))
        μ_sq = μ^(2n)/(p+1)^(2n)
        μ_c = coeffs(μ_sq)
    elseif ty == "ave"
        μ = Polynomial(ones(2p+1))
        μ_sq = μ^(n)/(2p+1)^(n)
        μ_c = coeffs(μ_sq)
    else
        μ_c = ones(2n*p+1)
    end
    μ_c
end

function _window(L; win = "Par",two_sided = true)
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
    two_sided ? [lam[L+1:-1:2]; lam] : lam
end

"""
    my_crosscov
I don't remember why I wrote this or if it has any advantage over some builtin
function.
"""
function _crosscov_con(x::AbstractVector{<:Number},
                      y::AbstractVector{<:Number},
                      lags)
    lx = size(x,1)
    ly = size(y,1)
    lx == ly || throw(DimensionMismatch("series must be same length"))

    if maximum(lags) > lx
        println("lag cannot be greater than lenght of series")
        lags = filter(x -> abs(x) < lx, lags)
    end

    x .-= mean(x)
    y .-= mean(y)
    C = conv(x,conj(reverse(y)))/lx
    C = [C[k + lx] for k in lags]
end

function _crosscov_dot(x::AbstractVector{<:Number},
                      y::AbstractVector{<:Number},
                      lags)
    L = min(length(x),length(y))
    m = length(lags)

    zx = x .- mean(x)
    zy = y .- mean(y)

    C = zeros(Complex,m)
    for k = 1:m
        l = lags[k]
        C[k] = ( l >= 0 ? dot(zx[1+l : L],zy[1 : L-l]) : dot(zx[1 : L+l],zy[1-l : L]))/L
    end
    C
end

function my_crosscov(x::AbstractVector{<:Number},
                     y::AbstractVector{<:Number},
                     lags)
    length(lags) > 1000 ? _crosscov_con(x,y, lags) : _crosscov_dot(x,y, lags)
end

function my_crosscor(x::AbstractVector{<:Number},
                     y::AbstractVector{<:Number},
                     lags)
    my_crosscov(x,y,lags)/my_crosscov(x,y,0:0)[1]
end


"""
    z_crossspect_fft(sig::Array{Complex{Float64},2},pred::Array{Complex{Float64},2};
n = 3, p = 2500, win = "Par")

The output has length nfft = nextfastfft(steps)
"""
function z_crossspect_fft(
    sig,
    pred::Array{T,2} where T <: Number;
    nfft = 0,
    n = 3,
    p = 2500,
    win = "Par")

    ## sig = d x steps, pred = nu x steps
    d, stepsx = size(sig)
    nu, stepsy = size(pred)

    stepsx == stepsy || print("sig and pred are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])
    nfft = nfft == 0 ? nextfastfft(steps) : nfft
    steps == nfft || println("adjusted no. of steps from $steps to $nfft")
    steps = nfft

    z_spect_mat = zeros(Complex, d, nu, nfft)
    for i = 1 : d
        for j = 1 : nu
            z_spect_mat[i,j,:] = z_crsspect_scalar(sig[i,:],pred[j,:],
                                                  nfft = nfft, n = n, p = p)
        end
    end
    z_spect_mat
end

function z_spect_scalar(sig; n = 3, p=100, ty = "ave")
    μ = _smoother(n,p,ty = "ave")

    siz = length(sig)
    nfft = nextfastfft(siz)
    sig_pad = [sig; zeros(nfft - siz)]

    peri = abs.(fft(sig_pad)).^2/nfft
    peri_pad = [peri[end - p*n + 1 : end]; peri; peri[1:p*n]]
    z_spect_smoothed = conv(μ,peri_pad)[2n*p:2n*p+nfft-1]
end

"""
z_crsspect_scalar has output of size nfft
"""
function z_crsspect_scalar(sig,pred; nfft = 0, n = 3, p=100, ty = "ave")
    μ = _smoother(n,p,ty = "ave")

    l_sig = length(sig)
    l_pred = length(pred)
    l_sig == l_pred || println("sizes must be the same, taking min and truncating")
    l = min(l_sig,l_pred)

    nfft = nfft == 0 ? nfft = nextfastfft(l) : nfft

    # nfft == l || println("adjusted size from $l to $nfft")
    sig_pad = l < nfft ? [sig[1:l]; zeros(nfft - l)] : sig[1:nfft]
    pred_pad = l < nfft ? [pred[1:l]; zeros(nfft - l)] : pred[1:nfft]

    fftsig = fft(sig_pad)
    fftpred = conj(fft(pred_pad))

    peri = fftsig .* fftpred / nfft
    peri_pad = [peri[end - p*n + 1 : end]; peri; peri[1:p*n]]
    z_crsspect_smoothed = conv(μ,peri_pad)[2n*p:2n*p+nfft-1]
end



"""
    vector_wiener_filter_fft

"""
function vector_wiener_filter_fft(
    sig,
    pred::Array{T,2} where T <: Number,
    M_out = 20;
    par::Int64 = 1500,
    win = "Par",
    n = 3,
    p = 1500,
    PI = true,
    rtol = 1e-6
    )

    d, stepsy = size(sig)
    nu, stepsx = size(pred)

    stepsx == stepsy || print("X and Y are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])
    nfft = nextfastfft(steps)
    nffth = Int(floor(nfft/2))

    L = par
    lags = -L:L

    # Smoothed viewing window
    lam = _window(L, win = win, two_sided = false)

    R_pred_smoothed = zeros(Complex,nu,nu,length(0:L))
    for i = 1 : nu
        for j = 1 : nu
            temp = my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
            temp = .5*(temp[L+1:end] + conj(reverse(temp[1:L+1])))
            R_pred_smoothed[i,j,:] = lam .* temp
        end
    end

    # Compute coefficients of spectral factorization of z-spect-pred
    l = PI ? spectfact_matrix_CKMS_pinv(R_pred_smoothed,rtol = rtol) :
             spectfact_matrix_CKMS(R_pred_smoothed)

    l_pad_minus = nfft >= L+1 ? cat(dims = 3,l,zeros(nu,nu,nfft - L - 1)) :
                               l[:,:,1:nfft]

    z_spect_pred_minus_num_fft = fft(l_pad_minus,3)
    z_spect_pred_plus_num_fft = zeros(Complex,nu,nu,nfft)
    for i = 1 : nfft
        z_spect_pred_plus_num_fft[:,:,i] = z_spect_pred_minus_num_fft[:,:,i]'
    end

    # Compute z-cross-spectrum of sigpred
    z_crossspect_sigpred_num_fft = z_crossspect_fft(sig, pred,
                        nfft = nfft, n = n, p = p, win = "Par");

    # This computes the impule response (coefficeints of z) for S_{yx}{S_x^+}^{-1}
    S_sigpred_overS_plus_fft_num = complex(zeros(d,nu,nfft))

    for i = 1 : nfft
        S_sigpred_overS_plus_fft_num[:,:,i] = z_crossspect_sigpred_num_fft[:,:,i]/
                                              z_spect_pred_plus_num_fft[:,:,i]
    end

    S_sigpred_overS_plus_fft_num_fft = ifft(S_sigpred_overS_plus_fft_num,3)

    # Extracts causal part coefficinets of S_{yx}{S_x^+}^{-1}, {S_{yx}{S_x^+}^{-1}}_+
    S_sigpred_overS_plus_fft_plus_num_fft = cat(dims = 3,
                    S_sigpred_overS_plus_fft_num_fft[:,:,1: nffth],
                    zeros(d,nu,nfft - nffth))

    # Computes causal part of S_{yx}/S_x^+, {S_{yx}/S_x^+}_+
    S_sigpred_overS_plus_plus_num_fft = fft(S_sigpred_overS_plus_fft_plus_num_fft,3);

    # Obtain transfer function H by dividing {S_{yx}/S_x^+}_+ by S_x^-

    H_num = complex(zeros(d,nu,nfft))
    for i = 1: nfft
        H_num[:,:,i] = S_sigpred_overS_plus_plus_num_fft[:,:,i]/
                       z_spect_pred_minus_num_fft[:,:,i]
    end

    # Extrct tranferfunction coeffifcients (impulse responce of Weiner filter)
    h_num_raw = ifft(H_num, 3)

    # Truncate
    M_out > nfft && println("M_out > nfft, taking min")
    M = min(M_out, nfft)
    h_num_fft = h_num_raw[:,:,1:M]
end