using FFTW
using LinearAlgebra: dot, I
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
    Nex = 1024,
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

    h_wf = rl ? real(vector_wiener_filter_fft(pred, sig, M_out, Nex = Nex,
                        PI = PI, rtol = rtol)) :
                     vector_wiener_filter_fft(pred, sig, M_out, Nex = Nex,
                        PI = PI, rtol = rtol)
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

function spectfact_matrix_CKMS(P)
    d = size(P)[1];
    m = size(P)[3] - 1

    NN = reverse(P[:,:,2:end],dims = 3)
    Re = Rr = p0 = P[:,:,1]

    F = [[zeros(d,d*(m-1)); I] zeros(d*m,d)]
    h = [zeros(d,d*(m-1)) I]

    K = complex(zeros(d*m,d))
    for i = 0 : m-1
        K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
    end
    L = K

    for i = 1:200
        hL = h*L; FL = F*L

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
    re = Re

    sqrt_re = sqrt(re)

    l = complex(zeros(d,d,m+1))
    l[:,:,1] = sqrt_re;
    for i = m-1:-1:0
        l[:,:,m-i+1] = k[d*i + 1: d*(i+1),:]*sqrt_re
    end
    l
end

function spectfact_matrix_CKMS_pinv(P, rtol = 1e-6)
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

    for i = 1:200
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
                      lags::UnitRange{Int})
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
                      lags::UnitRange{Int})
    L = min(length(x),length(y))
    m = length(lags)

    zx = x .- mean(x)
    zy = y .- mean(y)

    C = complex(zeros(m))
    for k = 1:m
        l = lags[k]
        C[k] = ( l >= 0 ? dot(zx[1+l : L],zy[1 : L-l]) : dot(zx[1 : L+l],zy[1-l : L]))/L
    end
    C
end

function my_crosscov(x::AbstractVector{<:Number},
                     y::AbstractVector{<:Number},
                     lags::UnitRange{Int})
    length(lags) > 1000 ? _crosscov_con(x,y, lags) : _crosscov_dot(x,y, lags)
end

function my_crosscor(x::AbstractVector{<:Number},
                     y::AbstractVector{<:Number},
                     lags::UnitRange{Int})
    my_crosscov(x,y,lags)/my_crosscov(x,y,0:0)[1]
end

function z_crossspect_fft(
    sig::Array{Complex{Float64},2},
    pred::Array{Complex{Float64},2};
    L = 50,
    Nex = 2^10,
    win = "Par")

    ## sig = d x steps, pred = nu x steps
    d, stepsx = size(sig)
    nu, stepsy = size(pred)

    Nexh = Int(floor(Nex/2))
    lags = -L:L;

    stepsx == stepsy || print("sig and pred are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])

    # Smoothed viewing window
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

    C_smoothed = complex(zeros(d,nu,length(lags)))
    for i = 1 : d
        for j = 1 : nu
            C_smoothed[i,j,:] = Lam .* my_crosscov(sig[i,1:steps],pred[j,1:steps],lags)
        end
    end

    ## C_smoothed = d x nu x 2L+1

    ## Pad with zeros in preparation for fft
    C_padded = cat(dims = 3, zeros(d,nu,Nex - Nexh - L), C_smoothed, zeros(d,nu,Nexh - L - 1))
    C = fftshift(C_padded,3)

    z_crossspect_num_fft = fft(C,3);
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

function z_crsspect_scalar(sig,pred; n = 3, p=100, ty = "ave")
    μ = _smoother(n,p,ty = "ave")

    l_sig = length(sig)
    l_pred = length(pred)
    l_sig == l_pred || println("sizes must be the same, taking min and truncating")
    l = min(l_sig,l_pred)

    nfft = nextfastfft(l)
    sig_pad = [sig[1:l]; zeros(nfft - l)]
    pred_pad = [pred[1:l]; zeros(nfft - l)]

    fftsig = fft(sig_pad)
    fftpred = conj(fft(pred_pad))

    peri = fftsig .* fftpred / nfft
    peri_pad = [peri[end - p*n + 1 : end]; peri; peri[1:p*n]]
    z_crsspect_smoothed = conv(μ,peri_pad)[2n*p:2n*p+l-1]
end



"""
    vector_wiener_filter_fft

"""
function vector_wiener_filter_fft(pred::Array{Complex{Float64},2}, sig,
    M_out = 20;
    par::Int64 = 55,
    Nex::Int64 = 1024,
    win = "Par",
    PI = true,
    rtol = 1e-6
    )

    d, stepsy = size(sig)
    nu, stepsx = size(pred)

    stepsx == stepsy || print("X and Y are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])

    Nexh = Int(floor(Nex/2))

    L = par
    lags = 0:L;

    # Smoothed viewing window
    lam = _window(L, win = win, two_sided = false)

    R_pred_smoothed = complex(zeros(nu,nu,length(lags)))
    for i = 1 : nu
        for j = 1 : nu
            R_pred_smoothed[i,j,:] = lam .* my_crosscov(pred[i,1:steps],pred[j,1:steps],lags)
        end
    end

    # Compute coefficients of spectral factorization of z-spect-pred
    l = PI ? spectfact_matrix_CKMS_pinv(R_pred_smoothed,rtol = rtol) :
             spectfact_matrix_CKMS(R_pred_smoothed)

    l_pad_minus = Nex >= L+1 ? cat(dims = 3,l,zeros(nu,nu,Nex - L - 1)) :
                               l[:,:,1:Nex]

    z_spect_pred_minus_num_fft = fft(l_pad_minus,3)
    z_spect_pred_plus_num_fft =complex(zeros(nu,nu,Nex))
    for i = 1 : Nex
        z_spect_pred_plus_num_fft[:,:,i] = z_spect_pred_minus_num_fft[:,:,i]'
    end

    # Compute z-cross-spectrum of sigpred
    z_crossspect_sigpred_num_fft = z_crossspect_fft(sig, pred, L = par, Nex = Nex, win = "Par");

    # This computes the impule response (coefficeints of z) for S_{yx}{S_x^+}^{-1}
    S_sigpred_overS_plus_fft_num = complex(zeros(d,nu,Nex))

    for i = 1: Nex
        S_sigpred_overS_plus_fft_num[:,:,i] = z_crossspect_sigpred_num_fft[:,:,i]/
                                              z_spect_pred_plus_num_fft[:,:,i]
    end

    S_sigpred_overS_plus_fft_num_fft = ifft(S_sigpred_overS_plus_fft_num,3)

    # Extracts causal part coefficinets of S_{yx}{S_x^+}^{-1}, {S_{yx}{S_x^+}^{-1}}_+
    S_sigpred_overS_plus_fft_plus_num_fft = cat(dims = 3,
                    S_sigpred_overS_plus_fft_num_fft[:,:,1: Nexh],
                    zeros(d,nu,Nex - Nexh))

    # Computes causal part of S_{yx}/S_x^+, {S_{yx}/S_x^+}_+
    S_sigpred_overS_plus_plus_num_fft = fft(S_sigpred_overS_plus_fft_plus_num_fft,3);

    # Obtain transfer function H by dividing {S_{yx}/S_x^+}_+ by S_x^-

    H_num = complex(zeros(d,nu,Nex))
    for i = 1: Nex
        H_num[:,:,i] = S_sigpred_overS_plus_plus_num_fft[:,:,i]/
                       z_spect_pred_minus_num_fft[:,:,i]
    end

    # Extrct tranferfunction coeffifcients (impulse responce of Weiner filter)
    h_num_raw = ifft(H_num,3)

    # Truncate
    M_out > Nex && println("M_out > Nex, taking min")
    M = min(M_out, Nex)
    h_num_fft = h_num_raw[:,:,1:M]
end
