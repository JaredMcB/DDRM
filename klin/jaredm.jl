## This is taken from DDMR/Tools/AnalysisToolbox.jl

module jaredm

using Statistics
using FFTW
using LinearAlgebra
using DSP: conv, nextfastfft
using Polynomials
using StatsBase

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
    # length(lags) > 1000 ? _crosscov_con(x,y, lags) : _crosscov_dot(x,y, lags)
    _crosscov_con(x,y, lags)
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
        μ_c = ones(2n*p+1)/(2n*p+1)
    end
    round(sum(μ_c);digits = 5) == 1.0 || println("bad smoother")
    μ_c
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
    ty = "bin")

    ## sig = d x steps, pred = nu x steps
    d, stepsx = size(sig)
    nu, stepsy = size(pred)

    stepsx == stepsy || print("sig and pred are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])
    nfft = nfft == 0 ? nextfastfft(steps) : nfft
    # steps == nfft || println("adjusted no. of steps from $steps to $nfft")
    steps = nfft

    z_spect_mat = zeros(Complex, d, nu, nfft)
    for i = 1 : d
        for j = 1 : nu
            z_spect_mat[i,j,:] = z_crossspect_scalar_ASP(sig[i,:],pred[j,:];
                                                  nfft, n, p,ty)
        end
    end
    z_spect_mat
end

"""
z_crsspect_scalar has output of size nfft
"""
function z_crossspect_scalar(sig,pred; nfft = 0, n = 3, p=100, ty = "ave")
    μ = _smoother(n,p;ty)

    # Of cousre we need these to be mean zero
    sig .-= mean(sig)
    pred .-= mean(pred)

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
    z_crsspect_smoothed = conv(μ,peri_pad)[2n*p+1:2n*p+nfft]
end

function z_crossspect_scalar_ASP(
    sig,
    pred;
    nfft = 2^10, # The length of each subseries
    n = 3,
    p = 10,
    ty = "bin")

    # Of cousre we need these to be mean zero
    sig .-= mean(sig)
    pred .-= mean(pred)

    # Check length of series
    l_sig = length(sig)
    l_pred = length(pred)
    l_sig == l_pred || println("sizes must be the same, taking min and truncating")
    l = min(l_sig,l_pred)

    # The total nuber of subseries
    R = floor(Int,l/nfft)
    # Computation of the average periodogram
    aperi = complex(zeros(nfft))
    for r = 1:R
        fftsig = fft(sig[(r-1)*nfft+1:r*nfft])
        fftpred = conj(fft(pred[(r-1)*nfft+1:r*nfft]))
        aperi .+= fftsig .* fftpred
    end
    aperi ./= nfft*R

    # Smoothing it too.
    if ty != "none"
        aperi_pad = [aperi[end - p*n + 1 : end]; aperi; aperi[1:p*n]]
        μ = _smoother(n,p; ty)
        aperi = conv(μ,aperi_pad)[2n*p+1:2n*p+nfft]
    end
    aperi
end

rowmatrix(x) = reshape(x,1,length(x))
z_crossspect_dm(sig,pred; flags...) = z_crossspect_fft_old(rowmatrix(sig), rowmatrix(pred); flags...)[1:end]

function z_crossspect_fft_old(
    sig,
    pred;
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

end#module
