module AnalysisToolbox

using Statistics
using FFTW
using LinearAlgebra
using DSP: conv, nextfastfft
using Polynomials
using StatsBase
using SparseArrays



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

function my_crosscor(x::AbstractVector{<:Number},
                     y::AbstractVector{<:Number},
                     lags)
    my_crosscov(x,y,lags)/my_crosscov(x,y,0:0)[1]
end

function my_autocov(x::AbstractVector{<:Number},
                     lags)
    # length(lags) > 1000 ? _crosscov_con(x,x, lags) : _crosscov_dot(x,x, lags)
    _crosscov_con(x,x, lags)
end

function my_autocor(x::AbstractVector{<:Number},
                     lags)
    my_crosscov(x,x,lags)/my_crosscov(x,x,0:0)[1]
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

function smoother_plot(n,p,ty)
    μ = _smoother(n,p; ty)
    plot(-n*p:n*p,μ)
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
    z_crossspect_fft(sig::Array{Complex{Float64},2},pred::Array{Complex{Float64},2};
n = 3, p = 2500, win = "Par")

The output has length nfft = nextfastfft(steps)
"""
function z_crossspect_fft(
    sig,
    pred::Array{T,2} where T <: Number;
    nfft = 0,
    n = 2,
    p = 5,
    ty = "bin")

    ## sig = d x steps, pred = nu x steps
    d, stepsx = size(sig)
    nu, stepsy = size(pred)

    stepsx == stepsy || print("sig and pred are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])
    nfft = nfft == 0 ? nextfastfft(steps) : nfft
    # steps == nfft || println("adjusted no. of steps from $steps to $nfft")

    # blks = ceil(Int, steps/nfft)

    z_spect_mat = zeros(Complex, d, nu, nfft)
    for i = 1 : d
        for j = 1 : nu
            z_spect_mat[i,j,:] = z_crossspect_scalar_ASP(sig[i,:],pred[j,:];
                                                  nfft, n, p,ty)
        end
    end
    z_spect_mat
end

function z_spect_scalar(sig; n = 3, p=100, ty = "ave")
    μ = _smoother(n,p;ty)

    siz = length(sig)
    nfft = nextfastfft(siz)
    sig_pad = [sig; zeros(nfft - siz)]

    peri = abs.(fft(sig_pad)).^2/nfft
    peri_pad = [peri[end - p*n + 1 : end]; peri; peri[1:p*n]]
    z_spect_smoothed = conv(μ,peri_pad)[2n*p+1:2n*p+nfft]
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
    n = 2,
    p = 5,
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

function auto_times(x::AbstractVector{<:Real};plt = false)
    lx = size(x,1)
    L = minimum([lx - 1, 10^6])

    lags = 0:L
    A = real(my_autocov(x,lags))

    end_int = try
                findall(A.<0)[1] - 1
            catch e
                if isa(e, BoundsError)
                    L
                end
            end
    end_exp = Int64(round(end_int/3)) #what if int_int < 2? we get and error

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

function visual_test_ckms(P,l,nfft;semilog = false)
    d  = size(P,1)
    lp = size(P,3)
    ll = size(l,3)
    S_fun(z)    = P[:,:,1] + sum(P[:,:,i]*z^(-i+1) + P[:,:,i]'*z^(i-1) for i = 2:lp)
    S_fun_minus(z) = sum(l[:,:,i]*z^(-i+1) for i = 1:ll)
    S_fun_plus(z) = sum(l[:,:,i]'*z^(i-1) for i = 1:ll)

    Θ = 2π*(0:nfft-1)/nfft
    Z = exp.(im*Θ)
    S = complex(zeros(d,d,nfft))
    S_l = complex(zeros(d,d,nfft))
    for i = 1:nfft
        S[:,:,i] = S_fun(Z[i])
        S_l[:,:,i] = S_fun_minus(Z[i])*S_fun_plus(Z[i])
    end


    for i = 1:d
        for j = i:d
            semilog ? semilogy(Θ,real(S[i,j,:]), label = "S ($i,$j)") :
                      plot(Θ,real(S[i,j,:]), label = "S ($i,$j)")

            semilog ? semilogy(Θ,real(S_l[i,j,:]), label = "S_l ($i,$j)") :
                      plot(Θ,real(S_l[i,j,:]), label = "S_l ($i,$j)")
        end
    end
    legend()
end

function emp_cdf(series;
    bn = 0,
    plt = true)
    l = length(series)
    series = reshape(series,l)
    sort!(series)

    if bn != 0
        bw = (series[end] - series[1])/bn
    else
        bw = 2*iqr(series)/l^(1/3)
        bn = Int64(ceil((series[end] - series[1])/bw))
    end

    b_pts = bw*(0:bn) .+ series[1]

    cdf = zeros(bn+1)
    for i = 1:bn
        cdf[i+1] = sum(series .< b_pts[i+1])
    end
    cdf /= l

    plt ? [cdf,b_pts,plot(b_pts,cdf)] : [cdf,b_pts]
end

function emp_pdf(series;
    bn = 0,
    plt = true)

    cdf, b_pts = emp_cdf(series; bn, plt = false)
    bn = length(cdf) - 1
    bw = b_pts[2] - b_pts[1]

    pdf = zeros(bn)
    b_midpts = zeros(bn)
    for i=1:bn
        pdf[i] = cdf[i+1] - cdf[i]
        b_midpts[i] = (b_pts[i] + b_pts[i+1])/2 # Midpoint as bin location
    end
    pdf /= bw

    plt ? [pdf,b_midpts,plot(b_midpts,pdf,label = "Emp pdf")] :
          [pdf,b_midpts]
end

"""
Title: ARNA_Generator.jl
Author: Jared McBride (Nov 13, 2019)

Given the inputs it produces a realization of an ARMA(p,q) process of length
steps.

The inputs are
  l the list autoregressive coefficients (a_0, a_1, a_2, ... a_p)
  w the list of moving average coefficinets (b_0, b_1, b_2, ... b_q)
  r the variance of the white noise process
The kwargs are
  steps the length of the output time series
  Zeros the set of zeros for the tranfer function. This may be used in place
      w, so it prescribes the moving avearge behavior.
  Poles the set of poles for the tranfer function. This may be used in place
      l, so it prescribes the autoregressive behavior.
  e is a input signal if desired, if e is a generale time series the result
      this filtered by H(z) = (b_0 + b_1*z^(-1) + b_2*z^(-2) + ... + b_q*z^(-q))
                              /(1 + a_1*z^(-1) + a_2*z^(-2) + ... + a_p*z^(-p)).

 To ensure the stability of the AR part the roots of
      (1 + a_1*z^(1) + a_2*z^(2) + ... + a_p*z^(p))
      need to be outside of the unit circle. So, The roots input above are
      actually the roots of
      (1 + a_1*z^(-1) + a_2*z^(-2) + ... + a_p*z^(-p))

The output is just the time series x
"""


function ARMA_gen(; l = [1, -5/4, 3/8],
                    w = [1],
                    r::Float64 = 1.0,
                    steps::Int64 = 10^4,
                    Zeros = [],
                    Poles = [],
                    e = [],
                    discard::Int64 = 10^3)

    p, q = length(l) - 1, length(w) - 1

    if Poles != []
        p = length(Poles)
        P = prod([Polynomial([1]); [Polynomial([1,-z]) for z in Poles]])
        # Produces a poly with roots: Poles.^(-1)
        l = coeffs(P);
    end

    if Zeros != []
        q = length(Zeros)
        Q = prod([Polynomial([1]); [Polynomial([1,-z]) for z in Zeros]])
        w = coeffs(Q);
    end
    steps_tot = steps + discard

    x = complex(zeros(steps_tot));
    if e == []
        e = sqrt(r) * randn(steps_tot);
    end

    pvq = maximum([p q])
    if length(l) == 0
        for i = pvq + 1 : steps_tot
            x[i] = dot(reverse(w),e[i - q:i])
        end
    else
        for i = pvq + 1 : steps_tot
            x[i] = dot(-reverse(l)[1:p],x[i - p:i-1]) + dot(reverse(w),e[i - q:i])
        end
    end
    x[discard + 1:end]
end
end # module
