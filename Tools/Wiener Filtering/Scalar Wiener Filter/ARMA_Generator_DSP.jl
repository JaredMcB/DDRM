##############################################################################
# Title: ARNA_Generator.jl
# Author: Jared McBride (Nov 13, 2019)
#
# Given the inputs it produces a realization of an ARMA(p,q) process of length
# steps.
#
# The inputs are
#   l the list autoregressive coefficients (a_0, a_1, a_2, ... a_p)
#   w the list of moving average coefficinets (b_0, b_1, b_2, ... b_q)
#   r the variance of the white noise process
# The kwargs are
#   steps the length of the output time series
#   Zeros the set of zeros for the tranfer function. This may be used in place
#       w, so it prescribes the moving avearge behavior.
#   Poles the set of poles for the tranfer function. This may be used in place
#       l, so it prescribes the autoregressive behavior.
#   e is a input signal if desired, if e is a generale time series the result
#       this filtered by H(z) = (b_0 + b_1*z^(-1) + b_2*z^(-2) + ... + b_q*z^(-q))
#                               /(1 + a_1*z^(-1) + a_2*z^(-2) + ... + a_p*z^(-p)).
#
#  To ensure the stability of the AR part the roots of
#       (1 + a_1*z^(1) + a_2*z^(2) + ... + a_p*z^(p))
#       need to be outside of the unit circle. So, The roots input above are
#       actually the roots of
#       (1 + a_1*z^(-1) + a_2*z^(-2) + ... + a_p*z^(-p))
#
# The output is just the time series x
#
##############################################################################

using Polynomials
using LinearAlgebra
using StatsBase

function ARMA_gen(  l = [1, -5/4, 3/8],
                    w = [1];
                    r::Float64 = 1.0,
                    steps::Int64 = 10^4,
                    Zeros = [],
                    Poles = [],
                    e = [],
                    discard::Int64 = 10^3,
                    out_poly = false)

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
    out_poly ? [x[discard + 1:end], P, Q] : x[discard + 1:end]
end

function z_spect(x,L; win = "Bart")
    lags = 0:L;
    A = autocov(x,lags)
    A[1] = A[1]/2

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

    z_spect_num(z) = sum([lam[i+1]*(A[i+1]'*z^(-i) +
                        A[i+1]*z^(i)) for i = 0 : L])
end

dB(s) = 10*log.(s)./log(10)
