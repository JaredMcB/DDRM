module WFMR_lasso
using DSP: conv # Solely for Psi

using GLMNet
import StatsBase: zscore!

zscore!(X) = zscore!(X, mean_and_std(X)[1], mean_and_std(X)[2])

at = include("AnalysisToolbox.jl")

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

function get_pred(signal, Psi)
    d, steps = size(signal)
    nu = size(Psi(zeros(d,1)),1)

    pred = zeros(ComplexF64, nu, steps)
    for n = 1:steps
        pred[:,n] = Psi(@view signal[:,n])
    end
    pred
end

function Get_PRED(sig,pred)
    d, stepsy = size(sig)
    nu, stepsx = size(pred)
    stepsx == stepsy || print("X and Y are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])

    sig = Array(transpose(sig))
    pred = Array(transpose(pred))

    Pred = zeros(ComplexF64,steps-M_out+1,M_out*nu)
    for m = 1:M_out
        Pred[:,nu*(M_out-m)+1:nu*(M_out-m+1)] = @view pred[m:(steps + m - M_out),:]
    end

    PRED = [real(Pred) -imag(Pred); imag(Pred) real(Pred)]
    SIG = [real(sig[M_out:steps,:]); imag(sig[M_out:steps,:])]
    SIG, PRED
end

function mkfilter(beta; l = length(beta) ÷ 2, d = 1, M_out, nu = l ÷ M_out)
    l = length(beta) ÷ 2
    h = zeros(ComplexF64,l,1)
    for i = 1:l
        h[i] = complex(beta[i],beta[i+l])
    end
    h_wfls = zeros(ComplexF64,d,nu,M_out)
    for m = 1:M_out
        h_wfls[:,:,m] = (@view h[(nu*(m-1) + 1):nu*m,:])'
    end
    h_wfls 
end

function get_wf_ls(signal,Psi::Function; M_out)
    d, steps = size(signal)
    nu = size(Psi(zeros(d,1)),1)

    sig = @view signal[:,2:steps]   # sig is now one a head of signal
    steps -= 1                      # this makes steps the length of sig

    pred = zeros(ComplexF64, nu, steps)
    for n = 1:steps
    pred[:,n] = Psi(@view signal[:,n])
    end # pred is now even with signal and therefore one step
    # behind sig. I.e. pred[:,n] = Psi(sig[:,n-1])
    # which is what we want so as to ensure the reduced
    # model can run explicitly.

    get_wf_ls(sig,pred; M_out)
end

function get_wf_ls(sig,pred; M_out, lambda = Array(1e-7*(1:500:30000)))
    d, stepsy = size(sig)
    nu, stepsx = size(pred)
    stepsx == stepsy || print("X and Y are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])

    SIG, PRED = Get_PRED(sig,pred; M_out)

    h = zeros(ComplexF64,M_out*nu,d)
    for i = 1:d
        cv = glmnetcv(PRED,SIG[:,i])
        locbetamin = argmin(cv.stdloss)
        h_temp = cv.path.betas[:,locbetamin]
        for j = 1:M_out*nu
            h[j,i] = complex(h_temp[j],h_temp[j+M_out*nu])
        end
    end

    h_wfls = zeros(ComplexF64,d,nu,M_out)
    for m = 1:M_out
        h_wfls[:,:,m] = (@view h[(nu*(m-1) + 1):nu*m,:])'
    end
    h_wfls
end

function get_Psi_2017(h)
    # Build PSI
    function InvBurgRK4_1step(x)
       lx = length(x)
       function F(x)
           𝑥 = [conj(@view x[lx:-1:1]) ;0; x]
           conv(𝑥,𝑥)[2*lx+2:3*lx+1]
       end
       k1 = F(x)
       k2 = F(x .+ h*k1/2)
       k3 = F(x .+ h*k2/2)
       k4 = F(x .+ h*k3)
       A = @. x + h/6*(k1 + 2k2 + 2k3 + k4)
    end

    function Inertialman_part(x)
       lx = length(x)
       𝑥(j) = ( j <= lx ? x[j] : im*sum(x[l]*x[j-l] for l = j-lx:lx) )

       L = complex(zeros(lx^2))
       for j = 1:lx
          for k = 1:lx
             L[ (j-1)*lx+k] = 𝑥(j+lx)*𝑥(j+lx-k)
          end
       end
       L
    end

    Psi(x) = [x; InvBurgRK4_1step(x); Inertialman_part(x)]
    return Psi
end

end #Module
