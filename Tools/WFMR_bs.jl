module WFMR_bs


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

function get_wf_bs(signal,Psi::Function; M_out)
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

    get_wf_bs(sig,pred; M_out)
end

function get_wf_bs(sig,pred; M_out)
    d, stepsy = size(sig)
    nu, stepsx = size(pred)
    stepsx == stepsy || print("X and Y are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])

    sig = Array(transpose(sig))
    pred = Array(transpose(pred))

    PRED = zeros(ComplexF64,steps-M_out+1,M_out*nu)
    for m = 1:M_out
        PRED[:,nu*(M_out-m)+1:nu*(M_out-m+1)] = @view pred[m:(steps + m - M_out),:]
    end

    h = PRED \ sig[M_out:steps,:]

    h_wfbs = zeros(ComplexF64,d,nu,M_out)
    for m = 1:M_out
        h_wfbs[:,:,m] = (@view h[(nu*(m-1) + 1):nu*m,:])'
    end
    h_wfbs
end

end #Module
