module WFMR_bs

using DSP: conv # Solely for Psi


at = include("AnalysisToolbox.jl")

"""
    `get_wf` provides the causal wiener filter that best approximates `signal`
in terms of `Psi(signal[:,i-1])`. It has built in the one step prediction. First,
it generates the predictor series. Then it forms a design matrix which houses the
time-lagged predictors. We then use backslash to solve for the linear least squares 
estimate of the filter.

Here is an example of it's implementations. `sig` is the signal ypu are using `pred`
to approximate. `M_out = 60` is the number of coefficients we get out.

M_out = 60
h_wf = mrb.get_wf_bs(sig1, pred1; M_out);

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




"""
sig and pred should 2D arrays
"""

function get_wf_bs(sig::Array{<: Number,2},pred::Array{<: Number,2}; M_out)
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


###########################
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
