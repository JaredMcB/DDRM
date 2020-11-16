#####################################################################
"""
Here are just some tools for analyzing the results of model reduction
by the wiener filter.
"""

using Distributions


function plot_WF(h_wf;rl = true, plotter = plot)
    d, nu, M_out = size(h_wf)

    H_wf_plotable = zeros(M_out,d*nu)

    for i=1:d
        for j=1:nu
            plotter(real(h_wf[i,j,:]),label = "El ($i,$j) real")
            rl || plotter(imag(h_wf[i,j,:]),":",label = "El ($i,$j) imag")
        end
    end
    legend()
    xlabel("Lags")
    title("Elements of h_wf over lags")
end


function redmodrun(X, h_wf, Psi;
    steps = size(X,2),
    noise = true,
    noise_dist = MvNormal(I + zeros(size(X,1),size(X,1)))
    )

    d, nu, M_out = size(h_wf)

    bust = 0

    PSI  = complex(zeros(nu,steps))
    X_rm = complex(zeros(d,steps))
    X_rm[:,1:M_out] = X[:,1:M_out]

    # Initialize the PSI sequence
    for i = 1:M_out
        PSI[:,i] = Psi(X_rm[:,i])
    end

    Noise = noise ? rand(noise_dist,steps) : zeros(d,steps)

    for i = M_out+1:steps
        X_rm[:,i] = sum(h_wf[:,:,k]*PSI[:,i-k]
            for k = 1:M_out,dims = 2) + Noise[:,i]
                PSI[:,i] = Psi(X_rm[:,i])
        if isnan(X_rm[1,i])
            bust = i
            break
        end
    end

    if bust != 0
        println("reduced model blewup at step $bust")
    else
        println("reduced model did not blowup")
    end
    X_rm
end


function one_step_pred(X_sig, h_wf, pred)
    # Observe that if we are to have the predictors one step behind the signal
    # then we need to discard the first signal point(after it has been included
    # into the predictors series. This way the series are the same length.
    d, nu, M_out = size(h_wf)
    steps = size(X_sig,2)

    X_hat = complex(zeros(d,steps))
    X_hat[:,1:M_out] = X_sig[:,1:M_out]

    for i = M_out+1:steps
        X_hat[:,i] = sum(h_wf[:,:,k]*pred[:,i-k+1]
            for k = 1:M_out,dims = 2)
        # The '+1' in the pred is important to make sure
        # These the preds and sigs line up. In the
    end
    X_hat
end

function causal_test(X_hat, X_sig, pred;
    lags = -100:10)

    X_err = X_sig - X_hat
    C = my_crosscov(X_err[:],pred[:],-100:10)
end
