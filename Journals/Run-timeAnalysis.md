# Full run-time and Memory Usage Analysis of The Wiener Filtering Algorithm

The usual call when we want the Wiener filter (WF) is
```julia
h_wf = mr.get_wf(signal, Psi; M_out)
```
The signal here is the true data, 𝑑 by 𝑁. the function `Psi` takes in vectors of length 𝑑 and returns a vector of length ν. `M_out`(= 𝑚) specifies the number of coefficients of the resulting WF approximation, so that `h_wf` has dimensions 𝑑 × ν × 𝑚.

### `get_wf`
The algorithm can be organized as follows:
1. Offset `signal` and call it `sig`.
2. Extract parameters `d`, `steps`, and `nu`.
3. Generate predictors (`pred`)
   1. allocate space
   2. populate array with images of `sig` by `Psi`.
4. Call `vector_wiener_filter_fft`
5. Sort out output based on flags.

Now, we look at the Wiener filtering algorithm itself:

### `vector_wiener_filter_fft`
This algorithm is organized as follows
1. Get parameters: `d`, `nu`, `steps`, `nfft`, `nffth`, and `L`.
2. Compute smoothed autocovariance sequence (`R_pred_smoothed`) of `pred` of length `L+1`
   1. This calls the function `matrix_autocov_seq`
3. Compute coefficients `l` of spectral factorization of `R_pred_smoothed`.
   1. This calls the function `spectfact_matrix_CKMS`. This is perhaps the most costly part of the      whole process.
4. Pad `L` and compute the point-wise approximations of 𝑆ₚᵣ⁺(𝑧) (`z_spect_pred_plus_num_fft`) and 𝑆ₚᵣ⁻(𝑧) (`z_spect_pred_minus_num_fft`),
   1. The padding is done useing the `cat` function.
   2. The pointwise approximations are  formed using `fft`
5. Compute z-cross-spectrum `z_crossspect_sigpred_num_fft` of the `sig` and `pred`
   This calls `at.z_crossspect_fft_old`
6. Point-wise divide  `z_crossspect_sigpred_num_fft` by `z_spect_pred_plus_num_fft` to form `S_sigpred_overS_plus_fft_num`
7. Take `ifft` of `S_sigpred_overS_plus_fft_num` to form `S_sigpred_overS_plus_fft_num_fft`
8. Forget half of `S_sigpred_overS_plus_fft_num_fft` and replace with zeros to form `S_sigpred_overS_plus_fft_plus_num_fft`
9. Take `fft` of `S_sigpred_overS_plus_fft_plus_num_fft` to form `S_sigpred_overS_plus_plus_num_fft`
10. point-wise devide `S_sigpred_overS_plus_plus_num_fft` by `z_spect_pred_minus_num_fft` to arrive at `H_num`
11. take 'fft' of `H_num` to get the coefficients of the wiener filter.

```julia
function vector_wiener_filter_fft(
    sig,
    pred::Array{T,2} where T <: Number;
    M_out = 20,
    par::Int64 = 1500,
    nfft = 0,
    win = "Par",
    n = 3,
    p = 1500,
    ty = "bin",
    xspec_est = "old",
    PI = true,
    rtol = 1e-6,
    info = false
    )

    d, stepsy = size(sig)
    nu, stepsx = size(pred)

    stepsx == stepsy || print("X and Y are not the same length. Taking min.")
    steps = minimum([stepsx stepsy])
    nfft = nfft == 0 ? nextfastfft(steps) : nfft
    nffth = nfft ÷ 2
    L = par

    R_pred_smoothed = matrix_autocov_seq(pred; L, steps, nu, win)

    # Compute coefficients of spectral factorization of z-spect-pred
    l = @time PI ? spectfact_matrix_CKMS_pinv(R_pred_smoothed,rtol = rtol) :
             spectfact_matrix_CKMS(R_pred_smoothed)

    l_pad_minus = nfft >= L+1 ? cat(dims = 3,l,zeros(nu,nu,nfft - L - 1)) :
                               @view l[:,:,1:nfft]

    z_spect_pred_minus_num_fft = fft(l_pad_minus,3)
    z_spect_pred_plus_num_fft = complex(zeros(nu,nu,nfft))
    for i = 1 : nfft
        z_spect_pred_plus_num_fft[:,:,i] = @view z_spect_pred_minus_num_fft[:,:,i]'
    end

    # Compute z-cross-spectrum of sigpred
    z_crossspect_sigpred_num_fft = xspec_est == "SP" ? at.z_crossspect_fft(sig, pred;
                        nfft, n, p, ty) : at.z_crossspect_fft_old(sig, pred; L, Nex = nfft);

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

    if info
        h_num_fft = [h_num_raw[:,:,1:M],
                 z_crossspect_sigpred_num_fft,
                 z_spect_pred_minus_num_fft,
                 z_spect_pred_plus_num_fft,
                 S_sigpred_overS_plus_fft_num,
                 S_sigpred_overS_plus_plus_num_fft,
                 H_num] ###
    else
        h_num_fft = h_num_raw[:,:,1:M]
    end
    h_num_fft
end
```
