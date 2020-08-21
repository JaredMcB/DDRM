sig1 = sig[1,:]; pred1 = pred[1,:]


μ = _smoother(n,p,ty = "ave")

l_sig1 = length(sig1)
l_pred1 = length(pred1)
l_sig1 == l_pred1 || println("sizes must be the same, taking min and truncating")
l = min(l_sig1,l_pred1)

nfft = nfft == 0 ? nfft = nextfastfft(l) : nfft

nfft == l || println("adjusted size from $l to $nfft")
sig1_pad = l < nfft ? [sig1[1:l]; zeros(nfft - l)] : sig1[1:nfft]
pred1_pad = l < nfft ? [pred1[1:l]; zeros(nfft - l)] : pred1[1:nfft]

fftsig1 = fft(sig1_pad)
fftpred1 = conj(fft(pred1_pad))

peri = fftsig1 .* fftpred1 / nfft
peri_pad = [peri[end - p*n + 1 : end]; peri; peri[1:p*n]]
z_crsspect_smoothed = conv(μ,peri_pad)[n*p:n*p+nfft-1]

z_crsspect_scalar(sig[1,:],pred[1,:],
                                      nfft = nfft, n = n, p = p)
