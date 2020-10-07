z_spect_pred = z_crossspect_fft(pred,
                    pred; nfft, n = 0, p = 0, win,ty = "Ave");

n = 3
p = 2000
ty = "ave"
μ = _smoother(n,p;ty)
sum(μ)
# μ

# z_spect_pred = z_crossspect_fft(pred[1:1,:],
#                     pred[1:1,:]; nfft, n, p, win,ty);
#
# plot(real(z_spect_pred[1,1,:]))

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

z_spect_pred = z_crsspect_smoothed

# conv(μ,peri_pad)

# peri_pad - z_crsspect_smoothed

semilogy(real(z_spect_pred))

semilogy(peri)
