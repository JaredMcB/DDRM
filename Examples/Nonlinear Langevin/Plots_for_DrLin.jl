fig, axs = subplots(4,2, sharex = "all", sharey ="row")
# xspect
axs[1].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_dm[1,1,:]),label = "dm1")
axs[1].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_dm[1,2,:]),label = "dm2",linestyle = "dashed")
axs[1].set_title("Direct Method (real part)")
axs[1].set_ylabel("S_YX")
# axs[1].axis([.1, 6.2,-1e-1,1e-1])
axs[1].grid("on")
axs[1].legend()

axs[2].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1")
axs[2].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,2,:]),label = "dm2",linestyle = "dashed")
axs[2].set_ylabel("S_YX/S_X^+")
axs[2].grid("on")
axs[2].legend()

axs[3].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1")
axs[3].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,2,:]),label = "dm2",linestyle = "dashed")
axs[3].set_ylabel("{S_YX/S_X^+}_+")
axs[3].grid("on")
axs[3].legend()

axs[4].plot(2pi*(0:nfft-1)/nfft,real(H_num_dm[1,1,:]),label = "dm1")
axs[4].plot(2pi*(0:nfft-1)/nfft,real(H_num_dm[1,2,:]),label = "dm2",linestyle = "dashed")
axs[4].set_ylabel("{S_YX/S_X^+}_+/S_X^-")
axs[4].grid("on")
axs[4].set_xlabel("frequencies")
axs[4].legend()

axs[5].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_sp[1,1,:]),label = "sp1")
axs[5].semilogy(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_sp[1,2,:]),label = "sp2",linestyle = "dashed")
axs[5].set_title("Periodogram (real part)")
# axs[5].axis([.1, 6.2,-1e-1,1e-1])
axs[5].grid("on")
axs[5].legend()

axs[6].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1")
axs[6].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,2,:]),label = "sp2",linestyle = "dashed")
axs[6].grid("on")
axs[6].legend()

axs[7].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1")
axs[7].semilogy(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,2,:]),label = "sp2",linestyle = "dashed")
axs[7].grid("on")
axs[7].legend()

axs[8].plot(2pi*(0:nfft-1)/nfft,real(H_num_sp[1,1,:]),label = "sp1")
axs[8].plot(2pi*(0:nfft-1)/nfft,real(H_num_sp[1,2,:]),label = "sp2",linestyle = "dashed")
axs[8].grid("on")
axs[8].set_xlabel("frequencies")
axs[8].legend()




SVD1s = [svd(S_pred_plus[:,:,n]).S[1] for n = 1:nfft]
SVD2s = [svd(S_pred_plus[:,:,n]).S[2] for n = 1:nfft]

fig, (ax1, ax2) = subplots(1,2)
ax1.plot(2pi*(0:nfft-1)/nfft,SVD1s,label = "1")
ax1.plot(2pi*(0:nfft-1)/nfft,SVD2s,label = "2",linestyle = "dashed")
ax1.set_title("Singular values of S_X^+ plot")
ax1.set_xlabel("frequency")
ax1.set_ylabel("fisrt and second singular values")
ax1.legend()
ax2.semilogy(2pi*(0:nfft-1)/nfft,SVD1s,label = "1")
ax2.semilogy(2pi*(0:nfft-1)/nfft,SVD2s,label = "2",linestyle = "dashed")
ax2.set_title("Singular values of S_X^+ semilogy")
ax2.legend()
ax2.set_xlabel("frequency")



semilogx(SVD1s ./SVD2s)
