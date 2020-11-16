z_crossspect_sigpred_num_fft_dm = h_wf_packs[2]
z_crossspect_sigpred_num_fft_sp = h_wf_packs[9]

S_sigpred_overS_plus_fft_num_dm = h_wf_packs[5]
S_sigpred_overS_plus_fft_num_sp = h_wf_packs[12]

S_sigpred_overS_plus_plus_num_fft_dm = h_wf_packs[6]
S_sigpred_overS_plus_plus_num_fft_sp = h_wf_packs[13]

H_num_dm = h_wf_packs[7]
H_num_sp = h_wf_packs[14]

fig, axs = subplots(4,2, sharey ="row")
# xspect
axs[1].semilogx(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_sp[1,1,:]),label = "sp1",".-")
axs[1].semilogx(2pi*(0:nfft-1)/nfft,real(z_crossspect_sigpred_num_fft_dm[1,1,:]),label = "dm1",linestyle = "dashed",".-")
axs[1].set_title("First component DM and Periodogram overlayed (real part)")
axs[1].set_ylabel("S_YX")
# axs[1].axis([.1, 6.2,-1e-1,1e-1])
axs[1].grid("on")
axs[1].legend()

axs[2].semilogx(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1",".-")
axs[2].semilogx(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1",linestyle = "dashed",".-")
axs[2].set_ylabel("S_YX/S_X^+")
axs[2].grid("on")
axs[2].legend()

axs[3].semilogx(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1",".-")
axs[3].semilogx(2pi*(0:nfft-1)/nfft,real(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1",linestyle = "dashed",".-")
axs[3].set_ylabel("{S_YX/S_X^+}_+")
axs[3].grid("on")
axs[3].legend()

axs[4].plot(2pi*(0:nfft-1)/nfft,real(H_num_sp[1,1,:]),label = "sp1")
axs[4].plot(2pi*(0:nfft-1)/nfft,real(H_num_dm[1,1,:]),label = "dm1",linestyle = "dashed")
axs[4].set_ylabel("{S_YX/S_X^+}_+/S_X^-")
axs[4].grid("on")
axs[4].set_xlabel("frequencies")
axs[4].legend()

axs[5].semilogx(2pi*(0:nfft-1)/nfft,imag(z_crossspect_sigpred_num_fft_sp[1,1,:]),label = "sp1",".-")
axs[5].semilogx(2pi*(0:nfft-1)/nfft,imag(z_crossspect_sigpred_num_fft_dm[1,1,:]),label = "dm1",linestyle = "dashed",".-")
axs[5].set_title("First component DM and Periodogram overlayed (Imag part)")
# axs[1].axis([.1, 6.2,-1e-1,1e-1])
axs[5].grid("on")
axs[5].legend()

axs[6].semilogx(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_sp[1,1,:]),label = "sp1",".-")
axs[6].semilogx(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_fft_num_dm[1,1,:]),label = "dm1",linestyle = "dashed",".-")
axs[6].grid("on")
axs[6].legend()

axs[7].semilogx(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_sp[1,1,:]),label = "sp1",".-")
axs[7].semilogx(2pi*(0:nfft-1)/nfft,imag(S_sigpred_overS_plus_plus_num_fft_dm[1,1,:]),label = "dm1",linestyle = "dashed",".-")
axs[7].grid("on")
axs[7].legend()

axs[8].plot(2pi*(0:nfft-1)/nfft,imag(H_num_sp[1,1,:]),label = "sp1")
axs[8].plot(2pi*(0:nfft-1)/nfft,imag(H_num_dm[1,1,:]),label = "dm1",linestyle = "dashed")
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
ax2.semilogx(2pi*(0:nfft-1)/nfft,SVD1s,label = "1")
ax2.semilogx(2pi*(0:nfft-1)/nfft,SVD2s,label = "2",linestyle = "dashed")
ax2.set_title("Singular values of S_X^+ semilogx")
ax2.legend()
ax2.set_xlabel("frequency")



semilogx(SVD1s ./SVD2s)

A = 0.0:4.7936899621426287e-5:10^-2.5

plot(A,real(z_crossspect_sigpred_num_fft_sp[1,1,1:66]),label = "sp1",".-")
plot(A,real(z_crossspect_sigpred_num_fft_dm[1,1,1:66]),label = "dm1",linestyle = "dashed",".-")

sum(real(z_crossspect_sigpred_num_fft_sp[1,1,1:66]))
