using FFTW

Nex = 2^10
Theta = 2*pi*(0 : Nex-1)/(Nex - 1)
Z = exp.(im*Theta)

S(z) = ([1 2; 3 4] + [1 0; 0 2]*z^(-1) + [1 6; 0 2]*z)

D1 = complex(zeros(2,2,1024))
for i = 1:1024
    D1[:,:,i] = S(Z[i])
end

D2 = complex(zeros(2,2,1024))
for i = 1:1024
    D2[:,:,i] = S(Z[i])*Z[i]
end



coef1 = ifft(D2,3)
coef2 = circshift(ifft(D1,3),(0,0,-1))
