using LinearAlgebra
using FFTW

function Scalar_CKMS_f(s_ana; N = 100)
    # N is the number of point used in the fourier transorm it is the number of
    # terms of the approximating lauren polynomial.

    S = [s_ana(exp(im*2*pi*j/N)) for j = 0 : N-1];

    S_fft = fft(S)/N;

    Ne = N - Int64(floor(N/2)); # This is to get only the casual coefficients

    S_fft_1 = S_fft[1:Ne]; # These are only the casual coefficeints

    m = Ne - 1
    F = [[zeros(1,m-1); I] zeros(m)]
    h = [zeros(1,m-1) 1]
    NN = reverse(S_fft_1[2:end]);

    K = reshape(NN,m,1)
    L = reshape(NN,m,1)
    Re = reshape([S_fft_1[1]],1,1)
    Rr = reshape([S_fft_1[1]],1,1)
    for i = 1:100
        hL = h*L
        FL = F*L

        K_new = K - FL/Rr*hL'
        L_new = FL - K/Re*hL
        Re_new = Re - hL/Rr*hL'
        Rr_new = Rr - hL'/Re*hL

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end
    k = K/Re
    re = Re[1]

    l = sqrt(re)*[k[n] for n = m:-1:1]

    s_plus_num(z) = sqrt(re) + sum([l[n]*z^(-n) for n = 1:m]);

    return s_plus_num, l
end

function Scalar_CKMS_c(R)

    m = length(R) - 1
    F = [[zeros(1,m-1); I] zeros(m)]
    h = [zeros(1,m-1) 1]
    NN = reverse(R[2:end]);

    K = reshape(NN,m,1)
    L = reshape(NN,m,1)
    Re = reshape([R[1]],1,1)
    Rr = reshape([R[1]],1,1)
    for i = 1:100
        hL = h*L
        FL = F*L

        K_new = K - FL/Rr*hL'
        L_new = FL - K/Re*hL
        Re_new = Re - hL/Rr*hL'
        Rr_new = Rr - hL'/Re*hL

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end
    k = K/Re
    re = Re[1]

    l = [sqrt(abs(re))]
    [l; sqrt(abs(re))*[k[n] for n = m:-1:1]]
end
