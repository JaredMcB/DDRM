using FFTW
using LinearAlgebra
using DSP: conv, nextfastfft
using Polynomials
using StatsBase
using SparseArrays


function spectfact_matrix_CKMS_SC(P; ϵ = 1e-10)
    d = size(P)[1];
    m = size(P)[3] - 1

    NN = reverse(P[:,:,2:end],dims = 3)
    Re = Rr = p0 = P[:,:,1]

    F = sparse([[zeros(d,d*(m-1)); I] zeros(d*m,d)])
    h = sparse([zeros(d,d*(m-1)) I])

    K = complex(zeros(d*m,d))
    for i = 0 : m-1
        K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
    end
    L = K

    # spectfactLog = zeros(4,N_ckms)
    i = 0
    errK = errR = 1
    Err = zeros(0,2)
    while (errK > ϵ || errR > ϵ) && i <= 10^5
        hL = h*L; FL = F*L

        # Stopping criteria stuff
        i += 1
        FL_RrhLt = FL/Rr*hL'
        hL_RrhLt = hL/Rr*hL'
        errK = norm(FL_RrhLt)
        errR = norm(hL_RrhLt)
        Err = [Err; errK errR]
        i % 10 == 0 && println("err : $errK and $errR and i : $i" )


        K_new = K - FL_RrhLt
        L_new = FL - K/Re*hL
        Re_new = Re - hL_RrhLt
        Rr_new = Rr - hL'/Re*hL

        K = K_new
        L = L_new
        Re = Re_new
        Rr = Rr_new
    end

    println("i : $i")
    k = K/Re
    re = Re

    sqrt_re = sqrt(re)

    l = complex(zeros(d,d,m+1))
    l[:,:,1] = sqrt_re;
    for i = m-1:-1:0
        l[:,:,m-i+1] = k[d*i + 1: d*(i+1),:]*sqrt_re
    end

    # save("Data\\CKMS_dat.jld",
    #     "spectfactLog",
    #     spectfactLog)

    l, Err
end
