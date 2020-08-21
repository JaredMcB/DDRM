using StatsBase



function Autocov(X,lags = 0:50)
    d, steps = size(X)
    A = zeros(d,d,length(lags))
    for i = 1:d
        for j = 1:d
            A[i,j,:] = crosscov(X[i,:],X[j,:],lags)
        end
    end
    A
end

function z_spect(X,L = 50; win = "Bart")
    lags = 0:L;
    A = Autocov(X,lags)
    A[:,:,1] = A[:,:,1]/2

    if win == "Bar"
        lam = 1 .- (0:L)/L
    elseif win == "Tuk"
        lam = .5*(1 .+ cos.(pi/L*(0:L)))
    elseif win == "Par"
        LL = Int(floor(L/2))
        lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
        lam2 = 2*(1 .- (LL+1:L)/L).^3
        lam = [lam1; lam2]
    else
        lam = ones(L+1)
    end

    z_spect_num(z) = sum([lam[i + 1]*(A[:,:,1+i]*z^(-i) + A[:,:,1+i]'*z^(i))  for i = 0 : L])
end



function Crosscov(X,Y,lags = 0:50)
    d, stepsx = size(X)
    nu, stepsy = size(Y)

    lags = unique([-reverse(lags) lags])

    if stepsx != stepsy
        print("X and Y are not the same length. Taking min.")
    end

    steps = minimum([stepsx stepsy])

    C = zeros(d,nu,length(lags))
    for i = 1:d
        for j = 1:nu
            C[i,j,:] = crosscov(X[i,1:steps],Y[j,1:steps],lags)
        end
    end
    C
end

function z_crossspect(X,Y,L = 50; win = "Bart")
    lags = 0:L;
    C = Cutocov(X,Y,lags)

    if win == "Bar"
        lam = 1 .- (0:L)/L
    elseif win == "Tuk"
        lam = .5*(1 .+ cos.(pi/L*(0:L)))
    elseif win == "Par"
        LL = Int(floor(L/2))
        lam1 = 1 .- 6*((0:LL)/L).^2 .+ 6*((0:LL)/L).^3
        lam2 = 2*(1 .- (LL+1:L)/L).^3
        lam = [lam1; lam2]
    else
        lam = ones(L+1)
    end

    z_crossspect_num(z) = sum([lam[abs(i) + 1]*C[:,:,L+1+i]*z^(-i) for i = -L : L])
end

function Matrix_CKMS_c(P)
    d = size(P)[1];
    m = size(P)[3] - 1

    NN = reverse(P[:,:,2:end],dims = 3)
    Re = Rr = p0 = P[:,:,1]

    F = [[zeros(d,d*(m-1)); I] zeros(d*m,d)]
    h = [zeros(d,d*(m-1)) I]

    K = zeros(d*m,d)
    for i = 0 : m-1
        K[d*i + 1: d*(i+1),:] = NN[:,:,i+1]
    end
    L = K

    for i = 1:200
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
    re = Re

    sqrt_re = sqrt(re)

    l = [sqrt_re]
    append!(l,[k[d*i + 1: d*(i+1),:]*sqrt_re for i = m-1:-1:0])
    l
end
