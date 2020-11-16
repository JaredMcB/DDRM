function modgen_L63(;
    sig = 10,
    bet = 8/3,
    rho = 28,
    steps = 10^4 +1,
    t_start = 0,
    t_stop = 100,
    discard = 10^4)

    steps_tot = steps + discard

    t = range(t_start, t_stop, length = steps)
    h = t[2] - t[1]
    x_init = randn(3)

    # Vector field
    F(x) = [sig*(x[2] - x[1]);
        x[1]*(rho - x[3]) - x[2];
        x[1]*x[2] - bet*x[3]]

    # Numerical solution (Euler Maruyama)
    X = zeros(3,steps_tot)
    X[:,1] = x_init
    for n = 1: steps_tot-1
        X[: , n+1] = X[:, n] .+ h*F(X[:,n]) .+ sqrt(h)*randn(3)
    end

    X = X[:,discard + 1 : steps_tot]
end
