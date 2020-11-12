
module JMDataGen

using Distributions

function DataGen_DWOL(;
    #SDE parameters
    sigma    = [1],
    V_prime  = x -> -x.*(x.^2 .- 1),
    sig_init = [1.5],
    # Numerical estimate parameters
    scheme   = "FE",
    steps    = 10^6, # Nomber of time steps (not including those discarded)
    h        = .01,
    discard  = steps, # Number of time steps discarded
    gap      = 100, #1 + the number of time steps between observations
    ObsNoise = false, # if true keeps track of an returns noise e
    e        = ObsNoise ? randn(size(sigma,1),steps + discard) : 0
    )

    d = size(sigma,1)

    if discard == 0
        sig_init = [DWOL_dist_samp(1,σ = sigma)]
    end

    size(sig_init,1) == size(sigma,1) || print("Dimension of initial value "*
                                               "and sigma do not agree")

    steps_tot = steps + discard

    # So the multiplication is defined when d=1
    d == 1 && (sigma = reshape(sigma,1,1))

    # Here we genereate the signal process.
    signal = zeros(d,ceil(Int,steps/gap))
    tmp = sig_init
    if scheme == "FE"
        for n = 1 : steps_tot-1
            tmp = tmp + h*V_prime(tmp) +
                        sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "T2"
        for n = 1 : steps_tot-1
            tmp = tmp .+ h*(-tmp.^3 .+ tmp) .+
                            h^2/2*( tmp.*(tmp.^2 .- 1).*(tmp.^2 .- 1) ) .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "EM"
        for n = 1 : steps_tot-1
            k1 = tmp .+ h/2*V_prime(tmp)
            tmp = tmp .+ h*V_prime(k1) .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "ET"
        for n = 1 : steps_tot-1
            k0 = V_prime(tmp)
            k1 = tmp .+ h*k0
            tmp = tmp .+ h/2*V_prime(k1) .+ h/2*K0 .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    elseif scheme == "RK4"
        for n = 1 : steps_tot-1
            k1 = tmp
            k2 = tmp .+ h/2*V_prime(k1)
            k3 = tmp .+ h/2*V_prime(k2)
            k4 = tmp .+ h*V_prime(k3)
            tmp = tmp .+ h/6*( V_prime(k1) .+ 2*V_prime(k1) .+
                            2*V_prime(k1) .+ V_prime(k1) ) .+
                            sqrt(h)*sigma*( ObsNoise ? e[:,n+1] : randn(d) )
            if (n - discard > 0) && ((n - discard - 1) % gap == 0)
                signal[:,(n-discard-1)÷gap + 1] = tmp
            end
        end
    end

    result = ObsNoise ? [signal, e] : signal
end

function DWOL_dist_samp(
    N::Int64 = 1; # N is the number of samples you want
    σ = [.35] # The σ for the noise of the the DWOL
    )
    # sampling
    μ_pos = 1
    μ_neg = -1
    v = 1/(sqrt(2π)*1.05)

    D_neg = Normal(μ_pos,v)
    D_pos = Normal(μ_neg,v)

    q(x) = exp(-(x-μ_pos)^2/(2v^2))/(2v*sqrt(2\pi)) +
            exp(-(x-μ_neg)^2/(2v^2))/(2v*sqrt(2\pi))

    μ = σ[1]^2/2
    p(x) =  exp(-(x^2-1)^2/μ)
    c = 1.05*p(0)/q(0)
    cq(x) = c*q(x)
    Z = zeros(N)
    for n = 1:N
        z = (rand() < .5 ? rand(D_neg) : rand(D_pos))
        while rand()*cq(z) > p(z)
            z = (rand() < .5 ? rand(D_neg) : rand(D_pos))
        end
        Z[n] = z
    end
    N == 1 ? Z[1] : Z
end

end#module
