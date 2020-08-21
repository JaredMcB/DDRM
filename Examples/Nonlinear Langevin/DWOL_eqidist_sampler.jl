using Distributions
using JLD


## Preprosesing for DWOL_dist_samp

function DWOL_dist_samp(N::Int64 = 1, σ = .35)
    # sampling
    μ = σ^2/2
    p(z) =  exp(-(x^2-1)^2/μ)
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
