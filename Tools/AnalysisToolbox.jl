using DSP, Statistics

function autocov_con(x::AbstractVector{<:Number},lags::UnitRange{Int})
    lx = size(x,1)
    if maximum(lags) > lx
        print("lag cannot be greater than lenght of series")
        lags = filter(x -> abs(x) < lx, lags)
    end

    x .-= mean(x)
    A = conv(x,conj(reverse(x)))/lx
    A = [A[k + lx] for k in lags]
end

function autocor_con(x::AbstractVector{<:Number},lags::UnitRange{Int})
    A = autocov_con(x,lags)
    A ./= autocov_con(x,0:0)[1]
end

function auto_times(x::AbstractVector{<:Real};plt = false)
    lx = size(x,1)
    L = minimum([lx, 10^6])

    lags = 0:L
    A = autocov_con(x,lags)

    end_int = try
                findall(A.<0)[1] - 1
            catch e
                if isa(e, BoundsError)
                    L
                end
            end
    end_exp = Int64(round(end_int/3))

    A_mat = [ones(end_exp,1) reshape(1:end_exp,end_exp,1)]
    b = inv(A_mat'*A_mat)*A_mat'*log.(A[1:end_exp])
    τ_exp = -1/b[2]

    A ./= A[1]
    τ_int = .5 + sum(A[2:end_int])
    if plt
        P = plot(0:(end_int-1),log.(A[1:end_int]),
            ylabel = "Log of Autocov",
            xlabel = "Lags",
            label = "log(A)")
        P = plot!(P,0:end_exp-1,A_mat*b,
            label = "linear approx of log(A)")
        end
    plt ? [τ_exp, τ_int, P] : [τ_exp, τ_int]
end
