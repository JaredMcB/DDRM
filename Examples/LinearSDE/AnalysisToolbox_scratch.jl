x = X1[:]

lx = size(x,1)
L = minimum([lx - 1, 10^6])

lags = 0:L
A = my_autocov(x,lags)
y = X1[:]

lx = size(x,1)
ly = size(y,1)
lx == ly || throw(DimensionMismatch("series must be same length"))

if maximum(lags) > lx
    println("lag cannot be greater than lenght of series")
    lags = filter(x -> abs(x) < lx, lags)
end

x .-= mean(x)
y .-= mean(y)
C = conv(x,conj(reverse(y)))/lx
C = [C[k + lx] for k in lags]









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
