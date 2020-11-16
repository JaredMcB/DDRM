## From WienerROM.jl

module KLPowerSpec

using FFTW,Statistics

powerspec(x; flags...) = powerspec(x,x;flags...)

function powerspec(x::Union{Vector{T},SubArray{T}},
                   y::Union{Vector{T},SubArray{T}};
                   blks = 100,
                   err     = false,
                   ) where {T}
    
    ## subdivide time series into blks blocks of nrows each
    nrows = div(length(x),blks)

    ## windowing function
    periw = periwindow(nrows)
    
    sxyl  = map(j->periodogram(periw(view(x,((j-1)*nrows+1):(j*nrows))),
                               periw(view(y,((j-1)*nrows+1):(j*nrows)))),
                1:blks)

    sxym = hcat(sxyl...)

    ## make sure we return vectors
    if err
        return (mean(sxym, dims=2)[1:end],
                ## 1 standard error
                sqrt.(var(sxym, dims=2)[1:end]))
    else
        return mean(sxym, dims=2)[1:end]
    end
end

function periodogram(x::Union{Vector{T},SubArray{T}},
                     y::Union{Vector{T},SubArray{T}}) where {T}
    
    ## This is the standard periodogram.  It causes aliasing
    ## in the dc mode, but this actually helps us get lower
    ## freq components right without needing to interpolate.
    if x === y
        ## Optimize a common case:
        return map(x->x*conj(x), fft(x))
    else
        ## This is organized so that when x===y, it is
        ## guaranteed to produce a positive answer.  But is
        ## the compiler smart enough to reduce this code to
        ## the above case when x===y?
        return fft(x) .* conj(fft(y))
    end
end

## Bartlett window function
function periwindow(N)
    wsum = 1
    w(j) = (1 - abs((j-N/2)/(N/2)))/wsum
    wsum = sqrt(sum(j->w(j)^2, 0:(N-1)))
    function(x)
        map(j->x[j+1] * w(j), 0:(N-1))
    end
end
end#module
