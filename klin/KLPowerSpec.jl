##----------------------------------------------------------------------------
## tcorr.jl

## Compute time correlations and power spectra for time
## series.

@assert VERSION >= v"0.6"

module tcorr

using DSP, FFTW, Statistics

export time_corr,powerspec


##------------------------------------------------------
## Method for scalar time series, real or complex.

#=

TODO

1) The correlation function computed here is multiplied by a
"triangular" filter; this convolves the power spectrum
with the square of the sinc function.

2) This may be slow because it computes a big transform and
discards most of it, and reallocates lots of temporary
storage every time.

=#

function time_corr(x::AbstractVector{T}, ncorr::Int64=0; flags...) where {T}
    time_corr(x, x, ncorr; flags...)
end

function time_corr(x::AbstractVector{T},
                   y::AbstractVector{T},
                   ncorr::Int64=0;
                   shift    = true,
                   zeromean = true,
                   ) where {T}

    let y     = copy(y),
        nx    = length(x),
        ny    = length(y),
        ny_2  = div(ny,2),
        clen  = max(nx,ny)-min(nx,ny)+1+2*(min(nx,ny)-1)
        n0    = isodd(clen) ? div(clen,2)+1 : div(clen,2)

        if ncorr == 0
            ncorr = n0-1
        end

        ## Reverse and conjugate y in one loop, to reduce
        ## copying and loop overhead.
        for i=1:ny_2
            let tmp = y[i]
                y[i]      = conj(y[ny+1-i])
                y[ny+1-i] = conj(tmp)
            end
        end

        ## I really don't know how important this is...
        if isodd(ny)
            y[ny_2+1] = conj(y[ny_2+1])
        end

        ## This normalization implicitly applies a
        ## triangular windowing function, commonly known as
        ## Bartlett windowing.
        let mx   = mean(x),
            my   = mean(y),
            corr = (zeromean ? conv(x .- mx,(y .- my)) : conv(x,y)) ./ min(nx,ny) 

            ## sanity check
            @assert length(corr) == clen

            if isodd(clen)
                let corr = view(corr,(n0-ncorr):(n0+ncorr))
                    if shift
                        return mx,my,ifftshift(corr)
                    else
                        return mx,my,corr
                    end
                end
            else
                let corr = view(corr,(n0-ncorr):(n0+ncorr+1))
                    if shift
                        return mx,my,ifftshift(corr)
                    else
                        return mx,my,corr
                    end
                end
            end
        end
    end
end


##------------------------------------------------------
## Power spectra via periodograms.

## Subdivide the time series into =blks= segments and
## estimating a power spectrum from each, then average.
## This lets us estimate variance.

## TODO

#=

1) Avoid FFTs: currently powerspec() does a bfft() on
time_corr(), which itself calls conv(), which presumably
does 3 FFTs.  This can be sped up 2x easily.

2) Avoid extra copying.  'Nuf said.

# =#


## This is the original version.  It's slower than the NR
## version, and does less smoothing.

kpowerspec(x; flags...) = kpowerspec(x, x; flags...)

function kpowerspec(x::Union{Vector{T},SubArray{T}},
                    y::Union{Vector{T},SubArray{T}};
                    blks = 100,
                    err     = false,
                    ) where {T}
    
    nrows = div(length(x),blks)

    if x === y
        sxyl = map(j->kperiodogram(x[((j-1)*nrows+1):(j*nrows)]),
                   1:blks)
    else
        sxyl = map(j->kperiodogram(x[((j-1)*nrows+1):(j*nrows)],
                                   y[((j-1)*nrows+1):(j*nrows)]),
                   1:blks)
    end
    
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

function kperiodogram(x::Union{Vector{T},SubArray{T}}) where {T}
    spec = bfft(time_corr(x;zeromean=false)[3])
    map!(real, spec, spec)
end

function kperiodogram(x::Union{Vector{T},SubArray{T}},
                      y::Union{Vector{T},SubArray{T}}) where {T}
    bfft(time_corr(x,y;zeromean=false)[3])
end

## Version from Numerical Recipes

powerspec(x; flags...) = powerspec(x,x;flags...)

function powerspec(x::Union{Vector{T},SubArray{T}},
                   y::Union{Vector{T},SubArray{T}};
                   blks = 100,
                   err     = false,
                   ) where {T}
    
    ## subdivide time series into blks blocks of nrows each
    nrows = div(length(x),blks)

    ## windowing function
    window = PeriWindow(nrows)
    
    sxyl  = map(j->periodogram(window(view(x,((j-1)*nrows+1):(j*nrows))),
                               window(view(y,((j-1)*nrows+1):(j*nrows)))),
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
        return map(x->complex(abs2(x)), fft(x))
    else
        ## This is organized so that when x===y, it is
        ## guaranteed to produce a positive answer.  But is
        ## the compiler smart enough to reduce this code to
        ## the above case when x===y?
        return fft(x) .* conj(fft(y))
    end
end

## Bartlett window function
function PeriWindow(N)
    wsum = 1
    w(j) = (1 - abs((j-N/2)/(N/2)))/wsum
    ## normalize L^2 norm
    wsum = sqrt(sum(j->w(j)^2, 0:(N-1)))
    function(x)
        map(j->x[j+1] * w(j), 0:(N-1))
    end
end
end # module
