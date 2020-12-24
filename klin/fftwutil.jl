##----------------------------------------------------------------------------
module fftwutil

Complex128 = Complex{Float64}

using FFTW

export best_dim, best_odd_dim, best_even_dim, best_congruent_dim, realize, realize!, fourier

## ## compute current directory
## DIRPATH = dirname(Base.source_path())
## println("#fftwutil: DIRPATH = $DIRPATH")

## ## DIRPATH = string(ENV["HOME"], "/statmech/nmda/models")
## const LIBPATH = "$DIRPATH/fftwutil.dylib"

## run(`make --quiet -C $DIRPATH`)

## ## Note on @eval: ccall(), for whatever reason, expects its
## ## first argument to be a constant.

## function best_dim(n::Int64)
##     ccall((:fftw_best_dim, LIBPATH), Int64, (Int64,), n)
## end

## function best_odd_dim(n::Int64)
##     ccall((:fftw_best_odd_dim, LIBPATH), Int64, (Int64,), n)
## end

## function best_even_dim(n::Int64)
##     ccall((:fftw_best_even_dim, LIBPATH), Int64, (Int64,), n)
## end

## function best_congruent_dim(M::Int64, a::Int64, modulus::Int64)
##     ccall((:fftw_best_congruent_dim, LIBPATH),
##           Int64, (Int64,Int64,Int64),
##           M, a, modulus)
## end

##--------------------------------
## Julia versions

best_dim_table = [1]

function best_dim(M;
                  alist = [2, 3, 5, 7],
                  #blist = [11,13],
                  blist = [],  ## seems FFTW defaults have changed
                  )
    
    if M > length(best_dim_table)
        global best_dim_table = fill(0,2^(1+floor(Int,log2(M))))
        fill_dim_table!(best_dim_table, alist, blist)
    end
    
    return best_dim_table[M];
end

## This is probably the most efficient possible algorithm.
## For the meaning of the 'alist' and the 'blist', see FFTW
## documentation.

function fill_dim_table!(table, alist, blist)

    ## clear table first:
    fill!(table, 0)

    ## 1 is always in the list:
    table[1] = 1

    ## then take care of the alist:
    for a in alist
        for j=1:length(table)
            if table[j] > 0
                p = j*a
                if p <= length(table)
                    table[p] = 1
                else
                    break
                end
            end
        end
    end

    ## take care of the blist:
    for b in blist
        for j=1:length(table)
            if table[j] == 1
                p = j*b
                if p <= length(table)
                    table[p] = 2
                else
                    break
                end
            end
        end
    end

    ## convert the table into a useful form:
    let n=0
        for i=length(table):-1:1
            if table[i] > 0
                n = i
            end
            table[i] = n
        end
    end
end




## Sometimes it's useful to get an odd integer.  This is not
## the best algorithm, but it's the shortest.

function best_odd_dim(M)
    p = best_dim(M)
    if isodd(p)
        return p
    else
        return best_odd_dim(p+1)
    end
end


## Sometimes it's useful to get an even integer.

function best_even_dim(M)
    p = best_dim(M)
    if iseven(p)
        return p
    else
        return best_even_dim(p+1)
    end
end


## A slight generalization that is a little dangerous: it
## doesn't always terminate...

function best_congruent_dim(M, a, modulus;
                            maxiter = 1024)
    
    ## sanity check
    if  modulus == 0
        return M
    elseif a < 0 || a >= modulus 
        a %= modulus
    end
    
    ## don't let this run away
    n = best_dim(M)
    
    for count = 1:maxiter
        if n % modulus == a 
            return n
        else
            n = best_dim(n+1)
        end
    end

    @warn "best_congruent_dim() failed to terminate"
    return M
end



##--------------------------------
## given a vector of fourier coefficients, compute the
## corresponding real-valued function

function fourier(u)
    N = length(u)
    if isodd(N)
        fft(u)[2:div(N+1,2)]/N
    else
        fft(u)[2:div(N,2)]/N
    end
end

function realize(u::Array{Complex128,1};
                 NN = best_congruent_dim(2*length(u)+1,1,2))

    v = complex(zeros(NN))
    realize!(u,v)
    real(v)
end

function realize!(u::Array{Complex128,1},
                  v::Array{Complex128,1})
    N  = length(u)
    NN = length(v)

    v[1] = 0.0
    
    let i=2, ii=NN
        for k=1:N
            v[i] = u[k]
            v[ii] = conj(u[k])
            i = i+1
            ii = ii-1
        end

        while i <= ii
            v[i] = v[ii] = 0
            i = i+1
            ii = ii-1
        end
    end

    bfft!(v)
end

## 2d, spacetime version, assuming 2nd index is time
function realize(u::Array{Complex128,2};
                 NN = best_congruent_dim(2*size(u)[1]+1, 1, 2))

    N,M = size(u)
    v = zeros(Complex128,NN)
    w = zeros(NN,M)

    for i=1:M
        realize!(u[:,i],v)
        for j=1:NN
            w[j,i] = real(v[j])
        end
    end
    return w
end

##--------------------------------
## This is for backward compatibility.

export make_fft!
function make_fft!(args...; flags...)
    p = plan_fft!(args...; flags...)
    function doit!(A)
        p * A
    end
end

export make_bfft!
function make_bfft!(args...; flags...)
    p = plan_bfft!(args...; flags...)
    function doit!(A)
        p * A
    end
end

export make_r2r!
function make_r2r!(args...; flags...)
    p = FFTW.plan_r2r!(args...; flags...)
    function doit!(A)
        p * A
    end
end

end # module
