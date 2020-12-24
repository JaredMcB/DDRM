################################
## util.jl

## Note: would be better to split general utilities from
## mathematical ones.

using Polynomials


################################
## Take a structured numeric array (matrix, vector of
## vectors, etc) and "flatten" it into a numeric vector.  If
## possible, the new array is aliased to the old one, i.e.,
## changing an entry in the original object will mutate the
## corresponding entry in the flattened array.

## (This is not always possible to do, however, for example
## when flatting an array of arrays.)

function flatten(a::Matrix{Array{T,N}}) where {T,N}
    flatten(map(flatten,reshape(a,prod(size(a)))))
end

function flatten(a::Matrix{T}) where {T}
    reshape(a,prod(size(a)))
end

function flatten(v::Vector{Vector{T}}) where {T}
    if length(v) == 1
        return v[1]
    elseif isodd(length(v))
        return vcat(v[1],flatten(v[2:end]))
    else
        return vcat(flatten(v[1:div(end,2)]),
                    flatten(v[(div(end,2)+1):end]))
    end
end

function flatten(a::Array{T,N}) where {T,N}
    reshape(a,prod(size(a)))
end


################################
## more math utilities

function complex_randn()
    ## note variance normalized to 1
    complex(0.7071067811865476*randn(),
            0.7071067811865476*randn())
end

import Base.randn
function randn(::Type{Complex{Float64}})
    complex_randn()
end

function randn(::Type{Complex{Float64}}, n::Int)
    [complex_randn() for i=1:n]
end

function unif(a,b,args...)
    a+(b-a)*rand(args...)
end

function parity(n::Int)
    n&1
end

function length_mod2(a)
    parity(length(a))
end

function applyfilter(u,v)
    p = length(u)-1
    q = length(v)-1
    conv(u,v)[(1+min(q,p)):(1+max(q,p))]
end

function msqrt(A)
    U,d,Vt=svd(A)
    U*diagm(0 => sqrt.(d))*U'
end

function fftshift!(v)
    copyto!(v, fftshift(v))
end

function dist2(a::Number,b::Number)
    return abs2(a-b)
end

function dist2(A::Array{T,NT},
               B::Array{T,NT}) where {T,NT}
    sum = zero(T)
    foreach((a,b)->sum += dist2(a,b),A,B)
    return sum
end


################################
## BLAS & LAPACK wrappers

## Note BLAS and LAPACK generally expect the output variable
## to be the LAST argument, while Julia expects the output
## to be the FIRST.

## equivalent to sum += A*v
function linear_accum!(sum::AbstractVector{T},
                       A::AbstractMatrix{T},
                       v::AbstractVector{T}) where {T}
    e=one(T)
    BLAS.gemv!('N', e, A, v, e, sum)
end

## equivalent to sum += A*V
function linear_accum!(sum::AbstractMatrix{T},
                       A::AbstractMatrix{T},
                       v::AbstractMatrix{T}) where {T}
    e=one(T)
    BLAS.gemm!('N', 'N', e, A, v, e, sum)
end

## equivalent to sum += a*x
function linear_accum!(sum::AbstractVector{T},
                       a::Number, x::AbstractVector{T}) where {T}
    BLAS.axpy!(a, x, sum)
end
    
## equivalent to sum += a*x
function linear_accum!(sum::Matrix{T},
                       a::Number, x::Matrix{T}) where {T}
    
    BLAS.axpy!(a, flatten(x), flatten(sum))
end

## slower version
function linear_accum!(sum::AbstractMatrix{T},
                       a::Number, x::AbstractMatrix{T}) where {T}
    m,n = size(sum)
    for i=1:m
        for j=1:n
            sum[i,j] += a*x[i,j]
        end
    end
end

## equivalent to sum += a*x + b*y
function linear_accum!(sum, a, x, b, y)
    linear_accum!(sum, a, x)
    linear_accum!(sum, b, y)
end
    
function zero!(x)
    fill!(x,zero(x[1]))
end


## Solve linear least squares by calling LAPACK.gels!().
## Since gels!() is destructive, this wrapper lets us
## allocate memory just once and takes care of the copying.
## This is better than calling \, which reallocates memory
## every time.

make_lsqsolver(m,n) = make_lsqsolver(Float64,m,n)

## This version calls LAPACK, which (presumably) uses SVD.
function make_lsqsolver(T::DataType,m::Int,n::Int)
    let Acc = zeros(T,m,n),
        bcc = zeros(T,m)

        warn("make_lsqsolver(): using SVD")

        function solve(A,b)
            copyto!(Acc,A)
            copyto!(bcc,b)
            fact,x,sqerrs=LAPACK.gels!('N', Acc, bcc)
            return x,real(sum(sqerrs))
        end
    end
end

## This version solves the normal equations using LAPACK's
## positive-definite linear solver.
# function make_lsqsolver(T::DataType,m::Int,n::Int)
#     let Acc     = zeros(T,n,n),
#         bcc     = zeros(T,n),
#         resid   = zeros(T,m),
#         theone  = one(T),
#         thezero = zero(T)

#         warn("make_lsqsolver(): solving normal eqs")

#         function solve(A,b)
#             BLAS.gemm!('C','N',theone,A,A,thezero,Acc)
#             BLAS.gemv!('C',theone,A,b,thezero,bcc)
#             LAPACK.posv!('U', Acc, bcc)

#             let x = bcc
#                 copyto!(resid,b)
#                 BLAS.gemv!('N',theone,A,x,-theone,resid)
#                 return x,sum(abs2,resid)
#             end
#         end
#     end
# end


export print_to_file, print_to_mfile, timed_for, foreach, make_timer

##------------------------------------------------------

## this does not work
## export with_output_to_file

## function with_output_to_file(file, thunk)
##     let curr_stdout=stdout
##         let fp=open(file, "w")
##             redirect_stdout(fp)
##             thunk()
##             close(fp)
##         end
##         redirect_stdout(curr_stdout)
##     end
## end

##------------------------------------------------------

## note: these have been made obsolete by the MAT package

function print_to_file(file, x; append=true)
    if append
        fp = open(file, "a")
    else
        fp = open(file, "w")
    end
    print(fp,x)
    close(fp)
end

function print_to_mfile(file, x, xname; append=false, force=false)
    if isfile(file)
        if append
            fp = open(file, "a")
        elseif force
            fp = open(file, "w")
        else
            println("#print_to_mfile(): ", file, " exists")
            println("#print_to_mfile(): to proceed, set append or force to true")
            error()
        end
    else
        if append
            fp = open(file, "a")
        else
            fp = open(file, "w")
        end
    end
        
    # make a special rule for complex things
    if iseltype(x,Complex)
        print(fp, xname, "=[")
        print(fp,real(x))
        print(fp, "]+sqrt(-1)*[")
        print(fp,imag(x))
        println(fp, "];")
    else
        print(fp, xname, "=[")
        print(fp,x)
        println(fp, "];")
    end
    close(fp)
end

######################################################
## More, from lib/util.jl

import Base.foreach
function foreach(f::Function, range, tag::AbstractString; delay=1.0)
    let t0::Float64 = time(),
        tlp = 0.0,
        count::Int64 = 0,
        total::Int64 = length(range)

        let imin = minimum(range),
            imax = maximum(range)

            for i in range
                f(i)
                count = count + 1
                dt = time() - t0
                if dt >= tlp + delay || i == imin || i == imax
                    tlp = dt
                    println("#", tag, ": ", count, "/",
                            total,
                            " steps took ",
                            nicedate(dt), "; eta ",
                            nicedate((dt/count)*(total-count)))
                    flush(stdout)
                end
            end
        end
        return time()-t0
    end
end

function foreach(f::Function, range)
    t0::Float64 = time()
    for i in range
        f(i)
    end
    return time()-t0
end

##------------------------------------------------------

function timed_for(f::Function, range; delay=1.0)
    timed_for(f, range, "timed_for", delay=delay)
end

function timed_for(f::Function, range, tag; delay=1.0)
    foreach(f, range, tag, delay=delay)
end

##------------------------------------------------------

function make_timer(f::Function; delay=5.0)
    let t0    = time(),
        tlast = t0
        function()
            let tnow=time()
                if tnow-tlast >= delay
                    f(tnow-t0)
                    tlast = tnow
                end
                return tnow-t0
            end
        end
    end
end

function make_timer()
    make_timer(t->nothing)
end

##------------------------------------------------------

function nicedate(sec::Float64)
    if sec < 60
        sec = chop(sec)
        return "$sec sec"
    else
        min = int_floor(floor(sec / 60))
        sec = chop( sec - 60*min )
        if min < 60
            return "$min min $sec sec"
        else
            hrs = int_floor(floor(min / 60))
            min -= 60*hrs
            if hrs < 24
                return "$hrs hrs $min min $sec sec"
            else
                days = int_floor(floor(hrs / 24))
                hrs -= 24*days
                return "$days days $hrs hrs $min min $sec sec"
            end
        end
    end
end

function chop(t::Float64; n::Int64=3)
    int_floor(t*10^n)/10^n
end

function int_floor(x::Float64)
    Integer(floor(x))
end

##------------------------------------------------------

## Factor a real polynomial into a product of real quadratic
## factors, times one real linear factor if the degree is
## odd.  Issue: we don't know if Polynomial.roots() returns
## conjugate roots that are numerically exact.  So we trust,
## but verify.

quadfact(A::Vector; flags...) = quadfact(Polynomial(A); flags...)

function quadfact(A::Polynomial; tol=0.0)
    let p    = degree(A),

        rts  = roots(A),
        rrts = filter(z->imag(z)==0, rts),

        urts = filter(z->imag(z) > 0, rts)
        sort!(urts; by=z->abs(imag(z)), alg=InsertionSort)
        sort!(urts; by=z->real(z),      alg=InsertionSort)
        
        lrts = filter(z->imag(z) < 0, rts)
        sort!(lrts; by=z->abs(imag(z)), alg=InsertionSort)
        sort!(lrts; by=z->real(z),      alg=InsertionSort)

        if length(urts) != length(lrts)
            println("# urts: ", urts)
            println("# lrts: ", lrts)
            error("#urts != #lrts")
        end

        qfl = Polynomial[]

        for k=1:length(urts)
            let P = Polynomial([-lrts[k],1])*Polynomial([-urts[k],1])
                @assert maximum(imag.(P[0:end])) <= tol
                push!(qfl, Polynomial(real.(P[0:end])))
            end
        end
        
        let k=1
            while k <= length(rrts)
                if k+1 <= length(rrts)
                    push!(qfl, Polynomial([-real.(rrts[k]),1])*Polynomial([-real.(rrts[k+1]),1]))
                    k += 2
                else
                    push!(qfl, Polynomial([-real.(rrts[k]),1]))
                    k += 1
                end
            end
        end
        
        return sort(qfl,by=degree)
    end
end


##------------------------------------------------------
## Handy utilities for working with expressions

function parseval(s)
    println("$(s)")
    eval(Meta.parse(s))
end

function once_id(expr)
    Symbol("once_id$(hash(expr))")
end

macro once(expr)
    let id = once_id(expr)
        return :(
                 $(esc(id)) = if isdefined(Main,$(Expr(:quote,id)))
                          $(esc(id))
                       else
                          $(esc(expr))
                       end
                 )
    end
end

################################
## memoize expensive functions
function memoize(f)
    h = Dict{Any,Any}()
    function fmem(args...;flags...)
        key = (args,flags)
        if !haskey(h,key)
            h[key] = f(args...;flags...)
        end
        return h[key]
    end
end


##------------------------------------------------------
## For backward compatiblity
echo(x...) = println("# ", x...)
