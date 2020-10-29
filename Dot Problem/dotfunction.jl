
"""
    dot(x, y)
    x ⋅ y
Compute the dot product between two vectors. For complex vectors, the first
vector is conjugated.
`dot` also works on arbitrary iterable objects, including arrays of any dimension,
as long as `dot` is defined on the elements.
`dot` is semantically equivalent to `sum(dot(vx,vy) for (vx,vy) in zip(x, y))`,
with the added restriction that the arguments must have equal lengths.
`x ⋅ y` (where `⋅` can be typed by tab-completing `\\cdot` in the REPL) is a synonym for
`dot(x, y)`.
# Examples
```jldoctest
julia> dot([1; 1], [2; 3])
5
julia> dot([im; im], [1; 1])
0 - 2im
julia> dot(1:5, 2:6)
70
julia> x = fill(2., (5,5));
julia> y = fill(3., (5,5));
julia> dot(x, y)
150.0
```
"""
function dot end

function dot(x, y) # arbitrary iterables
    ix = iterate(x)
    iy = iterate(y)
    if ix === nothing
        if iy !== nothing
            throw(DimensionMismatch("x and y are of different lengths!"))
        end
        return dot(zero(eltype(x)), zero(eltype(y)))
    end
    if iy === nothing
        throw(DimensionMismatch("x and y are of different lengths!"))
    end
    (vx, xs) = ix
    (vy, ys) = iy
    s = dot(vx, vy)
    while true
        ix = iterate(x, xs)
        iy = iterate(y, ys)
        ix === nothing && break
        iy === nothing && break
        (vx, xs), (vy, ys) = ix, iy
        s += dot(vx, vy)
    end
    if !(iy === nothing && ix === nothing)
        throw(DimensionMismatch("x and y are of different lengths!"))
    end
    return s
end

dot(x::Number, y::Number) = conj(x) * y

function dot(x::AbstractArray, y::AbstractArray)
    lx = length(x)
    if lx != length(y)
        throw(DimensionMismatch("first array has length $(lx) which does not match the length of the second, $(length(y))."))
    end
    if lx == 0
        return dot(zero(eltype(x)), zero(eltype(y)))
    end
    s = zero(dot(first(x), first(y)))
    for (Ix, Iy) in zip(eachindex(x), eachindex(y))
        @inbounds s += dot(x[Ix], y[Iy])
    end
    s
end

foo(n) = randn(n) + im*randn(n)

dot(foo(10^7),foo(10^7))
