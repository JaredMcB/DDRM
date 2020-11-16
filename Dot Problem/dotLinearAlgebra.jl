using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
foo(n) = randn(n) + im*randn(n)

dot(foo(10001),foo(10001))

for n = 3.9:0.01:5
    dot(foo(floor(Int,10^n)),foo(floor(Int,10^n)))
    println(n)
end
